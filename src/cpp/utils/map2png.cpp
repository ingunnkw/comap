#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <spng.h>
#include <scolorizer.h>
#include <simplefont.h>
#include <limits>
#include <rotmatrix.h>
#include <pointing.h>
#include <stimer.h>
#include <stypes.h>
#include <hmap_io.h>
#include <healpix_map.h>
#include <cstring>

using namespace std;
using namespace skn;

//static const double pi=3.141592653589793238462643383279502884197;
//static const double degr2rad=pi/180.0;

typedef Grid<double,2> Projection;
enum { PRO_MOLL, PRO_GNOM };
enum { TRF_LIN,  TRF_HIST, TRF_LOG };

typedef HMap<double> Map;
string set_file_extension(const string & s, const string & ext);

Projection pro_mollw(const Map & map, int sub, int sig, int xsize, double lon0, double lat0)
{
	Healpix_Base hpix(map.get_nside(), map.get_ordering() == 1 ? RING : NEST, nside_dummy());
	int ysize=xsize/2;
	Projection img(xsize,ysize,Healpix_undef);
	double xc=(xsize-1)/2., yc=(ysize-1)/2.;
	double lon0rad = lon0*degr2rad;
	double lat0rad = lat0*degr2rad;

	rotmatrix rot;
	rot.Make_CPAC_Euler_Matrix(0,-lat0rad,-lon0rad);

	for (int i=0; i<img.size(0); ++i)
	for (int j=0; j<img.size(1); ++j)
	{
		double u = 2*(i-xc)/(xc/1.02);
		double v =  -(j-yc)/(yc/1.02);
		bool mask = ((u*u/4 + v*v) <= 1);
		if (mask)
		{
			pointing ptg (pi/2-(asin(2/pi*(asin(v) + v*sqrt((1-v)*(1+v))))),
				-pi/2*u/max(sqrt((1-v)*(1+v)),1e-6));
			vec3 pnt = rot.Transform(ptg.to_vec3());
			img(i,j) = map.pix(sub,sig,hpix.ang2pix(pnt));
		}
	}
	return img;
}

Projection pro_gno(const Map & map, int sub, int sig, int xsize, double lon0, double lat0, double res)
{
	Healpix_Base hpix(map.get_nside(), map.get_ordering() == 1 ? RING : NEST, nside_dummy());
	double lon0rad = lon0*degr2rad;
	double lat0rad = lat0*degr2rad;

	rotmatrix rot;
	rot.Make_CPAC_Euler_Matrix(lon0rad,-lat0rad,0);

	double delta=res*degr2rad/60.;
	double start=-(xsize/2.)*delta;
	Projection img(xsize, xsize, Healpix_undef);
	for (int i=0; i<img.size(0); ++i)
	for (int j=0; j<img.size(1); ++j)
	{
		vec3 point (1,-(start+i*delta), -(start+j*delta));
		point = rot.Transform(point);
		pointing pnt = point;
		int pix = hpix.ang2pix(pnt);
		img(i,j) = map.pix(sub,sig,pix);
	}
	return img;
}

bool is_nan(double v) { return v != v; }
bool is_inf(double v) { return v == numeric_limits<double>::infinity() || v == -numeric_limits<double>::infinity(); }
bool is_hpbad(double v) { return approx(v,Healpix_undef); }
bool ok_val(double v) { return !is_nan(v) && !is_inf(v) && !is_hpbad(v); }

Vector<double,2> get_maxmin(const Projection & p)
{
	Vector<double,2> range(0,1);
	int i; for(i = 0; i < p.size(); i++) if(ok_val(p[i])) break;
	if(i >= p.size()) return range;
	range = p[i];
	for(++i; i < p.size(); i++)
		if(ok_val(p[i]))
		{
			if(p[i] < range(0)) range(0) = p[i];
			if(p[i] > range(1)) range(1) = p[i];
		}
	return range;
}

/* Transforms normal valued back and forwards. Do not call for values
 * that aren't ok! */
class Transformation
{
public:
	virtual ~Transformation() {}
	virtual double forwards(double v) const = 0;
	virtual double backwards(double v) const = 0;
};

class LinTrans : public Transformation
{
public:
	LinTrans(double from, double to):range(from, to) {}
	virtual double forwards(double v) const
	{
		return v <= range(0) ? 0 : v >= range(1) ? 1 : range(1) == range(0) ? 0 : (v-range(0))/(range(1)-range(0));
	}
	virtual double backwards(double v) const
	{
		return v*(range(1)-range(0)) + range(0);
	}
protected:
	Vector<double,2> range;
};

class HistTrans : public Transformation
{
public:
	HistTrans(const vector<double> & values):vals(values) {}
	virtual double forwards(double v) const
	{
		if(vals.size() <= 1) return 0;
		vector<double>::const_iterator i = lower_bound(vals.begin(), vals.end(), v);
		if(i == vals.end()) return 1;
		else return double(i-vals.begin())/(vals.size()-1);
	}
	virtual double backwards(double v) const
	{
		if(vals.empty()) return v;
		else if(v < 0) return vals[0];
		else if(v >= 1) return vals.back();
		else return vals[int(v*vals.size())];
	}
protected:
	vector<double> vals;
};


Image make_num(int margin, const SimpleFont & font, int dig, double v)
{
	const int buflen = 0x400;
	char format[buflen];
	if(abs(v) >= 100 || abs(v) <= 0.01)
		snprintf(format, buflen, "%%%d.%de", dig, dig-6);
	else
		snprintf(format, buflen, "%%%d.%df", dig, dig-4);

	char buf[buflen];
	snprintf(buf, buflen, format, v);
	int len = strlen(buf);
	int wx = len*font.xpix;
	Image res(2*margin + wx, 2*margin + font.ypix, Pixel(255,255,255));
	for(int i = 0; i < len; i++)
	for(int x = 0; x < font.xpix; x++)
	for(int y = 0; y < font.ypix; y++)
	{
		byte v = 255*(1-font.data[(buf[i]*font.ypix+y)*font.xpix + x]);
		res(i*font.xpix + x+margin, y+margin) = Pixel(v,v,v);
	}
	return res;
}

// Copy image from into image to, starting at position x,y in to.
void copy(const Image & from, Image & to, int x, int y)
{
	int nx = min(to.size(0)-x, from.size(0)),
		ny = min(to.size(1)-y, from.size(1));
	int i0 = x < 0 ? -x : 0, j0 = y < 0 ? -y : 0;
	for(int i = i0; i < nx; i++)
	for(int j = j0; j < ny; j++)
		to(x+i,y+j) = from(i,j);
}

Image make_colorbar(Transformation * trf, const Colorizer & color,
	const SimpleFont & font, int dig, int nx, int ticks,
	const Vector<double,2> & drange, const Vector<double,2> & urange,
	const Vector<bool,2> & ugiven)
{
	int margin = 5;
	Image dleft = make_num(margin, font, dig, drange(0)),
		dright = make_num(margin, font, dig, drange(1));
	Image uleft = make_num(margin, font, dig, urange(0)),
		uright = make_num(margin, font, dig, urange(1));

	int ny = dleft.size(1);
	Image res(nx, ny, Pixel(255,255,255));
	int xl = 0, xr = nx;
	copy(dleft, res, xl, 0); xl += dleft.size(0);
	if(ugiven(0)) { copy(uleft,  res, xl, 0); xl += uleft.size(0); }
	copy(dright, res, xr-dright.size(0)-1, 0); xr -= dright.size(0);
	if(ugiven(1)) { copy(uright, res, xr-uright.size(0)-1, 0); xr -= uright.size(0); }

	double a = trf->backwards(0), b = trf->backwards(1);

	// Find number of thicks
	int tickpow = a == b ? 0 : int(floor(log(b-a)/log(10.))-ticks+1);
	int maxtick = min(ticks, ny/3);

	int n = xr-xl;
	double dv = b == a ? 0 : (b-a)/(n-1);
	for(int i = 0, x = xl; i < n; i++, x++)
	{
		double v = dv*i + a, vn = v+dv;
		int ticksize = 0;
		if(ticks) while(ticksize < maxtick && floor(v/pow(10.,tickpow+ticksize)) != floor(vn/pow(10.,tickpow+ticksize))) ticksize++;

		double w = trf->forwards(v);
		Pixel p = color(w);
		for(int y = 0; y < ny; y++)
			res(x,y) = p;
		for(int y = ny-3*ticksize; y < ny; y++) res(x,y) = Pixel(0,0,0);
	}
	return res;
}

double nearest_whole(double a, double b)
{
	return fmod(a,b) < b/2 ? floor(a/b)*b : (floor(a/b)+1)*b;
}

Pixel blend(const Pixel & bg, const Pixel & fg, double alpha)
{
	Pixel res;
	for(int i = 0; i < 3; i++) res(i) = max(0,min(0xff,int(
		alpha*fg(i)+(1-alpha)*bg(i))));
	return res;
}

double grid_val(const pointing & pnt, double lat_interval, double lon_interval, double linewidth)
{
	double lon_nearest = nearest_whole(pnt.phi, lon_interval);
	double lat_nearest = nearest_whole(pnt.theta, lat_interval);
	double a = pnt.phi-lon_nearest, b = pnt.theta-lat_nearest;
	if(abs(a/linewidth) > 5 && abs(b/linewidth) > 5) return 0;
	double v1 = exp(-a*a/linewidth/linewidth), v2 = exp(-b*b/linewidth/linewidth);
	return min(1.0,2*(v1+v2));
}

void grid_gno (Image & img, double lon0, double lat0, double res, double lonint, double latint, double linewidth = 1)
{
	double lon0rad = lon0*degr2rad;
	double lat0rad = lat0*degr2rad;
	double lon_interval = lonint*degr2rad;
	double lat_interval = latint*degr2rad;

	rotmatrix rot;
	rot.Make_CPAC_Euler_Matrix(lon0rad,-lat0rad,0);

	double delta=res*degr2rad/60.;

	int xsize = img.size(0), ysize = img.size(1);
	double start=-(xsize/2.)*delta;
	for (int i=0; i < xsize; ++i)
	for (int j=0; j < ysize; ++j)
	{
		vec3 point (1,-(start+i*delta), -(start+j*delta));
		double v = grid_val(rot.Transform(point), lat_interval, lon_interval, linewidth*delta);
		if(v) img(i,j) = blend(img(i,j), Pixel(0,0,0), v);
	}
}

void grid_moll(Image & img, double lon0, double lat0, double lonint, double latint, double linewidth = 1)
{
	double lon0rad = lon0*degr2rad;
	double lat0rad = lat0*degr2rad;
	double lon_interval = lonint*degr2rad;
	double lat_interval = latint*degr2rad;

	rotmatrix rot;
	rot.Make_CPAC_Euler_Matrix(0,-lat0rad, -lon0rad);
	int xsize = img.size(0), ysize = img.size(1);
	double xc=(xsize-1)/2., yc=(ysize-1)/2.;
	double delta = (2*pi)/xsize;
	for (int i=0; i < xsize; ++i)
	for (int j=0; j < ysize; ++j)
	{
		double u = 2*(i-xc)/(xc/1.02);
		double v =  -(j-yc)/(yc/1.02);
		bool mask = ((u*u/4 + v*v) <= 1);
		if (mask)
		{
			pointing ptg (pi/2-(asin(2/pi*(asin(v) + v*sqrt((1-v)*(1+v))))),
				-pi/2*u/max(sqrt((1-v)*(1+v)),1e-6));
			double v = grid_val(rot.Transform(ptg.to_vec3()), lat_interval, lon_interval, linewidth*delta);
			if(v) img(i,j) = blend(img(i,j), Pixel(0,0,0), v);
		}
	}
}

double quantile(const vector<double> & dist, double q)
{
	int n = dist.size();
	int i = max(0,min(n-1,int(q*n)));
	return dist.empty() ? 0 : dist[i];
}

bool begins(const string & a, const string & ref, int n)
{
	if(ref.size() < n || a.size() < n) return false;
	if(a.size() > ref.size()) return false;
	int i;
	for(i = 0; i < a.size(); i++)
		if(a[i] != ref[i]) return false;
	return true;
}

void help()
{
	fprintf(stderr, "map2png: Produce a png image from a healpix fits map.\n"
		"Usage: map2png [options] map.[fits|hdf] img.png\n"
		" -maximum NUM   Color scale maximum\n"
		" -minimum NUM   Color scale minimum\n"
		" -range NUM     Color scale symmetric range\n"
		" -norange       Color scale automatic range (default)\n"
		" -gnomonic      Gnominic projection\n"
		" -mollweide     Mollweide projection (default)\n"
		" -linear        Linear color scale (default)\n"
		" -logarithmic   Logarithimic color scale\n"
		" -histogram     Histogram (maximum contrast) color scale\n"
		" -bar           Color bar on\n"
		" -nobar         Color bar off\n"
		" -verbose       Verbose output on\n"
		" -quiet         Verbose output off\n"
		" -grid NUM      Grid with specified interval in degrees\n"
		" -glat NUM      Grid interval in latitud\n"
		" -glon NUM      Grid interval in longitude\n"
		" -nogrid        Grid off (default)\n"
		" -color wmap|gray|COLORSPEC Color scheme to use (default wmap)\n"
		" -lw NUM        Grid line width in pixels\n"
		" -latitude NUM  Center of map in latitude\n"
		" -longitude NUM Center of map in longitude\n"
		" -signal INT    Column in fits table to use, can be specified multiple times\n"
		" -map INT       Which map to use (only for hdf, multiple ok)\n"
		" -xsz INT       Width of map in pixels\n"
		" -resolution NUM Map resolution in arcmin per pixel (default 1)\n"
		" -ticks INT     Levels of subtics in color bar\n"
		" -digits INT    Digits of numbers in color bar range\n"
		" -auto NUM      Quantiles to use for automatic color range (default 0)\n"
		" -ncol INT      Number of columns to use in grid layout (default 0: unlimited)\n"
		" -bg RRGGBB     Background color of image\n");
}

int main(int argc, char ** argv)
{
	try{
		// Defaults
		int pro = PRO_MOLL, trans = TRF_LIN;
		int xsize = 1024, ncol = 0, ticklevel = 2, dig = 7;
		double autoquant = 0;
		double lat = 0, lon = 0, res = 1;
		double gridlat = 10, gridlon = 10, gridline = 0.75;
		Vector<double,2> urange(0,0);
		Vector<bool,2> ugiven(false,false);
		bool bar = false, verbose = false, grid = false;
		vector<int> cols, submaps;
		vector<char*> args;
		Colorizer colorize = Colorizer::wmap;
		SimpleFont font = font_medium_bold;
		Pixel bgcol(128,128,128);

		// Read options
		for(char ** i = argv+1; *i; i++)
			if(begins(*i,"-maximum",4)) { urange(1) = atof(*++i); ugiven(1) = true; }
			else if(begins(*i,"-minimum",4)) { urange(0) = atof(*++i); ugiven(0) = true; }
			else if(begins(*i,"-range",6)) { urange = atof(*++i); urange(0) *= -1; ugiven = true; }
			else if(begins(*i,"-norange",8)) ugiven = false;
			else if(begins(*i,"-gnomonic",4)) pro = PRO_GNOM;
			else if(begins(*i,"-mollweide",4)) pro = PRO_MOLL;
			else if(begins(*i,"-linear",4)) trans = TRF_LIN;
			else if(begins(*i,"-logarithmic",4)) trans = TRF_LOG;
			else if(begins(*i,"-histogram",5)) trans = TRF_HIST;
			else if(begins(*i,"-bar", 4)) bar = true;
			else if(begins(*i,"-nobar", 6)) bar = false;
			else if(begins(*i,"-verbose", 2)) verbose = true;
			else if(begins(*i,"-quiet", 2)) verbose = false;
			else if(begins(*i,"-grid", 5)) { gridlat = gridlon = atof(*++i); grid = true; }
			else if(begins(*i,"-glat", 5)) { gridlat = atof(*++i); grid = true; }
			else if(begins(*i,"-glon", 5)) { gridlon = atof(*++i); grid = true; }
			else if(begins(*i,"-nogrid", 4)) grid = false;
			else if(begins(*i,"-color", 5))
			{
				string s = *++i;
				if(s == "wmap") colorize = Colorizer::wmap;
				else if(s == "gray" || s == "grey") colorize = Colorizer("0:000000,1:ffffff");
				else colorize = Colorizer(s);
			}
			else if(begins(*i,"-lw", 3)) gridline = atof(*++i);
			else if(begins(*i,"-latitude", 4)) lat = atof(*++i);
			else if(begins(*i,"-longitude", 4)) lon = atof(*++i);
			else if(begins(*i,"-signal", 4)) cols.push_back(atoi(*++i)-1);
			else if(begins(*i,"-map", 4)) submaps.push_back(atoi(*++i)-1);
			else if(begins(*i,"-xsz", 4) || begins(*i,"-size",5)) xsize = atoi(*++i);
			else if(begins(*i,"-resolution", 4)) res = atof(*++i);
			else if(begins(*i,"-ticks", 5)) ticklevel = atoi(*++i);
			else if(begins(*i,"-digits", 4)) dig = atoi(*++i);
			else if(begins(*i,"-auto", 5)) autoquant = atof(*++i);
			else if(begins(*i,"-ncol", 5)) ncol = atoi(*++i);
			else if(begins(*i,"-bg", 3)) bgcol = parse_color(*++i);
			else if(begins(*i,"-", 1)) { help(); return 1; }
			else args.push_back(*i);

		if(cols.empty()) cols.push_back(0);
		if(submaps.empty()) submaps.push_back(0);
		if(!ncol) ncol = cols.size();
		double t1, t2;

		// Read the maps
		if(args.size() < 1) { help(); return 1; }
		t1 = wall_time(); if(verbose) fprintf(stderr, "Reading maps ");
		string filename = args[0],
			ofilename = args.size() > 1 ? args[1] : set_file_extension(args[0],"png");
		Map maps = read_hmap<double>(filename);
		t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);

		int max_sig = *max_element(cols.begin(), cols.end());
		if(max_sig >= maps.size(1))
			serror("Asked for signal %d, but only %d exist%s!", max_sig+1, maps.size(1), maps.size(1) == 1 ? "s" : "");
		int max_map = *max_element(submaps.begin(), submaps.end());
		if(max_map >= maps.size(0))
			serror("Asked for map %d, but only %d exist%s!", max_map+1, maps.size(0), maps.size(0) == 1 ? "s" : "");

		// Project them down
		t1 = wall_time(); if(verbose) fprintf(stderr, "Projecting maps ");
		vector<Projection> projs;
		for(int m = 0; m < submaps.size(); m++)
		for(int c = 0; c < cols.size(); c++)
		{
			if(pro == PRO_MOLL)
				projs.push_back(pro_mollw(maps,submaps[m],cols[c],xsize,lon,lat));
			else
				projs.push_back(pro_gno(maps,submaps[m],cols[c],xsize,lon,lat,res));
		}
		t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);

		// Measure max/min/quantiles
		t1 = wall_time(); if(verbose) fprintf(stderr, "Measuring maps ");
		vector<double> dvals;
		for(int i = 0; i < projs.size(); i++)
		for(int j = 0; j < projs[i].size(); j++)
			if(ok_val(projs[i][j]))
			{
				// HACK
				if(trans == TRF_LOG) projs[i][j] = log(abs(projs[i][j]));
				dvals.push_back(projs[i][j]);
			}
		// HACK
		if(trans == TRF_LOG) trans = TRF_LIN;

		sort(dvals.begin(), dvals.end());
		Vector<double,2> drange;
		if(!dvals.empty()) drange = Vector<double,2>(dvals.front(), dvals.back());

		// If autoqantization is used, simulate user caps.
		if(autoquant)
		{
			if(!ugiven(0)) urange(0) = quantile(dvals, autoquant);
			if(!ugiven(1)) urange(1) = quantile(dvals, 1-autoquant);
			ugiven = 1;
		}

		// Cap them by the user ranges
		Vector<double,2> range;
		for(int i = 0; i < size(range); i++)
			range(i) = ugiven(i) ? urange(i) : drange(i);
		vector<double> vals;
		for(int i = 0; i < dvals.size(); i++)
			if(dvals[i] >= range(0) && dvals[i] <= range(1))
				vals.push_back(dvals[i]);

		t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);

		// Transform to the range [0,1] using the selected transformation
		t1 = wall_time(); if(verbose) fprintf(stderr, "Transforming maps ");
		Transformation * trf;
		switch(trans) {
			case TRF_LIN:  trf = new LinTrans(range(0), range(1)); break;
			case TRF_HIST: trf = new HistTrans(vals); break;
			default: serror("Unrecognized transformation!"); break;
		}
		for(int i = 0; i < projs.size(); i++)
			for(int j = 0; j < projs[i].size(); j++)
				if(ok_val(projs[i][j]))
					projs[i][j] = (*trf).forwards(projs[i][j]);
		t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);

		// Go from doubles to pixels by using the colorizer
		t1 = wall_time(); if(verbose) fprintf(stderr, "Coloring maps ");
		vector<Image> images;
		for(int i = 0; i < projs.size(); i++)
		{
			Image img(projs[i].extent());
			for(int j = 0; j < img.size(); j++)
			{
				double v = projs[i][j];
				if(is_nan(v))        img[j] = Pixel(255,255,255);
				else if(is_inf(v))   img[j] = Pixel(0,0,0);
				else if(is_hpbad(v)) img[j] = Pixel(128,128,128);
				else                 img[j] = colorize(v);
			}
			images.push_back(img);
		}
		Image blank(images[0].extent(),bgcol);
		t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);

		// Overlay the grid
		if(grid)
		{
			t1 = wall_time(); if(verbose) fprintf(stderr, "Drawing grid ");
			for(int i = 0; i < images.size(); i++)
				if(pro == PRO_MOLL)
					grid_moll(images[i], lon, lat, gridlon, gridlat, gridline);
				else
					grid_gno(images[i], lon, lat, res, gridlon, gridlat, gridline);
			t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);
		}

		// Make the color bar
		t1 = wall_time(); if(verbose) fprintf(stderr, "Making color bar ");
		Image colorbar = make_colorbar(trf, colorize, font, dig, ncol*xsize, ticklevel, drange, urange, ugiven);
		t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);

		// Create the final canvas
		t1 = wall_time(); if(verbose) fprintf(stderr, "Composing image ");
		int nrow = (images.size()+ncol-1)/ncol;
		int nx = images[0].size(0), ny = images[0].size(1);
		Image result(ncol*nx,nrow*ny + (bar ? colorbar.size(1) : 0), bgcol);
		for(int i = 0; i < ncol; i++)
		for(int j = 0; j < nrow; j++)
		{
			int k = i+j*ncol;
			copy(k < images.size() ? images[k] : blank, result, i*nx, j*ny);
		}
		copy(colorbar, result, 0, ny*nrow);
		t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);

		// And finally output
		t1 = wall_time(); if(verbose) fprintf(stderr, "Writing image ", result.size(0), result.size(1));
		write_png(ofilename.c_str(), result);
		t2 = wall_time(); if(verbose) fprintf(stderr, "%.6f\n", t2-t1);
	} catch(const Error & e) { fprintf(stderr, "%s\n", e.msg.c_str()); }
}

string set_file_extension(const string & s, const string & ext)
{
	size_t i = s.find_last_of(".");
	if(i == string::npos) return s + "." + ext;
	else return s.substr(0,i) + "." + ext;
}
