#include <sslice.h>
#include <serror.h>
#include <stoken.h>
#include <limits>
#include <cstdlib>

using namespace std;

namespace skn {

const int Any = numeric_limits<int>::max();
typedef Vector<int,3> Vec;

Slice parse_slice(const string & spec)
{
	Slice res;
	vector<string> fields = tokenize(spec,",");
	for(int i = 0; i < fields.size(); i++)
	{
		vector<string> toks = tokenize_full(fields[i],":");
		vector<int> nums;
		for(int j = 0; j < toks.size(); j++)
			if(toks[j].size() == 0) nums.push_back(Any);
			else nums.push_back(atoi(toks[j].c_str()));
		Vec s;
		switch(nums.size())
		{
		case 0: s = Vec(Any,Any,1); break;
		case 1: s = nums[0] == Any ? Vec(Any,Any,1) : Vec(nums[0],1,1); break;
		case 2: s = Vec(nums[0],nums[1],1); break;
		case 3: s = Vec(nums[0],nums[1],nums[2]); break;
		}
		res.push_back(s);
	}
	return res;
}

Slice apply_slice(const Slice & old, const vector<int> lens)
{
	Slice res;
	for(int i = 0; i < lens.size(); i++)
	{
		if(i >= old.size()) res.push_back(Vec(0,lens[i],1));
		else
		{
			Vec o = old[i];
			int from = o[0], n = o[1], step = o[2];
			if(step == 0) serror("Zero step size is not allowed!");
			else if(step > 0)
			{
				if(from == Any) from = 0;
				if(n    == Any) n    = (lens[i]-from-1)/step+1;
			}
			else
			{
				if(from == Any) from = lens[i]-1;
				if(n    == Any) n    = from/step+1;
			}
			res.push_back(Vec(from,n,step));
		}
	}
	return res;
}

cut::cut(int from, int num, int step)
{ work.push_back(Vector<int,3>(from,num,step)); }
cut & cut::operator()(int from, int num, int step)
{
	work.push_back(Vector<int,3>(from,num,step));
	return *this;
}

cut::operator const Slice&() const { return work; }


}
