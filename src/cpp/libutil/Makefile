targets := libutil.a
objs    := sslice.o stoken.o spng.o scolorizer.o stimer.o serror.o seph.o hpixdata.o

all: $(targets)

libutil.a: $(objs)
	ar crs $@ $^

%.o: %.cpp
	$(CXX) -c -I. $(CFLAGS) $(HEALPIX_INCLUDE_CPP) -I../ephcom $<

clean:
	rm -f $(targets) $(objs)
