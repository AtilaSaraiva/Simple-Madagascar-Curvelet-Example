host=$(shell hostname)
ifeq ($(host),JurosComposto)
    LDFLAGS= -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsf++ -lrsf -lm -ltirpc -lfftw3f -lfftw3 -O3

endif
ifeq ($(host),marreca)
    LDFLAGS= -I$(RSFROOT)/include -L$(RSFROOT)/lib -lrsf++ -lrsf -lm -lfftw3f -lfftw3 -O3
endif

sfcurvelet: Mcurv.o
	g++ -o sfcurvelet Mcurv.o $(FDCT)/fdct_wrapping_cpp/src/libfdct_wrapping.a -fPIC -L$(FFTW)/lib -lfftw $(LDFLAGS)

Mcurv.o: Mcurv.cpp
	g++ -g -Wall -W -Wno-sign-compare -Wno-unused-label -MMD -fPIC -I$(FFTW)/include -I$(FDCT)/fdct_wrapping_cpp/src -O4 -DNDEBUG $(LDFLAGS) -c Mcurv.cpp

clean:
	rm *.d *.o sfcurvelet
