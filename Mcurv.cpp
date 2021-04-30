#include <rsf.hh>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"

using namespace std;
using namespace fdct_wrapping_ns;

void wedgesDist (int nbscales, vector<int> &nbangles, vector< vector<CpxNumMat> > &c1)
{
    /* Function to assert that the quantity of wedges per scales is correctly calculated */

    for(int sc=0; sc<nbscales; sc++)
    {
        assert(nbangles[sc] == c1[sc].size());
        //cerr<<sc<<"<-->";
        //cerr<<nbangles[sc]<<"<-->";
        //cerr<<c1[sc].size()<<endl;
    }
}

double abs2 (complex<double> valor)
{
    /* Function to calculate the square of the absolute value of a complex number */
    return pow(valor.real(),2) + pow(valor.imag(),2);
}

template <typename T> int sgn(T val) {
    /* Signal function for every type but complex<double> */
    return (T(0) < val) - (val < T(0));
}
complex<double> sgn(complex<double> val) {
    /* Signal function for comple<double> */
    return val / abs(val);
}



void copyToComplexDouble (int nz, int nx, float* d1, CpxNumMat &d1Cpx)
{
    /* Function to copy a float array to a complex<double> one, in order for it to be used in the
     * CurveLab suit */
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            d1Cpx(iz,ix) = cpx(d1[ix * nz + iz],0.0);
        }
    }
}

void copyFromComplexDouble (int nz, int nx, CpxNumMat &dataOut2Cpx, float* dataOut)
{
    /* Function to copy the real part from a complex<double> array in order to be able to write the
     * result in an output file */

    for (int ix = 0; ix < nx; ix++)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            dataOut[ix * nz + iz] = dataOut2Cpx(iz,ix).real();
        }
    }
}

void calc_nbangles (int nbscales, int ac, int nbangles_coarse, vector<int> &nbangles)
{
    /* Function to calculate the distribution of wedges (angles) per scale */
    int sc;
    if(ac == 1) // Use curvelet transform in the coarsest scale?
    {
        nbangles[0] = 1;
        for(sc = 1; sc < nbscales; sc++)
        {
            nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
        }
    }
    else if(ac == 0) // Or wavelet transform?
    {
        nbangles[0] = 1;
        for(sc = 1; sc < nbscales-1; sc++)
        {
            nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
        }
        nbangles[nbscales-1] = 1;
    }
    else
    {
        throw "Invalid ac parameter value (should be 0 or 1)";
    }
}


void filtDataBasic(int nbscales, vector<int> &nbangles,
                vector< vector<CpxNumMat> > &c1)
{
    /* Function to just calculate de filter with both the recursive division algorithm and the
     * normal division one. What differs this function from the pure filtData is that this one uses
     * an epsilon based on the max value of the whole curvelet dataset */
    cpx aux; // auxiliary variable to store the value of the filter before applying
    for (int is = 0; is < nbscales; is++)
    {
        for (int iw = 0;iw < nbangles[is]; iw++)
        {
            for (int ix = 0; ix < c1[is][iw].n(); ix++)
            {
                for (int iz = 0; iz < c1[is][iw].m(); iz++)
                {
                    // Soft-thresholding
                    c1[is][iw](iz,ix) = c1[is][iw](iz,ix) * cpx(20,0);
                }
            }
        }
    }
}


void checkingResidue (int nz, int nx, CpxNumMat d1Cpx, CpxNumMat dataOut2Cpx)
{
    /* Function to calculate the residue between the filtered and the migrated data */
    CpxNumMat e(nz,nx);
    for(int i=0; i<nz; i++)
       for(int j=0; j<nx; j++)
          e(i,j) = d1Cpx(i,j) - dataOut2Cpx(i,j);
    cerr<<"accuracy of inversion "<<sqrt(energy(e)/(nz*nx))<<endl;
}

int main(int argc, char* argv[])
{
    /* Main program that reads and writes data and read input variables */
    bool verb;
    sf_init(argc,argv); // init RSF
    if(! sf_getbool("verb",&verb)) verb=0;
    int ix,iz,is,iw;

    // Parameters
    int nbscales;
    int nbangles_coarse;
    int ac;
    if (!sf_getint("nbs",&nbscales)) sf_error("Need nbs=");
    if (!sf_getint("nba",&nbangles_coarse)) sf_error("Need nba=");
    if (!sf_getint("ac",&ac)) sf_error("Need ac=");

    // Setting up I/O files
    sf_file Fd1=NULL, Fout=NULL; // I/O files
    Fd1 = sf_input("d1");
    Fout = sf_output("out");

    // R/W axes
    sf_axis ax,az;
    sf_axis ax2,az2;
    int nx,nz,dz,dx;
    int nx2,nz2,dz2,dx2;
    az = sf_iaxa(Fd1,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fd1,2); nx = sf_n(ax); dx = sf_d(ax);

    sf_oaxa(Fout,az,1);
    sf_oaxa(Fout,ax,2);

    // Allocating I/O arrays
    float *d1 = new float[nz*nx]; sf_floatread(d1,nz*nx,Fd1);
    float *dataOut = new float[nz*nx];

    if(verb) cout<<stderr<<endl;

    // "Changing" type of array from float to complex
    CpxNumMat d1Cpx(nz,nx);
    copyToComplexDouble (nz, nx, d1, d1Cpx);

    // Forward transform
    vector< vector<CpxNumMat> > c1;  //vector<int> extra;
    fdct_wrapping(nz, nx, nbscales, nbangles_coarse, ac, d1Cpx, c1);

    // Calculating number of wedges per scale
    vector<int> nbangles(nbscales);
    try
    {
        calc_nbangles (nbscales, ac, nbangles_coarse, nbangles);
    }
    catch(const char* msg)
    {
        cerr<<msg<<endl;
        return 1;
    }

    // Printing scale wedges distribution
    wedgesDist (nbscales, nbangles, c1);

    //filtDataBasic(nbscales, nbangles, c1);

    // Inverse transform
    CpxNumMat dataOut2Cpx(nz,nx); clear(dataOut2Cpx);
    ifdct_wrapping(nz, nx, nbscales, nbangles_coarse, ac, c1, dataOut2Cpx);

    checkingResidue (nz, nx, d1Cpx, dataOut2Cpx);

    // "Changing" type of array from complex to float
    copyFromComplexDouble (nz, nx, dataOut2Cpx, dataOut);

    sf_floatwrite(dataOut,nz*nx,Fout);

    if(verb) cerr<<stderr<<endl;

    // Freeing allocated arrays
    delete d1;
    delete dataOut;

    sf_close();

    return 0;
}


