// Simple program to fill holes in the Era Interim SST dataset.
// It does use a simple distance weight interpolation to set a temperature on
// land points from neighbour sea points.

#include <netcdf.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>

using namespace std;

static int nlat = 121;
static int nlon = 240;
static int ntot = nlon*nlat;
static int nsize = ntot*sizeof(short);
static int ntimes = 29824;
//static int ntimes = 1;
static short missval = -32767;
static int spole = 109;
void fillholes(short *a, short *b, int ni, int nj);

int main(int argc, char *argv[])
{
  short *iobuf = new short[ntot];
  short *wkbuf = new short[ntot*3];

  NcFile *newf = new NcFile("sstERAIN.1989-2009.nc",NcFile::Write);
  NcVar *ncvtmp = newf->get_var("sst");

  for (int i = 0; i < ntimes; i ++)
  {
    // Read data
    ncvtmp->set_cur(i,0,0);
    ncvtmp->get(iobuf,1,nlat,nlon);

    // Triple data for longitude wrapping
    for (int j = 0; j < nlat; j ++)
    {
      memcpy(wkbuf+(j*3)*nlon+0*nlon,iobuf+j*nlon,nlon*sizeof(short));
      memcpy(wkbuf+(j*3)*nlon+1*nlon,iobuf+j*nlon,nlon*sizeof(short));
      memcpy(wkbuf+(j*3)*nlon+2*nlon,iobuf+j*nlon,nlon*sizeof(short));
    }

    // Fill holes
    fillholes(wkbuf,iobuf,nlat,nlon);

    // Write back data
    ncvtmp->put(iobuf,1,nlat,nlon);
    newf->sync();
    std::cout << "Time : " << i+1 << std::endl;
  }

  newf->close();

  delete [] iobuf;
  delete [] wkbuf;

  return 0;
}

void fillholes(short *a, short*b, int ni, int nj)
{
  int li1 = 0;
  int li2 = ni-1;
  int lj1 = -nj;
  int lj2 = 2*nj-1;

  for (int i = 0; i < spole; i ++)
  {
    // Copy the line
    memcpy(b+i*nj,a+i*3*nj+nj,nj*sizeof(short));
    for (int j = 0; j < nj; j ++)
    {
      int ib = i*3*nj+nj+j;
      if (a[ib] == missval)
      {
        // African lakes fix
//        if ( i > 56 && i < 72 && j > 16 && j < 26 )
//        {
//          double twgt = 0.0;
//          double val = 0.0;
//          for (int ii = 56; ii < 72; ii ++)
//          {
//            for (int jj = 16; jj < 26; jj ++)
//            {
//              int ib = ii*3*nj+nj+jj;
//              if ( a[ib] != missval )
//              {
//                double wgt = 1.0 / 
//                  sqrt(double((ii-i)*(ii-i))+double((jj-j)*(jj-j)));
//                val = val + wgt*double(a[ib]);
//                twgt = twgt + wgt;
//              }
//            }
//          }
//          b[i*nj+j] = short(val/twgt);
//          continue;
//        }
        // Great lakes fix
//        if ( i > 26 && i < 36 && j > 177 && j < 193 )
//        {
//          double twgt = 0.0;
//          double val = 0.0;
//          for (int ii = 26; ii < 36; ii ++)
//          {
//            for (int jj = 177; jj < 193; jj ++)
//            {
//              int ib = ii*3*nj+nj+jj;
//              if ( a[ib] != missval )
//              {
//                double wgt = 1.0 / 
//                  sqrt(double((ii-i)*(ii-i))+double((jj-j)*(jj-j)));
//                val = val + wgt*double(a[ib]);
//                twgt = twgt + wgt;
//              }
//            }
//          }
//          b[i*nj+j] = short(val/twgt);
//          continue;
//        }
        // Caspian & black sea fix
//        if ( i > 25 && i < 38 && j > 30 && j < 44 )
//        {
//          double twgt = 0.0;
//          double val = 0.0;
//          for (int ii = 25; ii < 38; ii ++)
//          {
//            for (int jj = 30; jj < 44; jj ++)
//            {
//              int ib = ii*3*nj+nj+jj;
//              if ( a[ib] != missval )
//              {
//                double wgt = 1.0 / 
//                  sqrt(double((ii-i)*(ii-i))+double((jj-j)*(jj-j)));
//                val = val + wgt*double(a[ib]);
//                twgt = twgt + wgt;
//              }
//            }
//          }
//          b[i*nj+j] = short(val/twgt);
//          continue;
//        }
        // China lakes fix
//        if ( i > 22 && i < 32 && j > 66 && j < 76 )
//        {
//          double twgt = 0.0;
//          double val = 0.0;
//          for (int ii = 22; ii < 32; ii ++)
//          {
//            for (int jj = 66; jj < 76; jj ++)
//            {
//              int ib = ii*3*nj+nj+jj;
//              if ( a[ib] != missval )
//              {
//                double wgt = 1.0 / 
//                  sqrt(double((ii-i)*(ii-i))+double((jj-j)*(jj-j)));
//                val = val + wgt*double(a[ib]);
//                twgt = twgt + wgt;
//              }
//            }
//          }
//          b[i*nj+j] = short(val/twgt);
//          continue;
//        }
        bool closed = false;
        double twgt , val;
        int id = 3;
        while (! closed)
        {
          int nd = (id-1)*4; // Perimetral points
          int ipp = 0;       // Number of perimetral points found
          val  = 0.0;
          twgt = 0.0;
          int hid = id/2;
          int il1 = i-hid;
          int il2 = i+hid;
          int jl1 = j-hid;
          int jl2 = j+hid;
          int minpp = (nd/5);  // At least more than one fourth of perim points
          for (int ii = il1; ii <= il2; ii ++)
          {
            if ( ii < li1 || ii > li2 ) // If outside of data
            {
              if (ii == il1 || ii == il2)
                ipp = ipp+id; // whole line of perim points
              else
                ipp = ipp+2;  // just two perim points
              continue; // Goto next line
            }
            for (int jj = jl1; jj <= jl2; jj ++)
            {
              if ( jj < lj1 || jj > lj2 ) // If outside of data
              {
                if (jj == jl1 || jj == jl2) ipp ++; // Perim point
                continue;
              }
              int iib = ii*3*nj+nj+jj;
              if (a[iib] != missval)
              {
                double wgt = 1.0 /
                     sqrt(double((ii-i)*(ii-i))+double((jj-j)*(jj-j)));
                val = val + (double(a[iib])*wgt);
                twgt = twgt + wgt;
                // If this is a perim point count it
                if (jj == jl1 || jj == jl2 || ii == il1 || ii == il2) ipp ++;
              }
            }
          }
          if (ipp > minpp) closed = true;
          id = id+2;
        }
        if (twgt > 0.0) b[i*nj+j] = short(val/twgt);
      }
    }
  }
  for (int i = spole; i < ni; i ++)
  {
    for (int j = 0; j < nj; j ++)
    {
      b[i*nj+j] = -27444;
    }
  }
  return;
}
