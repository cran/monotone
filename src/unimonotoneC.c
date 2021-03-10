
#include <stdlib.h>
#include <string.h>

void unimonotoneC( int* rn, double* x, double* w )
// Function unimonotone(),
// performs unimodal monotone regression
// Copyright (C) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// This function is free software:
// you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this function.
// If not, see <https://www.gnu.org/licenses/>.
{
  size_t n = *rn;

  double* lx = ( double* ) calloc( n + 1, sizeof( double ) );
  double* lw = ( double* ) calloc( n + 1, sizeof( double ) );
  int* ls = ( int* ) calloc( n + 1, sizeof( int ) );
  double* le = ( double* ) calloc( n + 1, sizeof( double ) );
  memcpy( &lx[1], x, n * sizeof( x[1] ) );
  memcpy( &lw[1], w, n * sizeof( w[1] ) );

  {
    ls[1] = 1;
    double xim1 = lx[1];
    double wim1 = lw[1];
    double wxx = wim1 * xim1 * xim1;
    double bxx = wxx;
    le[1] = wxx - bxx;
    for ( size_t i = 2; i <= n; i++ ) {
      double xi = lx[i];
      double wi = lw[i];
      size_t lastb = i - 1;
      int csize = 1;
      wxx += wi * xi * xi;
      if ( xim1 > xi ) {
        csize += ls[lastb];
        bxx -= lw[lastb] * lx[lastb] * lx[lastb];
        double sumx = wim1 * xim1 + wi * xi;
        wi += wim1;
        xi = sumx / wi;
        le[i] = 1.7976931348623158e+308;
        while ( i < n && xi > lx[i + 1] ) {
          i++;
          le[i] = 1.7976931348623158e+308;
          csize++;
          sumx += lw[i] * lx[i];
          wi += lw[i];
          xi = sumx / wi;
          wxx += lw[i] * lx[i] * lx[i];
        }
        lastb -= ls[lastb];
        while ( lastb > 0 && lx[lastb] > xi ) {
          bxx -= lw[lastb] * lx[lastb] * lx[lastb];
          csize += ls[lastb];
          sumx += lw[lastb] * lx[lastb];
          wi += lw[lastb];
          xi = sumx / wi;
          lastb -= ls[lastb];
        }
        bxx += sumx * sumx / wi;
      }
      else {
        bxx += wi * xi * xi;
        le[i] = wxx - bxx;
      }
      lx[i] = xim1 = xi;
      lw[i] = wim1 = wi;
      ls[i] = csize;
    }
  }

  double* rx = ( double* ) calloc( n + 1, sizeof( double ) );
  double* rw = ( double* ) calloc( n + 1, sizeof( double ) );
  int* rs = ( int* ) calloc( n + 1, sizeof( int ) );
  double* re = ( double* ) calloc( n + 1, sizeof( double ) );
  memcpy( &rx[1], x, n * sizeof( x[1] ) );
  memcpy( &rw[1], w, n * sizeof( w[1] ) );

  {
    rs[n] = 1;
    double xip1 = rx[n];
    double wip1 = rw[n];
    double wxx = wip1 * xip1 * xip1;
    double bxx = wxx;
    re[n] = wxx - bxx;
    for ( size_t i = n - 1; i >= 1; i-- ) {
      double xi = rx[i];
      double wi = rw[i];
      size_t lastb = i + 1;
      int csize = 1;
      wxx += wi * xi * xi;
      if ( xip1 > xi ) {
        csize += rs[lastb];
        bxx -= rw[lastb] * rx[lastb] * rx[lastb];
        double sumx = wip1 * xip1 + wi * xi;
        wi += wip1;
        xi = sumx / wi;
        re[i] = 1.7976931348623158e+308;
        while ( i > 1 && xi > rx[i - 1] ) {
          i--;
          re[i] = 1.7976931348623158e+308;
          csize++;
          sumx += rw[i] * rx[i];
          wi += rw[i];
          xi = sumx / wi;
          wxx += rw[i] * rx[i] * rx[i];
        }
        lastb += rs[lastb];
        while ( lastb <= n && rx[lastb] > xi ) {
          bxx -= rw[lastb] * rx[lastb] * rx[lastb];
          csize += rs[lastb];
          sumx += rw[lastb] * rx[lastb];
          wi += rw[lastb];
          xi = sumx / wi;
          lastb += rs[lastb];
        }
        bxx += sumx * sumx / wi;
      }
      else {
        bxx += wi * xi * xi;
        re[i] = wxx - bxx;
      }
      rx[i] = xip1 = xi;
      rw[i] = wip1 = wi;
      rs[i] = csize;
    }
  }

  int mode = 0;
  double work = 0.0;
  double error = 1.7976931348623158e+308;
  for ( int i = 1; i <= n; i++ ) {
    if ( ls[i] == 1 && rs[i] == 1 && error > ( work = le[i] + re[i] ) ) {
      mode = i;
      error = work;
    }
  }

  if ( mode != 0 ) {
    double* xx = &x[-1];
    int lk = mode - 1;
    while ( lk > 0 ) {
      const double xk = lx[lk];
      const size_t csize = ls[lk];
      for ( size_t i = 1; i <= csize; i++, lk-- ) xx[lk] = xk;
    }
    // x[mode] = 0.5 * ( lx[mode] + rx[mode] ); // same
    size_t rk = mode + 1;
    while ( rk < n ) {
      const double xk = rx[rk];
      const size_t csize = rs[rk];
      for ( size_t i = 1; i <= csize; i++, rk++ ) xx[rk] = xk;
    }
  }

  free( lx );
  free( lw );
  free( ls );
  free( le );

  free( rx );
  free( rw );
  free( rs );
  free( re );
} // unimonotone
