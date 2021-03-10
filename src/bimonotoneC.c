
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void bimonotoneC( int* rn, int* rm, double* rx, double* rw, int* maxiter, double* eps )
// Function bimonotone(),
// performs bivariate monotone regression
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
  size_t m = *rm;

  double** x = ( double** ) calloc( n + 1, sizeof( double ) );
  double* bx = ( double* ) calloc( ( n + 1 ) * ( m + 1 ), sizeof( double ) );
  for ( size_t i = 0, im1 = 0; i <= n; i++, im1++ ) x[i] = &bx[im1 * m];
  for ( size_t i = 1, k = 0; i <= n; i++ ) for ( size_t j = 1; j <= m; j++, k++ ) x[i][j] = rx[k];

  double** w = ( double** ) calloc( n + 1, sizeof( double ) );
  double* bw = ( double* ) calloc( ( n + 1 ) * ( m + 1 ), sizeof( double ) );
  for ( size_t i = 0, im1 = 0; i <= n; i++, im1++ ) w[i] = &bw[im1 * m];
  for ( size_t i = 1, k = 0; i <= n; i++ ) for ( size_t j = 1; j <= m; j++, k++ ) w[i][j] = rw[k];

  const size_t MAXITER = *maxiter;
  const double EPS = *eps;
  const double DELTA = 1.4901161193847656e-08;
  const double FRACT = 0.5;

  double wsum = 0.0;
  double wxsum = 0.0;
  double wmin = 1.7976931348623158e+308;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      double ww = w[i][j];
      if ( ww < DELTA ) continue;
      wsum += ww;
      wxsum += ww * x[i][j];
      if ( ww < wmin ) wmin = ww;
    }
  }
  if ( wsum < DELTA ) {
    free( x ); free( bx );
    free( w ); free( bw );
    return;
  }
  double wmean = wxsum / wsum;

  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      double ww = w[i][j];
      if ( ww < DELTA ) {
        w[i][j] = FRACT * wmin;
        x[i][j] = wmean;
      }
    }
  }

  double** g = ( double** ) calloc( n + 1, sizeof( double ) );
  double* bg = ( double* ) calloc( ( n + 1 ) * ( m + 1 ), sizeof( double ) );
  for ( size_t i = 0, im1 = 0; i <= n; i++, im1++ ) g[i] = &bg[im1 * m];

  double** c = ( double** ) calloc( n + 1, sizeof( double ) );
  double* bc = ( double* ) calloc( ( n + 1 ) * ( m + 1 ), sizeof( double ) );
  for ( size_t i = 0, im1 = 0; i <= n; i++, im1++ ) c[i] = &bc[im1 * m];

  double** resr = ( double** ) calloc( n + 1, sizeof( double ) );
  double* bresr = ( double* ) calloc( ( n + 1 ) * ( m + 1 ), sizeof( double ) );
  for ( size_t i = 0, im1 = 0; i <= n; i++, im1++ ) resr[i] = &bresr[im1 * m];

  double** resc = ( double** ) calloc( n + 1, sizeof( double ) );
  double* bresc = ( double* ) calloc( ( n + 1 ) * ( m + 1 ), sizeof( double ) );
  for ( size_t i = 0, im1 = 0; i <= n; i++, im1++ ) resc[i] = &bresc[im1 * m];

  double* xn = ( double* ) calloc( n + 1, sizeof( double ) );
  double* wn = ( double* ) calloc( n + 1, sizeof( double ) );
  double* xm = ( double* ) calloc( m + 1, sizeof( double ) );
  double* wm = ( double* ) calloc( m + 1, sizeof( double ) );

  double* hn = ( double* ) calloc( n + 1, sizeof( double ) );
  double* hm = ( double* ) calloc( m + 1, sizeof( double ) );

  size_t* idxn = ( size_t* ) calloc( n + 1, sizeof( size_t ) );
  size_t* idxm = ( size_t* ) calloc( m + 1, sizeof( size_t ) );

  size_t iflag = 0;
  size_t icount = 0;
  double rsqr = 0.0;
  double rsqc = 0.0;
  bool reset = true;

  for ( size_t iter = 1; iter <= MAXITER; iter++ ) {

    bool runrows = true;
    if ( reset ) {
      for ( size_t i = 1; i <= n; i++ ) {
        for ( size_t j = 1; j <= m; j++ ) {
          g[i][j] = x[i][j];
          resr[i][j] = 0.0;
          resc[i][j] = 0.0;
        }
      }
      icount = 0;
      if ( iflag == 1 ) runrows = false;
      else if ( iter == 3 && rsqc > rsqr ) runrows = false;
    }

    if ( runrows == true ) {
      size_t jcount = 0;
      for ( size_t i = 1; i <= n; i++ ) {
        for ( size_t j = 1; j <= m; j++ ) {
          xm[j] = hm[j] = g[i][j] - resr[i][j];
          wm[j] = w[i][j];
        }

        idxm[0] = 0;
        idxm[1] = 1;
        size_t b = 1;
        double xbm1 = xm[b];
        double wbm1 = wm[b];
        for ( size_t j = 2; j <= m; j++ ) {
          b++;
          double xb = xm[j];
          double wb = wm[j];
          if ( xbm1 > xb ) {
            b--;
            double sb = wbm1 * xbm1 + wb * xb;
            wb += wbm1;
            xb = sb / wb;
            while ( j < m && xb >= xm[j + 1] ) {
              j++;
              sb += wm[j] * xm[j];
              wb += wm[j];
              xb = sb / wb;
            }
            while ( b > 1 && xm[b - 1] > xb ) {
              b--;
              sb += wm[b] * xm[b];
              wb += wm[b];
              xb = sb / wb;
            }
          }
          xm[b] = xbm1 = xb;
          wm[b] = wbm1 = wb;
          idxm[b] = j;
        }
        size_t from = m;
        for ( size_t k = b; k > 0; k-- ) {
          const size_t to = idxm[k - 1] + 1;
          const double xk = xm[k];
          for ( size_t j = from; j >= to; j-- ) xm[j] = xk;
          from = to - 1;
        }

        size_t kcount = 0;
        for ( size_t j = 1; j <= m; j++ ) {
          double ord = xm[j];
          resr[i][j] = ord - hm[j];
          if ( fabs( ord - g[i][j] ) < EPS ) kcount = kcount + 1;
          g[i][j] = ord;
        }
        if ( kcount == m ) jcount = jcount + 1;
      }
      icount = icount + 1;
      if ( icount == 2 && iflag == 1 ) {
        rsqr = 0.0;
        for ( size_t i = 1; i <= n; i++ ) {
          for ( size_t j = 1; j <= m; j++ ) {
            rsqr += resr[i][j] * resr[i][j];
          }
        }
        iflag++;
      }
      if ( icount == 2 && iflag == 2 && iter == 2 ) {
        reset = true;
        continue;
      }
      if ( icount != 1 ) {
        if ( icount == MAXITER ) break;
        if ( jcount == n ) break;
      }
    }

    size_t lcount = 0;
    for ( size_t j = 1; j <= m; j++ ) {
      for ( size_t i = 1; i <= n; i++ ) {
        xn[i] = hn[i] = g[i][j] - resc[i][j];
        wn[i] = w[i][j];
      }

      idxn[0] = 0;
      idxn[1] = 1;
      size_t b = 1;
      double xbm1 = xn[b];
      double wbm1 = wn[b];
      for ( size_t i = 2; i <= n; i++ ) {
        b++;
        double xb = xn[i];
        double wb = wn[i];
        if ( xbm1 > xb ) {
          b--;
          double sb = wbm1 * xbm1 + wb * xb;
          wb += wbm1;
          xb = sb / wb;
          while ( i < n && xb >= xn[i + 1] ) {
            i++;
            sb += wn[i] * xn[i];
            wb += wn[i];
            xb = sb / wb;
          }
          while ( b > 1 && xn[b - 1] > xb ) {
            b--;
            sb += wn[b] * xn[b];
            wb += wn[b];
            xb = sb / wb;
          }
        }
        xn[b] = xbm1 = xb;
        wn[b] = wbm1 = wb;
        idxn[b] = i;
      }
      size_t from = n;
      for ( size_t k = b; k > 0; k-- ) {
        const size_t to = idxn[k - 1] + 1;
        const double xk = xn[k];
        for ( size_t i = from; i >= to; i-- ) xn[i] = xk;
        from = to - 1;
      }

      size_t mcount = 0;
      for ( size_t i = 1; i <= n; i++ ) {
        double ord = xn[i];
        resc[i][j] = ord - hn[i];
        if ( fabs( ord - g[i][j] ) < EPS ) mcount = mcount + 1;
        g[i][j] = ord;
      }
      if ( mcount == n ) lcount = lcount + 1;
    }
    icount = icount + 1;
    if ( icount == 2 && iflag == 0 ) {
      rsqc = 0.0;
      for ( size_t i = 1; i <= n; i++ ) {
        for ( size_t j = 1; j <= m; j++ ) {
          rsqc += resc[i][j] * resc[i][j];
        }
      }
      iflag++;
    }
    if ( icount == 2 && iflag == 1 ) {
      reset = true;
      continue;
    }
    if ( icount == 1 ) {
      reset = false;
      continue;
    }
    if ( lcount == m ) break;
    if ( icount == MAXITER ) break;
    reset = false;
  }

  free( x ); free( bx );
  free( w ); free( bw );

  free( c ); free( bc );
  free( resr ); free( bresr );
  free( resc ); free( bresc );

  free( xn );
  free( wn );
  free( xm );
  free( wm );

  free( hn );
  free( hm );

  free( idxn );
  free( idxm );

  for ( size_t i = 1, k = 0; i <= n; i++ ) for ( size_t j = 1; j <= m; j++, k++ ) rx[k] = g[i][j];

  free( g ); free( bg );
} // bimonotone
