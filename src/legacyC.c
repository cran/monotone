
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void fake( int* n, double* x, double* w ) 
// author(s): F.M.T.A. Busing
// language : C
{
}

void fitm( int* n, double* x, double* w ) 
// author(s): J.B. Kruskal
// origin   : subroutine FITM, MDSCAL, version 5 (October 1971), unchanged from version 4 (January 1968)
// language : Fortran
{
  x--;
  w--;
  int* lblock = ( int* ) calloc( *n + 1, sizeof( int ) );
  double* dhat = ( double* ) calloc( *n + 1, sizeof( double ) );
  for ( size_t i = 1; i <= *n; i++ ) lblock[i] = 1;
  for ( size_t i = 1; i <= *n; i++ ) dhat[i] = x[i] * w[i];
  int ma = 1;
  int mb = 0;
  while ( true ) {
    int lud = 2;
    int nsatis = 0;
    while ( true ) {
      int k = lblock[ma];
      mb = ma + k - 1;
      double wt = ( k - 1 == 0 ? w[mb] : dhat[mb] );
      double dav = dhat[ma] / wt;
      if ( lud == 1 ) {
        if ( ma - 1 == 0 ) nsatis = nsatis + 1;
        else {
          int mbd = ma - 1;
          int kd = lblock[mbd];
          int mad = mbd - kd + 1;
          double wtd = ( kd - 1 == 0 ? w[mbd] : dhat[mbd] );
          double davd = dhat[mad] / wtd;
          if ( davd - dav < 0.0 ) nsatis = nsatis + 1;
          else {
            int knew = k + kd;
            lblock[mad] = knew;
            lblock[mb] = knew;
            double dtonew = dhat[mad] + dhat[ma];
            dhat[mad] = dtonew;
            dhat[mb] = wt + wtd;
            nsatis = 0;
            ma = mad;
          }
        }
      }
      if ( lud == 2 ) {
        if ( mb - *n == 0 ) nsatis = nsatis + 1;
        else {
          int mau = mb + 1;
          int ku = lblock[mau];
          int mbu = mau + ku - 1;
          double wtu = ( ku - 1 == 0 ? w[mbu] : dhat[mbu] );
          double davu = dhat[mau] / wtu;
          if ( dav - davu < 0.0 ) nsatis = nsatis + 1;
          else {
            int knew = k + ku;
            lblock[ma] = knew;
            lblock[mbu] = knew;
            double dtonew = dhat[ma] + dhat[mau];
            dhat[ma] = dtonew;
            dhat[mbu] = wt + wtu;
            nsatis = 0;
          }
        }
      }
      lud = 3 - lud;
      if( nsatis - 1 > 0 ) break;
    }
    if( mb == *n ) break;
    ma = mb + 1;
  }
  ma = 1;
  while ( ( ma - *n - 1 ) < 0 ) {
    int k = lblock[ma];
    int mb = ma + k - 1;
    if ( k - 1 > 0 ) {
      double temp1 = dhat[ma] / dhat[mb];
      for ( int m = ma; m <= mb; m++ ) x[m] = temp1;
    }
    else dhat[ma] = x[ma];
    ma = mb + 1;
  }
  free( lblock );
  free( dhat );
} // fitm

void wmrmnh( int* n, double* x, double* w )
// author(s): Ernst van Waning and Ineke Stoop
// origin   : Ernst van Waning (1976). A set of programs to perform a Kruskal-type monotone regression.
// language : Fortran
{
  x--; 
  w--;
  size_t lovbkh = 0;
  long double wovbkh = 0.0L;
  size_t luph = 0;
  long double trialv = 0.0L;
  for ( size_t iup = 2; iup <= *n; iup++ ) {
    if ( x[iup] < x[iup - 1] ) {
      long double sds = x[iup] * w[iup];
      long double sw = w[iup];
      size_t idown = iup;
      do {
        idown -= 1;
        if ( luph == idown ) {
          sds += wovbkh * x[luph];
          sw += wovbkh;
          idown = idown - lovbkh;
        }
        else {
          sds += x[idown] * w[idown];
          sw += w[idown];
        }
        trialv = sds / sw;
        if ( idown == 1 ) break;
      } while ( x[idown - 1] > trialv );
      wovbkh = 0.0L;
      for ( size_t i = idown; i <= iup; i++ ) {
        x[i] = trialv;
        wovbkh += w[i];
      }
      lovbkh = iup - idown;
      luph = iup;
    }
  }
} // wmrmnh

void amalgm( int* n, double* x, double* w )
// author(s): G.W. Cran
// origin   : Algorithm AS 149: Amalgamation of Means in the Case of Simple Ordering
//            Journal of the Royal Statistical Society.Series ( Applied Statistics ), Vol. 29, No. 2 ( 1980 ), pp. 209 - 211
//            Wiley for the Royal Statistical Society
//            http://www.jstor.org/stable/2986312
// language : Fortran
{
  x--; 
  w--;
  const long double tol = 1.0E-6;
  double* xt = ( double* ) calloc( *n + 1, sizeof( double ) );
  double* wt = ( double* ) calloc( *n + 1, sizeof( double ) );
  for ( size_t i = 1; i <= *n; i++ ) {
    xt[i] = x[i];
    wt[i] = w[i];
  }
  size_t m = *n;
  size_t i = 1;
  while ( true ) {
    if ( i == m ) break;
    if ( i != m && xt[i] > xt[i + 1] ) {
      size_t i1 = i + 1;
      long double ww = wt[i] + wt[i1];
      xt[i] = ( wt[i] * xt[i] + wt[i1] * xt[i1] ) / ww;
      wt[i] = ww;
      size_t mm1 = m - 1;
      if ( i1 != m ) {
        for ( size_t j = i1; j <= mm1; j++ ) {
          size_t j1 = j + 1;
          xt[j] = xt[j1];
          wt[j] = wt[j1];
        }
      }
      m = mm1;
      if ( m == 1 ) break;
      while ( i != 1 && xt[i - 1] > xt[i] ) {
        size_t im1 = i - 1;
        long double ww = wt[im1] + wt[i];
        xt[im1] = ( wt[im1] * xt[im1] + wt[i] * xt[i] ) / ww;
        wt[im1] = ww;
        size_t mm1 = m - 1;
        if ( i != m ) {
          for ( size_t j = i; j <= mm1; j++ ) {
            size_t j1 = j + 1;
            xt[j] = xt[j1];
            wt[j] = wt[j1];
          }
        }
        i = im1;
        m = mm1;
        if ( m == 1 ) break;
      }
    }
    else i = i + 1;
  }
  size_t j = 0;
  size_t i1 = 1;
  for ( size_t i = 1; i <= m; i++ ) {
    long double s = 0.0L;
    for ( j = i1; j <= *n; j++ ) {
      s = s + w[j];
      x[j] = xt[i];
      if ( fabsl( s - wt[i] ) < tol ) break;
    }
    i1 = j + 1;
  }
  free( xt );
  free( wt );
} // amalgm

void pav( int* n, double* x, double* w )
// author(s): Bril (after Dykstra and Robertson)
// origin   : Algorithm AS 206.1: PAV() - apply pool adjacent violators theorem
//            Part of algorithm AS 206 Applied Statistics ( 1984 ) vol.33, no.3 
// language : Fortran
{
  x--;
  w--;
  size_t* nw = ( size_t* ) calloc( *n + 1, sizeof( size_t ) );
  double* fx = ( double* ) calloc( *n + 1, sizeof( double ) );
  double* pw = ( double* ) calloc( *n + 1, sizeof( double ) );
  double* w1 = ( double* ) calloc( *n + 1, sizeof( double ) );
  double* wt = ( double* ) calloc( *n + 1, sizeof( double ) );
  long double eps = 0.000001L;
  size_t nwc = *n;
  for ( size_t i = 1; i <= *n; i++ ) {
    nw[i] = 1;
    fx[i] = x[i];
    wt[i] = w[i];
    pw[i] = wt[i] * fx[i];
    w1[i] = w[i];
  }
  size_t icount = 0;
  size_t ibel = *n - 1;
  do {
    size_t i = 0;
    while ( true ) {
      i = i + 1;
      bool ct = false;
      while ( i <= ibel ) {
        size_t i1 = i + 1;
        ct = ( fx[i] - fx[i1] <= eps );
        if ( ct == true ) break;
        pw[i] += pw[i1];
        w1[i] += w1[i1];
        fx[i] = pw[i] / w1[i];
        nw[i] += nw[i1];
        nwc = nwc - 1;
        if ( i1 <= ibel ) {
          for ( size_t j = i1; j <= ibel; j++ ) {
            size_t j1 = j + 1;
            pw[j] = pw[j1];
            w1[j] = w1[j1];
            fx[j] = fx[j1];
            nw[j] = nw[j1];
          }
        }
        ibel = ibel - 1;
      }
      if ( ct == false ) break;
    }
    icount = 0;
    if ( ibel <= 0 ) break;
    for ( size_t l = 1; l <= ibel; l++ ) if ( fx[l] - fx[l + 1] <= eps ) icount++;
  }
  while ( icount != ibel );
  size_t j = 1;
  size_t jl = 1;
  size_t ju = nw[1];
  while ( true ) {
    for ( size_t l = jl; l <= ju; l++ ) x[l] = fx[j];
    j = j + 1;
    if ( j > nwc ) break;
    jl = ju + 1;
    ju += nw[j];
  }
  free( nw );
  free( fx );
  free( pw );
  free( w1 );
  free( wt );
} // pav

void isoreg( int *n, double* x, double* w )
// author(s): Gupta and Brian Ripley
// origin   : Isotonic regression, code simplified from VR_mds_fn() (part of MASS.c)
// language : C
{
  double* xc = ( double* ) calloc( *n + 1, sizeof( double ) );
  xc[0] = 0.0;
  long double tmp = 0.0L;
  for ( size_t i = 0; i < *n; i++ ) {
    tmp += x[i];
    xc[i + 1] = tmp;
  }
  size_t known = 0;
  size_t ip = 0;
  do {
    long double slope = 1.7E308L;
    for ( size_t i = known + 1; i <= *n; i++ ) {
      long double tmp = ( xc[i] - xc[known] ) / ( long double ) ( i - known );
      if ( tmp < slope ) {
        slope = tmp;
        ip = i;
      }
    }
    for ( size_t i = known; i < ip; i++ ) x[i] = ( xc[ip] - xc[known] ) / ( long double ) ( ip - known );
  }
  while ( ( known = ip ) < *n );
  free( xc );
} // isoreg

void iso_pava( int* n, double* x, double* w )
// author(s): Rolf Turner
// origin   : subroutine pava( y, w, kt, n )
//            R package Iso: Functions to Perform Isotonic Regression.
// language : Fortran
{
  x--;
  w--;
  int* kt = ( int* ) calloc( *n + 1, sizeof( int ) );
  bool same;
  for ( size_t i = 1; i <= *n; i++ ) kt[i] = i;
  while ( true ) {
    same = true;
    for ( size_t i = 2; i <= *n; i++ ) {
      if ( x[i - 1] > x[i] ) {
        size_t k1 = kt[i];
        size_t k2 = kt[i - 1];
        for ( size_t j = 1; j <= *n; j++ ) if ( kt[j] == k1 ) kt[j] = k2;
        long double wnew = w[i - 1] + w[i];
        long double ynew = ( w[i - 1] * x[i - 1] + w[i] * x[i] ) / wnew;
        for ( size_t j = 1; j <= *n; j++ ) {
          if ( kt[j] == k2 ) {
            x[j] = ynew;
            w[j] = wnew;
          }
        }
        same = false;
      }
    }
    if ( same ) break;
  }
  free( kt );
} // iso_pava

bool sorted( const size_t n, double* x )
{
  for ( size_t i = 2; i <= n; i++ ) if ( x[i - 1] > x[i] ) return false;
  return true;
} // sorted

void isotonic( int* n, double* x, double* w )
// author(s): Tom Kincaid
// origin   : Function isotonic in R Package mhweber / spsurvey
//            Home / GitHub / mhweber / spsurvey / R / isotonic.R
// language : C
{
  x--;
  w--;
  while ( !sorted( *n, x ) ) {
    size_t i = 1;
    while ( i < *n ) {
      size_t j = i;
      while ( x[j] >= x[j + 1] && j < *n ) j = j + 1;
      if ( i < j ) {
        long double mn = 0.0L;
        for ( size_t k = i; k <= j; k++ ) mn += x[k];
        mn /= ( long double ) ( j - i + 1 );
        for ( size_t k = i; k <= j; k++ ) x[k] = mn;
      }
      i = j + 1;
    }
  }
} // isotonic

void isomean( int *n, double* x, double* w )
// author(s): Kaspar Rufibach and Korbinian Strimmer
// origin   : isomean.c  (2007-07-06) by Korbinian Strimmer
//            ported from R code originally by Kaspar Rufibach / June 2004
//            part of fdrtool library for R and related languages.
// language : C
{
  size_t nn = *n;
  size_t* k = ( size_t* ) calloc( ( size_t ) nn, sizeof( size_t ) );
  double* gew = ( double* ) calloc( ( size_t ) nn, sizeof( double ) );
  double* ghat = ( double* ) calloc( ( size_t ) nn, sizeof( double ) );
  size_t c = 0;
  k[c] = 0;
  gew[c] = w[0];
  ghat[c] = x[0];
  for ( size_t j = 1; j < nn; j++ ) {
    c = c + 1;
    k[c] = j;
    gew[c] = w[j];
    ghat[c] = x[j];
    while ( ghat[c - 1] >= ghat[c] ) {
      long double neu = gew[c] + gew[c - 1];
      ghat[c - 1] = ghat[c - 1] + ( gew[c] / neu ) * ( ghat[c] - ghat[c - 1] );
      gew[c - 1] = neu;
      c = c - 1;
      if ( c == 0 ) break;
    }
  }
  while ( nn >= 1 ) {
    for ( size_t j = k[c]; j < nn; j++ ) {
      ghat[j] = ghat[c];
    }
    nn = k[c];
    c = c - 1;
  }
  for ( size_t j = 0; j < *n; j++ ) x[j] = ghat[j];
  free( k );
  free( gew );
  free( ghat );
} // isomean

void pooled_pava( int *n, double* x, double* w )
// author(s): Pedregosa and Tulloch
// origin   : wikipava, Tulloch, 2014, 2016
//            https://github.com/ajtulloch/Isotonic.jl/blob/master/src/pooled_pava.jl
// language : Julia
{
  x--;
  w--;
  size_t j = 1;
  size_t* S = ( size_t* ) calloc( *n + 1, sizeof( size_t ) );
  S[0] = 0;
  S[1] = 1;
  double* xdash = ( double* ) calloc( *n + 1, sizeof( double ) );
  double* wdash = ( double* ) calloc( *n + 1, sizeof( double ) );
  xdash[1] = x[1];
  wdash[1] = w[1];
  for ( size_t i = 2; i <= *n; i++ ) {
    j = j + 1;
    xdash[j] = x[i];
    wdash[j] = w[i];
    while ( j > 1 && xdash[j] < xdash[j - 1] ) {
      xdash[j - 1] = ( wdash[j] * xdash[j] + wdash[j - 1] * xdash[j - 1] ) / ( wdash[j] + wdash[j - 1] );
      wdash[j - 1] = wdash[j] + wdash[j - 1];
      j = j - 1;
    }
    S[j] = i;
  }
  for ( size_t k = 1; k <= j; k++ ) {
    for ( size_t l = S[k - 1] + 1; l <= S[k]; l++ ) x[l] = xdash[k];
  }
  free( S );
  free( xdash );
  free( wdash );
} // pooled_pava

void linear_pava( int* n, double* x, double* w )
// author(s): Tulloch and Varoquaux
// origin   : scikit - learn implementation by Andrew Tulloch, 2014
//            https://github.com/ajtulloch/Isotonic.jl/blob/master/src/pooled_pava.jl
// language : Julia and C
{
  size_t nm1 = *n - 1;
  while ( true ) {
    size_t i = 0;
    bool pooled = false;
    while ( i < nm1 ) {
      size_t k = i;
      while ( k < nm1 && x[k] >= x[k + 1] ) k++;
      if ( x[i] != x[k] ) {
        long double numerator = 0.0L;
        long double denominator = 0.0L;
        for ( size_t j = i; j <= k; j++ ) {
          numerator += x[j] * w[j];
          denominator += w[j];
        }
        for ( size_t j = i; j <= k; j++ ) x[j] = numerator / denominator;
        pooled = true;
      }
      i = k + 1;
    }
    if ( !pooled ) break;
  }
} // linear_pava

void inplace_pava( int* n, double* x, double* w ) 
// author(s): Nelle Varoquaux, Andrew Tulloch, Antony Lee (April 2013)
// origin   : https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/_isotonic.pyx
// language : Python
{
  size_t* target = ( size_t* ) calloc( *n, sizeof( size_t ) );
  for ( size_t i = 0; i < *n; i++ ) target[i] = i;
  size_t i = 0;
  while ( i < *n ) {
    size_t k = target[i] + 1;
    if ( k == *n ) break;
    if ( x[i] < x[k] ) {
      i = k;
      continue;
    }
    double sum_wx = w[i] * x[i];
    double sum_w = w[i];
    while ( true ) {
      double prev_x = x[k];
      sum_wx += w[k] * x[k];
      sum_w += w[k];
      k = target[k] + 1;
      if ( k == *n || prev_x < x[k] ) {
        x[i] = sum_wx / sum_w;
        w[i] = sum_w;
        target[i] = k - 1;
        target[k - 1] = i;
        if ( i > 0 ) i = target[i - 1];
        break;
      }
    }
  }
  i = 0;
  while ( i < *n ) {
    size_t k = target[i] + 1;
    for ( size_t j = i + 1; j < k; j++ ) x[j] = x[i];
    i = k;
  }
  free( target );
} // inplace_pava

void md_pava( int* n, double* x, double* w )
// author(s): Maximilien Danisch (April 2016)
// origin   : http://bit.ly/maxdan94
// language : C++
{
  size_t* nag = ( size_t* ) malloc( *n * sizeof( size_t ) );
  double* val = ( double* ) malloc( *n * sizeof( double ) );
  nag[0] = 1;
  val[0] = x[0];
  size_t j = 0;
  for ( size_t i = 1; i < *n; i++ ) {
    j += 1;
    val[j] = x[i];
    nag[j] = 1;
    while ( ( j > 0 ) && ( val[j] < val[j - 1] ) ) {
      val[j - 1] = ( nag[j] * val[j] + nag[j - 1] * val[j - 1] ) / ( nag[j] + nag[j - 1] );
      nag[j - 1] += nag[j];
      j--;
    }
  }
  for ( size_t i = 0, l = 0; i <= j; i++ )
    for ( size_t k = 1; k <= nag[i]; k++, l++ ) x[l] = val[i];
  free( nag );
  free( val );
} // md_pava

void reg_1d_l2( int *n, double* x, double* w )
// author(s): Zhipeng Xu, Chenkai Sun, Aman Karunakaran, and Quentin Stout
// origin   : Quentin F. Stout; Unimodal Regression via Prefix Isotonic Regression Computational Statistics and Data Analysis 53 (2008), pp. 289-297;
//            Spouge, J., Wan, H. & Wilbur, W. Journal of Optimization Theory and Applications (2003) 117: 585-605 doi.org/10.1023/A:1023901806339
//            https://github.com/xzp1995/UniIsoRegression
// language : C++
{
  double* mean_vec = ( double* ) calloc( *n, sizeof( double ) );
  double* sumwy_vec = ( double* ) calloc( *n, sizeof( double ) );
  double* sumw_vec = ( double* ) calloc( *n, sizeof( double ) );
  size_t* left_vec = ( size_t* ) calloc( *n, sizeof( size_t ) );
  size_t* right_vec = ( size_t* ) calloc( *n, sizeof( size_t ) );
  for ( size_t i = 0; i < *n; ++i ) {
    sumwy_vec[i] = x[i] * w[i];
    sumw_vec[i] = w[i];
    left_vec[i] = i;
    right_vec[i] = i;
  }
  size_t vec_back = 0;
  for ( size_t j = 1; j < *n; ++j ) {
    bool flag = false;
    while ( mean_vec[j] <= mean_vec[vec_back] ) {
      left_vec[j] = left_vec[vec_back];
      sumwy_vec[j] = sumwy_vec[j] + sumwy_vec[vec_back];
      sumw_vec[j] = sumw_vec[j] + sumw_vec[vec_back];
      mean_vec[j] = sumwy_vec[j] / sumw_vec[j];
      if ( !vec_back ) {
        flag = true;
        break;
      }
      vec_back = vec_back - 1;
    }
    if ( !flag ) vec_back = vec_back + 1;
    left_vec[vec_back] = left_vec[j];
    right_vec[vec_back] = right_vec[j];
    sumwy_vec[vec_back] = sumwy_vec[j];
    sumw_vec[vec_back] = sumw_vec[j];
    mean_vec[vec_back] = mean_vec[j];
  }
  for ( size_t k = 0; k <= vec_back; ++k ) {
    for ( size_t l = left_vec[k]; l <= right_vec[k]; ++l ) {
      x[l] = mean_vec[k];
    }
  }
  free( mean_vec );
  free( sumwy_vec );
  free( sumw_vec );
  free( left_vec );
  free( right_vec );
} // reg_1d_l2

struct block {
  double value;
  double weight;
  int size;
  int previous;
  int next;
}; // block

void jbkpava ( int *n, double *x, double *w )
// author(s): Jan de Leeuw
// origin   : Exceedingly Simple Monotone Regression (2017)
// language : C
{
  struct block *blocks = calloc ( ( size_t ) *n, sizeof( struct block ) );
  for ( int i = 0; i < *n; i++ ) {
    blocks[i].value = x[i];
    blocks[i].weight = w[i];
    blocks[i].size = 1;
    blocks[i].previous = i - 1;
    blocks[i].next = i + 1;
  }
  int active = 0;
  do {
    bool upsatisfied = false;
    int next = blocks[active].next;
    if ( next == *n ) upsatisfied = true;
    else if ( blocks[next].value > blocks[active].value ) upsatisfied = true;
    if ( !upsatisfied ) {
      double ww = blocks[active].weight + blocks[next].weight;
      int nextnext = blocks[next].next;
      double wxactive = blocks[active].weight * blocks[active].value;
      double wxnext = blocks[next].weight * blocks[next].value;
      blocks[active].value = ( wxactive + wxnext ) / ww;
      blocks[active].weight = ww;
      blocks[active].size += blocks[next].size;
      blocks[active].next = nextnext;
      if ( nextnext < *n ) blocks[nextnext].previous = active;
      blocks[next].size = 0;
    }
    bool downsatisfied = false;
    int previous = blocks[active].previous;
    if ( previous == -1 ) downsatisfied = true;
    else if ( blocks[previous].value < blocks[active].value ) downsatisfied = true;
    if ( !downsatisfied ) {
      double ww = blocks[active].weight + blocks[previous].weight;
      int previousprevious = blocks[previous].previous;
      double wxactive = blocks[active].weight * blocks[active].value;
      double wxprevious = blocks[previous].weight * blocks[previous].value;
      blocks[active].value = ( wxactive + wxprevious ) / ww;
      blocks[active].weight = ww;
      blocks[active].size += blocks[previous].size;
      blocks[active].previous = previousprevious;
      if ( previousprevious > -1 ) blocks[previousprevious].next = active;
      blocks[previous].size = 0;
    }
    if ( ( blocks[active].next == *n ) && downsatisfied ) break;
    if ( upsatisfied && downsatisfied ) active = next;
  } while ( true );
  int k = 0;
  for ( int i = 0; i < *n; i++ ) {
    int blksize = blocks[i].size;
    if ( blksize > 0.0 ) {
      for ( int j = 0; j < blksize; j++ ) {
        x[k] = blocks[i].value;
        k++;  
      }
    }
  }
  free ( blocks );
} // jbkpava
