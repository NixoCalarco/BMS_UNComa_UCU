# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "lambert_w.h"

/******************************************************************************/

double lambert_w ( double x, int nb, int l )

/******************************************************************************/
/*
  Purpose:

    lambert_w() approximates the W function.

  Discussion:

    The call will fail if the input value X is out of range.
    The range requirement for the upper branch is:
      -exp(-1) <= X.
    The range requirement for the lower branch is:
      -exp(-1) < X < 0.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    21 June 2023

  Author:

    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
    Patricia Culligan-Hensley.
    This version by John Burkardt.

  Reference:

    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
    Algorithm 743: WAPR - A Fortran routine for calculating real 
    values of the W-function,
    ACM Transactions on Mathematical Software,
    Volume 21, Number 2, June 1995, pages 172-181.

  Input:

    double x: the argument.

    int nb: indicates the desired branch.
    * 0, the upper branch;
    * nonzero, the lower branch.

    int l: indicates the interpretation of X.
    * 1, X is actually the offset from -(exp-1), so compute W(X-exp(-1)).
    * not 1, X is the argument; compute W(X);

  Output:

    double lambert_w: the approximate value of W(X).
*/
{
  double an2;
  double an3;
  double an4;
  double an5;
  double an6;
  double c13;
  double c23;
  double d12;
  double delx;
  double em;
  double em2;
  double em9;
  double eta;
  int i;
  int nbits;
  int niter;
  double reta;
  double s2;
  double s21;
  double s22;
  double s23;
  double t;
  double tb;
  double temp;
  double temp2;
  double ts;
  double value;
  double x0;
  double x1;
  double xx;
  double zl;
  double zn;

  value = 0.0;
/*
  Various mathematical constants.
*/
  nbits = 52;
  niter = 1;
  em = -exp ( -1.0 );
  em9 = -exp ( -9.0 );
  c13 = 1.0 / 3.0;
  c23 = 2.0 * c13;
  em2 = 2.0 / em;
  d12 = -em2;
  tb = pow ( 0.5, nbits );
  x0 = pow ( tb, 1.0 / 6.0 ) * 0.5;
  x1 = ( 1.0 - 17.0 * pow ( tb, 2.0 / 7.0 ) ) * em;
  an3 = 8.0 / 3.0;
  an4 = 135.0 / 83.0;
  an5 = 166.0 / 39.0;
  an6 = 3167.0 / 3549.0;
  s2 = sqrt ( 2.0 );
  s21 = 2.0 * s2 - 3.0;
  s22 = 4.0 - 3.0 * s2;
  s23 = s2 - 2.0;

  if ( l == 1 )
  {
    delx = x;

    if ( delx < 0.0 )
    {
      exit ( 1 );
    }

    xx = x + em;
  }
  else
  {
    if ( x < em )
    {
      exit ( 1 );
    }
    else if ( x == em )
    {
      value = -1.0;
      return value;
    }
    xx = x;
    delx = xx - em;
  }
/*
  Calculations for Wp.
*/
  if ( nb == 0 )
  {
    if ( fabs ( xx ) <= x0 )
    {
      value = xx / ( 1.0 + xx / ( 1.0 + xx 
        / ( 2.0 + xx / ( 0.6 + 0.34 * xx ))));
      return value;
    }
    else if ( xx <= x1 )
    {
      reta = sqrt ( d12 * delx );
      value = reta / ( 1.0 + reta / ( 3.0 + reta / ( reta 
        / ( an4 + reta / ( reta * an6 + an5 ) ) + an3 ) ) ) 
        - 1.0;
      return value;
    }
    else if ( xx <= 20.0 )
    {
      reta = s2 * sqrt ( 1.0 - xx / em );
      an2 = 4.612634277343749 * sqrt ( sqrt ( reta + 
        1.09556884765625 ));
      value = reta / ( 1.0 + reta / ( 3.0 + ( s21 * an2 
        + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0;
    }
    else
    {
      zl = log ( xx );
      value = log ( xx / log ( xx 
        / pow ( zl, exp ( -1.124491989777808 / 
        ( 0.4225028202459761 + zl )) ) ));
    }
  }
/*
  Calculations for Wm.
*/
  else
  {
    if ( 0.0 <= xx )
    {
      exit ( 1 );
    }
    else if ( xx <= x1 )
    {
      reta = sqrt ( d12 * delx );
      value = reta / ( reta / ( 3.0 + reta / ( reta / ( an4 
        + reta / ( reta * an6 - an5 ) ) - an3 ) ) - 1.0 ) - 1.0;
      return value;
    }
    else if ( xx <= em9 )
    {
      zl = log ( -xx );
      t = -1.0 - zl;
      ts = sqrt ( t );
      value = zl - ( 2.0 * ts ) / ( s2 + ( c13 - t 
        / ( 270.0 + ts * 127.0471381349219 )) * ts );
    }
    else
    {
      zl = log ( -xx );
      eta = 2.0 - em2 * xx;
      value = log ( xx / log ( -xx / ( ( 1.0 
        - 0.5043921323068457 * ( zl + 1.0 ) ) 
        * ( sqrt ( eta ) + eta / 3.0 ) + 1.0 )));
    }

  }

  for ( i = 1; i <= niter; i++ )
  {
    zn = log ( xx / value ) - value;
    temp = 1.0 + value;
    temp2 = temp + c23 * zn;
    temp2 = 2.0 * temp * temp2;
    value = value * ( 1.0 + ( zn / temp ) * ( temp2 - zn ) 
      / ( temp2 - 2.0 * zn ) );
  }

  return value;
}

