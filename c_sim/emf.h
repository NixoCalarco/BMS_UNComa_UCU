#ifndef EMF_H
#define EMF_H

#define N 1000 // cantidad de valores en cada LUT

static inline double linear_interp(double x, double x0, double x1, double y0,
                                   double y1);

static inline double binary_search(double x, double lut_x[static N],
                                   double lut_y[static N]);

double emf(double x);

double inverse_emf(double y);

double derivative_emf(double x);

#endif
