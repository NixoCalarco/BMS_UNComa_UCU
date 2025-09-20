#include "emf.h"
#include "lut_emf.h"
#include "lut_emf_df.h"

static inline double linear_interp(double x, double x0, double x1, double y0,
                                   double y1) {
    if (x0 == x1)
        return y0;

    return y0 + ((x - x0) / (x1 - x0)) * (y1 - y0);
}

static inline double binary_search(double x, double lut_x[static EMF_N],
                                   double lut_y[static EMF_N]) {
    if (x < lut_x[0])
        return emf_y[0];

    if (x > lut_x[EMF_N - 1])
        return lut_y[EMF_N - 1];

    int low = 0, high = EMF_N - 1;
    while (high - low > 1) {
        int mid = (low + high) / 2;
        if (x < lut_x[mid])
            high = mid;
        else
            low = mid;
    }

    return linear_interp(x, lut_x[low], lut_x[high], lut_y[low], lut_y[high]);
}

double emf(double x) { return binary_search(x, emf_x, emf_y); }

double inverse_emf(double y) { return binary_search(y, emf_y, emf_x); }

double derivative_emf(double x) { return binary_search(x, emf_df_x, emf_df_y); }
