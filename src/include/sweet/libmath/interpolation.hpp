/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#ifndef SWEET_INTERPOLATE_HPP
#define SWEET_INTERPOLATE_HPP


/**
 * Compute Lagrange interpolation for given interpolation points and interpolant values
 * 
 * https://de.wikipedia.org/wiki/Polynominterpolation#Lagrangesche_Interpolationsformel
 *
 * Here we need to provide the nonequidistantly coordinates of the interpolation points.
 */
template <int N>
double interpolation_lagrange_nonequidistant(
    double *x,    /// interpolation points
    double *y,    /// interpolation values
    double x_sample /// sample position
)
{
    double retval = 0;
    for (int i = 0; i < N; i++)
    {
        double denom = 1;
        double nom = 1;

        for (int j = 0; j < N; j++)
        {
            if (i == j)
                continue;

            nom *= x_sample - x[j];
            denom *= x[i] - x[j];
        }

        retval += y[i]*nom/denom;
    }

    return retval;
}




/**
 * Compute Lagrange interpolation for given interpolation points and interpolant values
 *
 * https://de.wikipedia.org/wiki/Polynominterpolation#Lagrangesche_Interpolationsformel
 *
 * Here we assume the equidistantly spaced points starting at 0
 */
template <int N>
double interpolation_lagrange_equidistant(
    double *y,    /// interpolation values
    double x_sample /// sample position
)
{
    double retval = 0;
    for (int i = 0; i < N; i++)
    {
        double denom = 1;
        double nom = 1;

        for (int j = 0; j < N; j++)
        {
            if (i == j)
                continue;

            nom *= x_sample - (double)j;
            denom *= (double)(i - j);
        }

        retval += y[i]*nom/denom;
    }

    return retval;
}

#endif
