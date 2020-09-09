#include <math.h>
#include <string.h>

#define TOL 1.0e-14

double gami(double s, double x)
{
    /*
     * The incomplete lower Gamma function
     * 
     * \gamma(s, x) = \int_0^x t^(s-1) e^-t dt
     *
     * Implementation details were drawn from the Boost documentation at:
     * https://www.boost.org/doc/libs/1_37_0/libs/math/doc/sf_and_dist/html/
     *      math_toolkit/special/sf_gamma/igamma.html
     */

    // We only compute for positive arguments
    if(s <= 0 || x <= 0)
        return 0.0;

    if(x > s && x > 1.1)
    {
        // Use the continued fraction for the upper
        // incomplete gamma function, then return Gamma(s) - upper(s, x)

        double a = 1.0;
        double b = 1 + x - s;
        double Am2 = 1.0;
        double Am1 = 0.0;
        double A = 1.0;
        double Bm2 = 0.0;
        double Bm1 = 1.0;
        double B = b;
        int k = 0;
        while(fabs(Bm1*A - Am1*B) > TOL * fabs(A*Bm1))
        {
            k++;
            Am2 = Am1;
            Am1 = A;
            Bm2 = Bm1;
            Bm1 = B;
            a = k*(s-k);
            b += 2;
            A = b*Am1 + a*Am2;
            B = b*Bm1 + a*Bm2;
        }

        return tgamma(s) - pow(x, s) * exp(-x) * A / B;
    }
    else
    {
        //Use the power series for lower incomplete gamma function directly
        double term = 1.0/s;
        double sum = term;
        int k = 0;

        while(term > TOL*sum)
        {
            k++;
            term *= x/(s+k);
            sum += term;
        }

        return pow(x, s) * exp(-x) * sum;
    }
}

double dgami_(double *s, double *x)
{
    // Fortran-esque interface hack.
    return gami(*s, *x);
}

double FGABND_12plog(char *el)
{
    if(strcmp(el, "H") == 0)
        return 12.0;
    else if(strcmp(el, "He") == 0)
        return 10.99;
    else if(strcmp(el, "C") == 0)
        return 8.38;
    else if(strcmp(el, "N") == 0)
        return 7.88;
    else if(strcmp(el, "O") == 0)
        return 8.69;
    else if(strcmp(el, "Ne") == 0)
        return 7.94;
    else if(strcmp(el, "Na") == 0)
        return 6.16;
    else if(strcmp(el, "Mg") == 0)
        return 7.40;
    else if(strcmp(el, "Al") == 0)
        return 6.33;
    else if(strcmp(el, "Si") == 0)
        return 7.27;
    else if(strcmp(el, "P") == 0)
        return 5.42;
    else if(strcmp(el, "S") == 0)
        return 7.09;
    else if(strcmp(el, "Cl") == 0)
        return 5.12;
    else if(strcmp(el, "Ar") == 0)
        return 6.41;
    else if(strcmp(el, "Ca") == 0)
        return 6.20;
    else if(strcmp(el, "Ti") == 0)
        return 4.81;
    else if(strcmp(el, "Cr") == 0)
        return 5.51;
    else if(strcmp(el, "Mn") == 0)
        return 5.34;
    else if(strcmp(el, "Fe") == 0)
        return 7.43;
    else if(strcmp(el, "Co") == 0)
        return 4.92;
    else if(strcmp(el, "Ni") == 0)
        return 6.05;
    else
        return -100.0;
}

float FGABND(char *el)
{
    double logAp12 = FGABND_12plog(el);

    return (float)pow(10.0, logAp12-12);
}

