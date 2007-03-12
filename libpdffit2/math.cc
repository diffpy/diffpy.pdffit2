/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Math functions used in formula parser.
*
* Comments:
*
* $Id$
*
***********************************************************************/

#include "pdffit.h"

double sind(double arg)
{
    return sin(arg*rad);
}

double cosd(double arg)
{
    return cos(arg*rad);
}

double tand(double arg)
{
    return tan(arg*rad);
}

double asind(double arg)
{
    return asin(arg)/rad;
}

double acosd(double arg)
{
    return acos(arg)/rad;
}

double atand(double arg)
{
    return atan(arg)/rad;
}

double neg(double arg)
{
    return -arg;
}

double dsin(double x)
{
    return cos(x);
}

double dcos(double x)
{
    return -sin(x);
}

double dtan(double x)
{
    return 1/sqr(cos(x));
}

double dsind(double x)
{
    return rad*cosd(x);
}

double dcosd(double x)
{
    return -rad*sind(x);
}

double dtand(double x)
{
    return rad/sqr(cosd(x));

}

double dasin(double x)
{
    return 1/sqrt(1-sqr(x));
}

double dacos(double x)
{
    return -1/sqrt(1-sqr(x));
}


double datan(double x)
{
    return 1/(1+sqr(x));
}

double dasind(double x)
{
    return dasin(x)/rad;
}

double dacosd(double x)
{
    return dacos(x)/rad;
}

double datand(double x)
{
    return datan(x)/rad;
}

double dexp(double x)
{
    return exp(x);
}

double dlog(double x)
{
    return 1/x;
}

double dlog10(double x)
{
    return 1/(x*M_LN10);
}

double dsqr(double x)
{
    return 2.0*x;
}

double dcube(double x)
{
    return 3.0*sqr(x);
}

double dsqrt(double x)
{
    return 1/(2.0*sqrt(x));
}

double dneg(double arg)
{
    return -1.0;
}

// End of file
