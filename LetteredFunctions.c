//
//  LetteredFunctions.c
//  BlueBookFunctions
//
//  Created by PeterCoolAssHuber on 2017-03-27.
//  Copyright © 2017 Peter Huber. All rights reserved.
//
// These functions are basically ported from Swift code - there may be more comments in the original code.

#include "LetteredFunctions.h"

double J0(const PCH_CoilSection theSection)
{
    double sectionHt = fabs(theSection.z2 - theSection.z1);
    
    return theSection.J * sectionHt / (theSection.windowHt * BBF_WINDHT_FACTOR);
}

double J(int n, const PCH_CoilSection theSection)
{
    double nPi = (double)n * M_PI;
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    
    return (2.0 * theSection.J / nPi) * (sin(nPi * (theSection.z2) / useWindHt) - sin(nPi * theSection.z1 / useWindHt));
}

double C(int n, const PCH_CoilSection theSection)
{
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    double m = (double)n * M_PI / useWindHt;
    
    double x1 = m * theSection.r1;
    double x2 = m * theSection.r2;
    
    return IntegralOfTK1(x1, x2);
}

double D(int n, const PCH_CoilSection theSection)
{
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    
    double xc = ((double)n * M_PI / useWindHt) * theSection.coreRadius;
    
    double Ri0 = gsl_sf_bessel_I0_scaled(xc);
    double Rk0 = gsl_sf_bessel_K0_scaled(xc);
    double eBase = exp(2.0 * xc);
    
    double result = eBase * (Ri0 / Rk0) * C(n, theSection);
    
    return result;
}

double ScaledD(int n, const PCH_CoilSection theSection)
{
    // returns Rd where D = exp(2.0 * xc - x1) * Rd (xc and x1 are functions of n)
    
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    double m = (double)n * M_PI / useWindHt;
    
    double x1 = m * theSection.r1;
    double x2 = m * theSection.r2;
    double xc = m * theSection.coreRadius;
    
    double Ri0 = gsl_sf_bessel_I0_scaled(xc);
    double Rk0 = gsl_sf_bessel_K0_scaled(xc);
    
    double scaledCn = ScaledIntegralOfTK1(x1, x2);
    
    return Ri0 / Rk0 * scaledCn;
}

double AlternateD(int n, const PCH_CoilSection theSection)
{
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    double m = (double)n * M_PI / useWindHt;
    
    double x1 = m * theSection.r1;
    double x2 = m * theSection.r2;
    double xc = m * theSection.coreRadius;
    
    double Ri0 = gsl_sf_bessel_I0_scaled(xc);
    double Rk0 = gsl_sf_bessel_K0_scaled(xc);
    
    double scaledCn = ScaledIntegralOfTK1(x1, x2);
    
    return exp(2.0 * xc - x1) * Ri0 / Rk0 * scaledCn;
}

double E(int n, const PCH_CoilSection theSection)
{
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    double m = (double)n * M_PI / useWindHt;
    
    double x2 = m * theSection.r2;
    
    return IntegralOfTK1From0to(x2);
}

double F(int n, const PCH_CoilSection theSection)
{
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    double m = (double)n * M_PI / useWindHt;
    
    double x1 = m * theSection.r1;
    
    double result = AlternateD(n, theSection) - IntegralOfTI1From0to(x1);
    
    return result;
}

double ScaledF(int n, const PCH_CoilSection theSection)
{
    // return Rf where F = exp(2.0 * xc - x1) * Rf (xc and x1 are functions of n)
    
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    double m = (double)n * M_PI / useWindHt;
    
    double x1 = m * theSection.r1;
    double xc = m * theSection.coreRadius;
    
    double exponent = 2.0 * xc - x1;
    
    double result = ScaledD(n, theSection) - exp(x1 - exponent) * ScaledIntegralOfTI1From0to(x1);
    
    return result;
}

double AlternateF(int n, const PCH_CoilSection theSection)
{
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    double m = (double)n * M_PI / useWindHt;
    
    double x1 = m * theSection.r1;
    double xc = m * theSection.coreRadius;
    
    double exponent = 2.0 * xc - x1;
    
    double result = exp(exponent) * (ScaledD(n, theSection) - exp(x1 - exponent) * ScaledIntegralOfTI1From0to(x1));
    
    return result;
}

double G(int n, const PCH_CoilSection theSection)
{
    double useWindHt = theSection.windowHt * BBF_WINDHT_FACTOR;
    double m = (double)n * M_PI / useWindHt;
    
    double x1 = m * theSection.r1;
    double x2 = m * theSection.r2;
    double xc = m * theSection.coreRadius;
    
    double Ri0 = gsl_sf_bessel_I0_scaled(xc);
    double Rk0 = gsl_sf_bessel_K0_scaled(xc);
    double eBase = exp(2.0 * xc);
    
    double result = eBase * (Ri0 / Rk0) * IntegralOfTK1(x1, x2) + IntegralOfTI1(x1, x2);
    
    return result;
}
