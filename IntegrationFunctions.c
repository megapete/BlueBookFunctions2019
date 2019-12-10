//
//  IntegrationFunctions.c
//  BlueBookFunctions
//
//  Created by PeterCoolAssHuber on 2017-03-27.
//  Copyright © 2017 Peter Huber. All rights reserved.
//
// These functions are basically ported from Swift code - there may be more comments in the original code.

#include "IntegrationFunctions.h"

double PartialScaledIntegralOfTL1(double a, double b, double *scaledResult)
{
    // This function splits up the result into two parts, unscaled and scaled (to exp(a)). The calling routine should take advantage (if possible) of this to increase accuracy of the result.
    
    double nonIntegralSum = a * M0(a) - b * M0(b) + (a*a - b*b) / M_PI;
    double m0integral = IntegralOfM0(a, b);
    // double scaledI1integral = ScaledIntegralOfTI1(a, b);
    
    if (scaledResult != NULL)
    {
        *scaledResult = ScaledIntegralOfTI1(a, b);
    }
    
    return nonIntegralSum + m0integral;
}

double AlternateIntegralOfTL1(double a, double b)
{
    // More "mathematical" way of calculating the integral
    double nonIntegralSum = a * M0(a) - b * M0(b) + (a*a - b*b) / M_PI;
    double m0integral = IntegralOfM0(a, b);
    double i1integral = exp(a) * ScaledIntegralOfTI1(a, b);
    
    return nonIntegralSum + m0integral + i1integral;
}

double IntegralOfTL1From0to(double b)
{
    return (-b * M0(b)) - (b * b / M_PI) + IntegralOfM0From0to(b) + IntegralOfTI1From0to(b);
}

double IntegralOfTL1(double a, double b)
{
    return IntegralOfTL1From0to(b) - IntegralOfTL1From0to(a);
}

double ScaledIntegralOfTI1From0to(double b)
{
    double Ri0 = gsl_sf_bessel_I0_scaled(b);
    double Ri1 = gsl_sf_bessel_I1_scaled(b);
    
    return M_PI_2 * b * (M1(b) * Ri0 - M0(b) * Ri1);
}

double ScaledIntegralOfTI1(double a, double b)
{
    double Ri0a = gsl_sf_bessel_I0_scaled(a);
    double Ri1a = gsl_sf_bessel_I1_scaled(a);
    double Ri0b = gsl_sf_bessel_I0_scaled(b);
    double Ri1b = gsl_sf_bessel_I1_scaled(b);
    
    double firstTerm = a * (M1(a) * Ri0a - M0(a) * Ri1a);
    double secondTerm = b * exp(b-a) * (M1(b) * Ri0b - M0(b) * Ri1b);
    
    return M_PI_2 * (secondTerm - firstTerm);
}

double IntegralOfTI1From0to(double b)
{
    double Ri0 = gsl_sf_bessel_I0_scaled(b);
    double Ri1 = gsl_sf_bessel_I1_scaled(b);
    double eBase = exp(b);
    
    return M_PI_2 * b * eBase * (M1(b) * Ri0 - M0(b) * Ri1);
}

double IntegralOfTI1(double a, double b)
{
    return IntegralOfTI1From0to(b) - IntegralOfTI1From0to(a);
}


double IntegralOfTK1From0to(double b)
{
    double Rk0 = gsl_sf_bessel_K0_scaled(b);
    double Rk1 = gsl_sf_bessel_K1_scaled(b);
    double eBase = exp(-b);
    
    double secondTerm = b * eBase * (M1(b) * Rk0 + M0(b) * Rk1);
    
    double result = M_PI_2 * (1.0 - secondTerm);
    
    return result;
}

double IntegralOfTK1(double a, double b)
{
    return IntegralOfTK1From0to(b) - IntegralOfTK1From0to(a);
}

double ScaledIntegralOfTK1From0to(double b)
{
    // return IntK1 where the actual integral = π / 2.0 * (1.0 - exp(-b) * IntK1)
    double Rk0 = gsl_sf_bessel_K0_scaled(b);
    double Rk1 = gsl_sf_bessel_K1_scaled(b);
    
    double result = b * (M1(b) * Rk0 + M0(b) * Rk1);
    
    return result;
}

double ScaledIntegralOfTK1(double a, double b)
{
    if (a == 0)
    {
        return ScaledIntegralOfTK1From0to(b);
    }
    
    // We assume that in the event a>b, the caller simply mixed the two up
    if (a > b)
    {
        fprintf(stderr, "Inverted bounds of integration in ScaledIntegralOfTK1. Fixing!");
        
        double newB = a;
        a = b;
        b = newB;
    }
    
    double Rk0a = gsl_sf_bessel_K0_scaled(a);
    double Rk1a = gsl_sf_bessel_K1_scaled(a);
    double Rk0b = gsl_sf_bessel_K0_scaled(b);
    double Rk1b = gsl_sf_bessel_K1_scaled(b);
    
    double firstTerm = a * (M1(a) * Rk0a + M0(a) * Rk1a);
    double secondTerm = b * exp(a-b) * (M1(b) * Rk0b + M0(b) * Rk1b);
    
    double result = M_PI_2 * (firstTerm - secondTerm);
    
    return result;
}

// Private function needed for the integral in IntegralOfM0From0to()
double IntM0_integrand(double theta, void *param)
{
    double *x = param;
    return (1.0 - exp(-(*x) * cos(theta))) / cos(theta);
}

double IntegralOfM0From0to(double b)
{
    double err = 0.0;
    size_t numEvals = 0;
    double result = 0.0;
    
    gsl_function F;
    double params = b;
    F.function = &IntM0_integrand;
    F.params = &params;
    
    int integrationResult = gsl_integration_qng(&F, 0.0, M_PI_2, 0.0, BBF_RELATIVE_ERROR, &result, &err, &numEvals);
    
    if (integrationResult > 0)
    {
        fprintf(stderr, "Error in IntM0 integration call: %d", integrationResult);
        result = 0.0;
    }
    
    return result * 2.0 / M_PI;
}

double IntegralOfM0(double a, double b)
{
    return IntegralOfM0From0to(b) - IntegralOfM0From0to(a);
}

// Private function needed for the integral in M1()
double M1X_integrand(double theta, void *param)
{
    double *x = param;
    return exp(-(*x) * cos(theta)) * cos(theta);
}

double M1(double x)
{
    double err = 0.0;
    size_t numEvals = 0;
    double result = 0.0;
    
    gsl_function F;
    double params = x;
    F.function = &M1X_integrand;
    F.params = &params;
    
    int integrationResult = gsl_integration_qng(&F, 0.0, M_PI_2, 0.0, BBF_RELATIVE_ERROR, &result, &err, &numEvals);
    
    if (integrationResult > 0)
    {
        fprintf(stderr, "Error in M1 integration call: %d", integrationResult);
        result = 0.0;
    }
    
    return (1.0 - result) * 2.0 / M_PI;
}

// Private function needed for the integral in M0()
double M0X_integrand(double theta, void *param)
{
    double *x = param;
    return exp(-(*x) * cos(theta));
}

double M0(double x)
{
    double err = 0.0;
    size_t numEvals = 0;
    double result = 0.0;
    
    gsl_function F;
    double params = x;
    F.function = &M0X_integrand;
    F.params = &params;
    
    int integrationResult = gsl_integration_qng(&F, 0.0, M_PI_2, 0.0, BBF_RELATIVE_ERROR, &result, &err, &numEvals);
    
    if (integrationResult > 0)
    {
        fprintf(stderr, "Error in M0 integration call: %d", integrationResult);
        result = 0.0;
    }
    
    return result * 2.0 / M_PI;
}
