//
//  IntegrationFunctions.c
//  BlueBookFunctions
//
//  Created by PeterCoolAssHuber on 2017-03-27.
//  Copyright © 2017 Peter Huber. All rights reserved.
//

#include "IntegrationFunctions.h"

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
