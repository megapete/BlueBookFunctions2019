//
//  IntegrationFunctions.h
//  BlueBookFunctions
//
//  Created by PeterCoolAssHuber on 2017-03-27.
//  Copyright © 2017 Peter Huber. All rights reserved.
//

#ifndef IntegrationFunctions_h
#define IntegrationFunctions_h

#include <stdio.h>

#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

// The "stopping" point for the integrals. This should be as small as possible without making calculations take too long
#define BBF_RELATIVE_ERROR (1.0E-8)
#define BBF_ABSOLUTE_ERROR (1.0E-12)

// Public function declarations
double M0(double x);
double M1(double x);

double IntegralOfM0From0to(double b);
double IntegralOfM0(double a, double b);

double IntegralOfTK1From0to(double b);
double IntegralOfTK1(double a, double b);

double IntegralOfTI1From0to(double b);
double IntegralOfTI1(double a, double b);

double IntegralOfTL1From0to(double b);
double IntegralOfTL1(double a, double b);

// Scaled versions of selected functions (used for accuracy and performance)

double ScaledIntegralOfTK1(double a, double b);
double ScaledIntegralOfTK1From0to(double b);
double ScaledIntegralOfTI1From0to(double b);
double ScaledIntegralOfTI1(double a, double b);

// Special, partially scaled function. The returned value is the unscaled part of the result and the parameter "scaledResult" is scaled.
double PartialScaledIntegralOfTL1(double a, double b, double *scaledResult);

// Alternate versions of selected functions (used for accuracy and performance)
double AlternateIntegralOfTL1(double a, double b);

#endif /* IntegrationFunctions_h */
