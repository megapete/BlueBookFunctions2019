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

// The "stopping" point for the integrals
#define BBF_RELATIVE_ERROR (1.0E-4)

// Public function declarations
double M0(double x);
double M1(double x);
double IntegralOfM0From0to(double b);
double IntegralOfM0(double a, double b);

#endif /* IntegrationFunctions_h */
