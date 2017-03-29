//
//  LetteredFunctions.h
//  BlueBookFunctions
//
//  Created by PeterCoolAssHuber on 2017-03-27.
//  Copyright © 2017 Peter Huber. All rights reserved.
//

#ifndef LetteredFunctions_h
#define LetteredFunctions_h

#include <stdio.h>

#include "IntegrationFunctions.h"

// The functions depend on a height factor. The Book says that a value of 3 works well, so we'll use that unless tetsing shows we should use something else.
#define BBF_WINDHT_FACTOR (3.0)

// All of the lettered functions depend on certain data that is associated with the coil section that we are considering. We define a struct for that data
typedef struct section_tag
{
    double r1; // innermost radius (m)
    double r2; // outermost radius (m)
    
    double z1; // dimension closest to bottom yoke (m)
    double z2; // dimension closest to top yoke (m)
    
    double N; // number of turns in the section
    double J; // current density in A/m2
    
    double coreRadius; // (m)
    double windowHt; // (m)
    
} PCH_CoilSection;

// Function declarations
double J0(const PCH_CoilSection theSection);
double J(int n, const PCH_CoilSection theSection);
double C(int n, const PCH_CoilSection theSection);
double D(int n, const PCH_CoilSection theSection);
double E(int n, const PCH_CoilSection theSection);
double F(int n, const PCH_CoilSection theSection);
double G(int n, const PCH_CoilSection theSection);

// Scaled versions of select functions
double ScaledD(int n, const PCH_CoilSection theSection);
double ScaledF(int n, const PCH_CoilSection theSection);

// Alternate versions of select functions
double AlternateD(int n, const PCH_CoilSection theSection);
double AlternateF(int n, const PCH_CoilSection theSection);

#endif /* LetteredFunctions_h */
