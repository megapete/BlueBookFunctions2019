//
//  CoilSectionData.hpp
//  BlueBookFunctions
//
//  Created by Peter Huber on 2019-12-08.
//  Copyright © 2019 Peter Huber. All rights reserved.
//

#ifndef CoilSectionData_hpp
#define CoilSectionData_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PCH_BasicStructs.h"

#ifndef CONVERGENCE_ITERATIONS
#define CONVERGENCE_ITERATIONS 300
#endif

#ifndef WINODW_HT_FACTOR
#define WINODW_HT_FACTOR 3.0
#endif

#ifndef MAX_COILS
#define MAX_COILS 8
#endif

typedef struct pch_ScaledDoublet {
    
    double unscaled;
    double scaled;
    
} PCH_ScaledDoublet;

typedef struct pch_RadialData {
    
    // A boolean that indicates that the struct has been properly initialized
    bool dataIsValid = false;
    
    // The arrays that hold the radial constants for all the values of n
    double *ScaledC = NULL;
    double *ScaledD = NULL;
    double *ScaledE = NULL;
    double *ScaledF = NULL;
    PCH_ScaledDoublet *PartialScaledIntL1 = NULL;
    double *ScaledIntI1 = NULL;
    
    // constructors & destructors
    pch_RadialData();
    pch_RadialData(double r1, double r2, double rc, double windHt);
    ~pch_RadialData();
    
} CoilRadialData;

class CoilSection
{
    
private:
    // an array that holds the radial constants for each coil
    static CoilRadialData *radialDataArray;
    
public:
    
    // an index into the radialDataArray
    int coilRef = -1;
    
    // section data
    double N = 0;
    double J = 0;
    double windHt = 0;
    double coreRadius = 0;
    PCH_Rect sectionRect = PCH_Rect();
    
    // constructors & destructors
    CoilSection();
    CoilSection(int coilRef, PCH_Rect sectionRect, double N, double J, double windHt, double coreRadius);
    ~CoilSection();
};

#endif /* CoilSectionData_hpp */
