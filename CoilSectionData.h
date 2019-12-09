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

typedef struct pch_Doublet {
    
    double d1;
    double d2;
    
} PCH_Doublet;

typedef struct pch_RadialData {
    
    double *ScaledC = NULL;
    double *ScaledD= NULL;
    double *ScaledE= NULL;
    double *ScaledF= NULL;
    PCH_Doublet *PartialScaledIntL1= NULL;
    double *ScaledIntI1= NULL;
    
    pch_RadialData();
    pch_RadialData(double r1, double r2, double rc, double windHt);
    ~pch_RadialData();
    
} CoilRadialData;

class CoilSection
{
    static CoilRadialData radialDataArray[MAX_COILS];
};

#endif /* CoilSectionData_hpp */
