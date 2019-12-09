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
#include <string>
#include <map>
#include "PCH_BasicStructs.hpp"

#ifndef CONVERGENCE_ITERATIONS
#define CONVERGENCE_ITERATIONS 300
#endif

#ifndef WINDOW_HT_FACTOR
#define WINDOW_HT_FACTOR 3.0
#endif

#ifndef MAX_COILS
#define MAX_COILS 8
#endif

#ifndef PCH_MU0
#define PCH_MU0 (M_PI * 4.0E-7)
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

typedef struct pch_SectionAttributes {
    
    std::string sectionID = "";
    
    int serialNumber = -1;
    int inNode;
    int outNode;
    
    double seriesCapacitance;
    std::map<std::string, double> shuntCapacitances;
    
    double resistance;
    
    double selfInductance;
    std::map<std::string, double> mutualInductances;
    
    // constructors & destructors
    pch_SectionAttributes();
    pch_SectionAttributes(std::string sectionID, int serialNumber, int inNode, int outNode);
    ~pch_SectionAttributes();
    
} PCH_SectionAttributes;

class CoilSection
{
    
private:
    // an array that holds the radial constants for each coil (class property)
    static CoilRadialData *radialDataArray;
    
public:
    
    // an index into the radialDataArray
    int coilRef = -1;
    
    // immutable section data
    double N = 0;
    double J = 0;
    double windHt = 0;
    double coreRadius = 0;
    PCH_Rect sectionRect = PCH_Rect();
    
    PCH_SectionAttributes attributes;
    
    // constructors & destructors
    CoilSection();
    CoilSection(int coilRef, PCH_Rect sectionRect, double N, double J, double windHt, double coreRadius);
    ~CoilSection();
    
    // Inductance calculators & associated functions
    double Jn(int n);
    double SelfInductance();
};

#endif /* CoilSectionData_hpp */
