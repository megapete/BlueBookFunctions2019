//
//  CoilSectionData.cpp
//  BlueBookFunctions
//
//  Created by Peter Huber on 2019-12-08.
//  Copyright © 2019 Peter Huber. All rights reserved.
//

#include "CoilSectionData.hpp"
#include "IntegrationFunctions.h"

#include "PCH_C_Logging.h"

#pragma mark CoilRadialData - Radial data field initializer
static CoilRadialData InitialRadialData[] = {CoilRadialData(), CoilRadialData(), CoilRadialData(), CoilRadialData(), CoilRadialData(), CoilRadialData(), CoilRadialData(), CoilRadialData()};

CoilRadialData *CoilSection::radialDataArray = InitialRadialData;

#pragma mark CoilRadialData - Constructors
CoilRadialData::pch_RadialData()
{
    // dummy constructor, does nothing
    
}

CoilRadialData::pch_RadialData(double r1, double r2, double rc, double windHt)
{
    this->ScaledC = new double[CONVERGENCE_ITERATIONS];
    this->ScaledD = new double[CONVERGENCE_ITERATIONS];
    this->ScaledE = new double[CONVERGENCE_ITERATIONS];
    this->ScaledF = new double[CONVERGENCE_ITERATIONS];
    this->PartialScaledIntL1 = new PCH_ScaledDoublet[CONVERGENCE_ITERATIONS];
    this->ScaledIntI1 = new double[CONVERGENCE_ITERATIONS];
    
    for (int n=0; n<CONVERGENCE_ITERATIONS; n++) {
        
        double useWindht = WINODW_HT_FACTOR * windHt;
        double m = (n + 1) * M_PI / useWindht;
        
        double x1 = m * r1;
        double x2 = m * r2;
        double xc = m * rc;
        
        double Ri0 = gsl_sf_bessel_I0_scaled(xc);
        double Rk0 = gsl_sf_bessel_K0_scaled(xc);
        
        // ScaledC[n-1] = ScaledIntegralOf_tK1_from(x1, toB: x2)
        this->ScaledC[n] = ScaledIntegralOfTK1(x1, x2);
        
        // ScaledD[n-1] = Ri0 / Rk0 * ScaledC[n-1]
        this->ScaledD[n] = Ri0 / Rk0 * this->ScaledC[n];
        
        // ScaledE[n-1] = ScaledIntegralOf_tK1_from0_to(x2)
        this->ScaledE[n] = ScaledIntegralOfTK1From0to(x2);
        
        // ScaledF[n-1] = (ScaledD[n-1] - exp(2.0 * (x1 - xc)) * ScaledIntegralOf_tI1_from0_to(x1))
        this->ScaledF[n] = this->ScaledD[n] - exp(2.0 * (x1 - xc)) * ScaledIntegralOfTI1From0to(x1);
        
        // PartialScaledIntL1[n-1] = PartialScaledIntegralOf_tL1_from(x1, toB: x2)
        double scaledL1Result = 0.0;
        double unscaledL1Result = PartialScaledIntegralOfTL1(x1, x2, &scaledL1Result);
        PCH_ScaledDoublet combinedResult;
        combinedResult.unscaled = unscaledL1Result;
        combinedResult.scaled = scaledL1Result;
        this->PartialScaledIntL1[n] = combinedResult;
        
        // ScaledIntI1[n-1] = ScaledIntegralOf_tI1_from(x1, toB: x2)
        this->ScaledIntI1[n] = ScaledIntegralOfTI1(x1, x2);
    }
    
    this->dataIsValid = true;
}

#pragma mark CoilRadialData - Destructor
pch_RadialData::~pch_RadialData()
{
    if (this->dataIsValid == false)
    {
        return;
    }
    
    delete[] this->ScaledC;
    delete[] this->ScaledD;
    delete[] this->ScaledE;
    delete[] this->ScaledF;
    delete[] this->PartialScaledIntL1;
    delete[] this->ScaledIntI1;
}

#pragma mark CoilSection - Constructors
CoilSection::CoilSection()
{
    // dummy constructor, does nothing
}

CoilSection::CoilSection(int coilRef, PCH_Rect sectionRect, double N, double J, double windHt, double coreRadius)
{
    if ((coilRef < 0) || (coilRef >= MAX_COILS))
    {
        DLog("Tried to pass an illegal number as a coil ref - returning!");
        return;
    }
    
    this->coilRef = coilRef;
    this->sectionRect = sectionRect;
    this->N = N;
    this->J = J;
    this->windHt = windHt;
    this->coreRadius = coreRadius;
    
    // check if the radial data entry for this coil ref has been created and if not, create it
    if (!CoilSection::radialDataArray[coilRef].dataIsValid)
    {
        CoilSection::radialDataArray[coilRef] = CoilRadialData(sectionRect.origin.x, sectionRect.origin.x + sectionRect.size.width, coreRadius, windHt);
    }
    
}
