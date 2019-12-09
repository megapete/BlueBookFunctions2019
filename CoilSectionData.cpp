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
        
        double useWindht = WINDOW_HT_FACTOR * windHt;
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

#pragma mark Section Attributes - Constructors & destructors

pch_SectionAttributes::pch_SectionAttributes()
{
    
}

pch_SectionAttributes::pch_SectionAttributes(std::string sectionID, int serialNumber, int inNode, int outNode)
{
    this->sectionID = sectionID;
    this->serialNumber = serialNumber;
    this->inNode = inNode;
    this->outNode = outNode;
}

pch_SectionAttributes::~pch_SectionAttributes()
{
    
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

#pragma mark CoilSection functions

double CoilSection::Jn(int n)
{
    double nPi = (double)n * M_PI;
    double useWindHt = this->windHt * WINDOW_HT_FACTOR;
    
    float z1 = this->sectionRect.origin.y;
    float z2 = z1 + this->sectionRect.size.height;
    
    return (2.0 * this->J / nPi) * (sin(nPi * z2 / useWindHt) - sin(nPi * z1 / useWindHt));
}

double CoilSection::SelfInductance()
{
    double I1 = this->J * (this->sectionRect.size.width * this->sectionRect.size.height) / this->N;
    
    double N1 = this->N;
    
    double r1 = this->sectionRect.origin.x;
    double r2 = r1 + this->sectionRect.size.width;
    double rc = this->coreRadius;
    
    double result = (M_PI * PCH_MU0 * N1 * N1 / (6.0 * WINDOW_HT_FACTOR * this->windHt)) * (gsl_pow_2(r2 + r1) + 2.0 * gsl_pow_2(r1));
    
    double multiplier = M_PI * PCH_MU0 * WINDOW_HT_FACTOR * this->windHt * N1 * N1 / gsl_pow_2(N1 * I1);
    
    int convergenceIterations = CONVERGENCE_ITERATIONS;
    
    for (int i=0; i<convergenceIterations; i++)
    {
        int n = i + 1;
        
        double m = n * M_PI / (WINDOW_HT_FACTOR * this->windHt);
        double mPow4 = m * m * m * m;
        
        double x1 = m * r1;
        double x2 = m * r2;
        double xc = m * rc;
        
        int coilIndex = this->coilRef;
        
        double scaledFn = CoilSection::radialDataArray[coilIndex].ScaledF[i];
        
        double exponentCnDn = 2.0 * (xc - x1);
        
        double scaledCn = CoilSection::radialDataArray[coilIndex].ScaledC[i];
        
        PCH_ScaledDoublet combinedL1n = CoilSection::radialDataArray[coilIndex].PartialScaledIntL1[i];
        double IntL1TermUnscaled = combinedL1n.unscaled;
        double scaledI1 = combinedL1n.scaled;
        
        double scaledEn = CoilSection::radialDataArray[coilIndex].ScaledE[i];
        
        double newWay = M_PI / 2.0 * exp(x1 - x2) * (-scaledEn) * scaledI1;
        newWay -= (M_PI / 2.0) *  IntL1TermUnscaled;
        newWay += exp(exponentCnDn) * (scaledFn * scaledCn);
        
        result += multiplier * (gsl_pow_2(this->Jn(n)) / mPow4) * newWay;
    }
    
    return result;
}

double CoilSection::MutualInductanceTo(CoilSection& otherSection)
{
    bool isSameRadialPosition = fabs(this->sectionRect.origin.x - otherSection.sectionRect.origin.x) <= 0.001;
    
    double N1 = this->N;
    double I1 = this->J * this->sectionRect.Area() / N1;
    
    double N2 = otherSection.N;
    double I2 = otherSection.J * otherSection.sectionRect.Area() / N2;
    
    double r1 = this->sectionRect.origin.x;
    double r2 = r1 + this->sectionRect.size.width;
    double r3 = otherSection.sectionRect.origin.x;
    double rc = this->coreRadius;
    
    double result = 0.0;
    
    if (isSameRadialPosition)
    {
        result = (M_PI * PCH_MU0 * N1 * N2 / (6.0 * WINDOW_HT_FACTOR * this->windHt)) * (gsl_pow_2(r2 + r1) + 2.0 * gsl_pow_2(r1));
    }
    else
    {
        result = (M_PI * PCH_MU0 * N1 * N2 / (3.0 * WINDOW_HT_FACTOR * this->windHt)) * (gsl_pow_2(r1) + r1 * r2 + gsl_pow_2(r2));
    }
    
}
