//
//  CoilSectionData.cpp
//  BlueBookFunctions
//
//  Created by Peter Huber on 2019-12-08.
//  Copyright © 2019 Peter Huber. All rights reserved.
//

#include "CoilSectionData.h"
#include "IntegrationFunctions.h"

pch_RadialData::pch_RadialData()
{
    // dummy constructor, does nothing
}

pch_RadialData::pch_RadialData(double r1, double r2, double rc, double windHt)
{
    this->ScaledC = new double[CONVERGENCE_ITERATIONS];
    this->ScaledD = new double[CONVERGENCE_ITERATIONS];
    this->ScaledE = new double[CONVERGENCE_ITERATIONS];
    this->ScaledF = new double[CONVERGENCE_ITERATIONS];
    this->PartialScaledIntL1 = new PCH_Doublet[CONVERGENCE_ITERATIONS];
    this->ScaledIntI1 = new double[CONVERGENCE_ITERATIONS];
    
    for (int n=0; n<CONVERGENCE_ITERATIONS; n++) {
        
        double useWindht = WINODW_HT_FACTOR * windHt;
        double m = n * M_PI / useWindht;
        
        double x1 = m * r1;
        double x2 = m * r2;
        double xc = m * rc;
        
        double Ri0 = gsl_sf_bessel_I0_scaled(xc);
        double Rk0 = gsl_sf_bessel_K0_scaled(xc);
        
        this->ScaledC[n] = ScaledIntegralOfTK1(x1, x2);
        // ScaledD[n-1] = Ri0 / Rk0 * ScaledC[n-1]
        this->ScaledD[n] = Ri0 / Rk0 * this->ScaledC[n];
    }
}

pch_RadialData::~pch_RadialData()
{
    delete this->ScaledC;
    delete this->ScaledD;
    delete this->ScaledE;
    delete this->ScaledF;
    delete this->PartialScaledIntL1;
    delete this->ScaledIntI1;
}
