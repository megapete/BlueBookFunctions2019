//
//  Inductance.h
//  BlueBookFunctions
//
//  Created by PeterCoolAssHuber on 2017-03-29.
//  Copyright © 2017 Peter Huber. All rights reserved.
//

#ifndef Inductance_h
#define Inductance_h

#include <stdio.h>
#include <stdbool.h>


#include "LetteredFunctions.h"
#include "IntegrationFunctions.h"

#define MU0 (4.0E-7 * M_PI)


double SelfInductance(const PCH_CoilSection aSection);
double MutualInductance(const PCH_CoilSection aFromSection, const PCH_CoilSection aToSection);

#endif /* Inductance_h */
