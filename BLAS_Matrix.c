//
//  BLAS_Matrix.c
//  BlueBookFunctions
//
//  Created by Peter Huber on 2018-01-11.
//  Copyright © 2018 Peter Huber. All rights reserved.
//

#include "BLAS_Matrix.h"
#include "PCH_C_Logging.h"

// Allocate memory for a double precision matrix of the given type
__CLPK_doublereal *AllocateMatrix(MatrixType type, MatrixPrecision precision, int rows, int columns)
{
    void *result = NULL;
    
    size_t elementSize = 0;
    
    if (precision == doublePrecisionMatrix)
    {
        elementSize = sizeof(__CLPK_doublereal);
    }
    else if (precision == complexPrecisionMatrix)
    {
        elementSize = sizeof(__CLPK_doublecomplex);
    }
    else
    {
        
    }
    
    return result;
}
