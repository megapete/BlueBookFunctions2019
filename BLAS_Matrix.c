//
//  BLAS_Matrix.c
//  BlueBookFunctions
//
//  Created by Peter Huber on 2018-01-11.
//  Copyright © 2018 Peter Huber. All rights reserved.
//

// NOTE: Most of the logic in this class comes from my Objective-C-based (and subsequent Swift-based) BLAS-matrix classes. See those projects for more notes on the indexing and storage methods.

#include "BLAS_Matrix.h"
#include "PCH_C_Logging.h"

// Calling routines of the Get...Value functions should test the real value returned for the GET_VALUE_ERROR error code.

__CLPK_doublereal GetDoubleValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col)
{
    if (theMatrix->precision != doublePrecisionMatrix)
    {
        DLog("Matrix is not double precision");
        return GET_VALUE_ERROR;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Illegal index!");
        return GET_VALUE_ERROR;
    }
    
    int index = -1;
    
    switch (theMatrix->type) {
            
        case vectorMatrix:
            
            index = row + col;
            break;
            
        case generalMatrix:
            
            index = col * theMatrix->numRows + row
            break;
            
        case symmetricMatrix:
            
            int useRow = row;
            int useCol = col;
            
            if (col < row)
            {
                useRow = col;
                useCol = row;
            }
            
            index = useRow + useCol * (useCol + 1) / 2];
            break;
            
        case diagonalMatrix:
            if (row == col)
            {
                index = row;
            }
            break;
            
        case upperTriangularMatrix:
        case lowerTriangularMatrix:
        case bandedMatrix:
            
            kl = theMatrix->numSubDiags;
            ku = theMatrix->numSuperDiags;
            
            if (row <= col + kl) && (col <= row + ku)
            {
                index = col * (2 * kl + ku + 1) + kl + ku + row - col
            }
            
            break;
            
        default:
            break;
    }
    
    if (index < 0)
    {
        return 0.0;
    }
    
    return theMatrix->buffer[index];
}

BLAS_Matrix *CreateVector(MatrixType type, MatrixPrecision precision, unsigned int numElements)
{
    return CreateMatrix(type, precision, 0, numElements, 0, 0);
}

BLAS_Matrix *CreateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals)
{
    // Do a bunch of error-checking to start
    if (rows == 0 && columns == 0)
    {
        DLog("Both rows and columns are 0");
        return NULL;
    }
    
    if (type != vectorMatrix && (rows == 0 || columns == 0))
    {
        DLog("Only vectors can be specified with 0 rows or 0 columns!");
        return NULL;
    }
    
    if (type == vectorMatrix && rows > 0 && columns > 0)
    {
        DLog("Vectors must be defined with at least one of 'rows' or 'columns' defined as 0");
        return NULL;
    }
    
    if (type == symmetricMatrix || type == lowerTriangularMatrix || type == upperTriangularMatrix)
    {
        if (rows != columns)
        {
            DLog("Expected square dimensions for this type of matrix!");
            return NULL;
        }
    }
    
    if (type == bandedMatrix)
    {
        // I think banded matrices need to be square (see https://en.wikipedia.org/wiki/Band_matrix, in the first sentence under the heading "Bandwidth", which says: "Formally, consider an n×n matrix...".
        
        if (rows != columns)
        {
            DLog("Expected square dimensions for this type of matrix!");
            return NULL;
        }
        
        // If the sub- or superDiagonals arguments are greater than the maximum allowed, we set it to the maximum.
        if (subDiagonals >= columns)
        {
            DLog("Specified subDiagonals too high - setting to max (consider using lowerTriangularMatrix instead");
            subDiagonals = columns - 1;
        }
        
        if (superDiagonals >= columns)
        {
            DLog("Specified subDiagonals too high - setting to max (consider using upperTriangularMatrix instead");
            superDiagonals = columns - 1;
        }
    }
    
    BLAS_Matrix *matrix = malloc(sizeof(BLAS_Matrix));
    
    matrix->type = type;
    matrix->precision = precision;
    
    matrix->numRows = rows;
    matrix->numCols = columns;
    
    if (type == vectorMatrix && columns == 0)
    {
        // We default vectors to being "single-row". Manipulation routines will have to consider this.
        matrix->numRows = 0;
        matrix->numCols = rows;
    }
    
    if (type == upperTriangularMatrix)
    {
        matrix->numSuperDiags = columns - 1;
        matrix->numSubDiags = 0;
    }
    else if (type == lowerTriangularMatrix)
    {
        matrix->numSuperDiags = 0;
        matrix->numSubDiags = columns - 1;
    }
    else if (type == bandedMatrix)
    {
        matrix->numSubDiags = subDiagonals;
        matrix->numSuperDiags = superDiagonals;
    }
    else
    {
        matrix->numSubDiags = 0;
        matrix->numSuperDiags = 0;
    }
    
    void *buff = AllocateMatrix(type , precision, rows, columns, subDiagonals, superDiagonals);
    
    if (buff == NULL)
    {
        free(matrix);
        DLog("Error creating the matrix buffer!");
        return NULL;
    }
    
    matrix->buffer = buff;
    
    return matrix;
}

// Allocate memory for a matrix of the given type & precision
__CLPK_doublereal *AllocateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals)
{
    void *result = NULL;
    
    size_t elementSize = 0;
    
    // Take care of precision first
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
        DLog("Illegal precision specified!");
        return NULL;
    }
    
    size_t numElements = 0
    // Take care of the different matrix types
    if (type == vectorMatrix)
    {
        // Vectors only have a single row (or column). We will take the higher of rows or columns as the dimension of the vector.
        numElements = (rows > columns ? rows : columns)
    }
    else if (type == generalMatrix)
    {
        numElements = rows * columns;
    }
    else if (type == lowerTriangularMatrix || type == upperTriangularMatrix)
    {
        unsigned int diagonals = rows - 1;
        numElements = columns * (2 * diagonals + 1);
    }
    else if (type == diagonalMatrix)
    {
        // see https://en.wikipedia.org/wiki/Main_diagonal
        numElements = (rows < columns ? rows : columns);
    }
    else if (type == bandedMatrix)
    {
        numElements = columns * (2 * subDiagonals + superDiagonals + 1);
    }
    else
    {
        DLog("Unimplemented matrix type specified!");
        return NULL;
    }
    
    if (numElements == 0)
    {
        DLog("Not able to allocate 0 bytes");
        return NULL;
    }
    
    result = calloc(numElements, elementSize);
    
    if (result == NULL)
    {
        DLog("Could not allocate memory");
    }
    
    return result;
}
