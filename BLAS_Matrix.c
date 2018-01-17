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

// Private routines
__CLPK_doublereal *AllocateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals, size_t *buffSize);

char *ElementString(BLAS_Matrix *matrix, unsigned int row, unsigned int col);

bool isVector(const BLAS_Matrix *matrix);

// Calling routines for the Set...Value should check the return value - if false, something went wrong

bool SetDoubleValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublereal value)
{
    if (theMatrix->precision != doublePrecisionMatrix)
    {
        DLog("Matrix is not double precision");
        return false;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Index out of range!");
        return false;
    }
    
    int index = -1;
    
    switch (theMatrix->type) {
            
        case vectorMatrix:
            
            index = row + col;
            break;
            
        case generalMatrix:
            
            index = col * theMatrix->numRows + row;
            break;
            
        case symmetricMatrix:
        {
            int useRow = row;
            int useCol = col;
            
            if (col < row)
            {
                useRow = col;
                useCol = row;
            }
            
            index = useRow + useCol * (useCol + 1) / 2;
            break;
        }
        case diagonalMatrix:
            if (row == col)
            {
                index = row;
            }
            break;
            
        case upperTriangularMatrix:
        case lowerTriangularMatrix:
        case bandedMatrix:
        {
            int kl = theMatrix->numSubDiags;
            int ku = theMatrix->numSuperDiags;
            
            if ((row <= col + kl) && (col <= row + ku))
            {
                index = col * (2 * kl + ku + 1) + kl + ku + row - col;
            }
            
            break;
        }
        default:
            break;
    }
    
    if (index < 0)
    {
        DLog("Illegal index for this matrix type!");
        return false;
    }
    
    theMatrix->buffer[index] = value;
    
    return true;
}

bool SetComplexValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublecomplex value)
{
    if (theMatrix->precision != complexPrecisionMatrix)
    {
        DLog("Matrix is not complex precision");
        return false;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Index out of range!");
        return false;
    }
    
    int index = -1;
    
    switch (theMatrix->type) {
            
        case vectorMatrix:
            
            index = row + col;
            break;
            
        case generalMatrix:
            
            index = col * theMatrix->numRows + row;
            break;
            
        case symmetricMatrix:
        {
            int useRow = row;
            int useCol = col;
            
            if (col < row)
            {
                useRow = col;
                useCol = row;
            }
            
            index = useRow + useCol * (useCol + 1) / 2;
            break;
        }
        case diagonalMatrix:
            if (row == col)
            {
                index = row;
            }
            break;
            
        case upperTriangularMatrix:
        case lowerTriangularMatrix:
        case bandedMatrix:
        {
            int kl = theMatrix->numSubDiags;
            int ku = theMatrix->numSuperDiags;
            
            if ((row <= col + kl) && (col <= row + ku))
            {
                index = col * (2 * kl + ku + 1) + kl + ku + row - col;
            }
            
            break;
        }
        default:
            break;
    }
    
    if (index < 0)
    {
        DLog("Illegal index for this matrix type!");
        return false;
    }
    
    __CLPK_doublecomplex *compBuff = (__CLPK_doublecomplex *)theMatrix->buffer;
    compBuff[index] = value;
    
    return true;
}

// Calling routines of the Get...Value functions should test the real value returned for the GET_VALUE_ERROR error code.

__CLPK_doublereal GetDoubleValue(const BLAS_Matrix *theMatrix, unsigned int row, unsigned int col)
{
    if (theMatrix->precision != doublePrecisionMatrix)
    {
        DLog("Matrix is not double precision");
        return GET_VALUE_ERROR;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Index out of range!");
        return GET_VALUE_ERROR;
    }
    
    int index = -1;
    
    switch (theMatrix->type) {
            
        case vectorMatrix:
            
            index = row + col;
            break;
            
        case generalMatrix:
            
            index = col * theMatrix->numRows + row;
            break;
            
        case symmetricMatrix:
        {
            int useRow = row;
            int useCol = col;
            
            if (col < row)
            {
                useRow = col;
                useCol = row;
            }
            
            index = useRow + useCol * (useCol + 1) / 2;
            break;
        }
        case diagonalMatrix:
            if (row == col)
            {
                index = row;
            }
            break;
            
        case upperTriangularMatrix:
        case lowerTriangularMatrix:
        case bandedMatrix:
        {
            int kl = theMatrix->numSubDiags;
            int ku = theMatrix->numSuperDiags;
            
            if ((row <= col + kl) && (col <= row + ku))
            {
                index = col * (2 * kl + ku + 1) + kl + ku + row - col;
            }
            
            break;
        }
        default:
            break;
    }
    
    if (index < 0)
    {
        return 0.0;
    }
    
    return theMatrix->buffer[index];
}

__CLPK_doublecomplex GetComplexValue(const BLAS_Matrix *theMatrix, unsigned int row, unsigned int col)
{
    __CLPK_doublecomplex result = {0.0, 0.0};
    
    if (theMatrix->precision != complexPrecisionMatrix)
    {
        DLog("Matrix is not complex precision");
        result.r = GET_VALUE_ERROR;
        return result;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Index out of range!");
        result.r = GET_VALUE_ERROR;
        return result;
    }
    
    int index = -1;
    
    switch (theMatrix->type) {
            
        case vectorMatrix:
            
            index = row + col;
            break;
            
        case generalMatrix:
            
            index = col * theMatrix->numRows + row;
            break;
            
        case symmetricMatrix:
        {
            int useRow = row;
            int useCol = col;
            
            if (col < row)
            {
                useRow = col;
                useCol = row;
            }
            
            index = useRow + useCol * (useCol + 1) / 2;
            break;
        }
        case diagonalMatrix:
            if (row == col)
            {
                index = row;
            }
            break;
            
        case upperTriangularMatrix:
        case lowerTriangularMatrix:
        case bandedMatrix:
        {
            int kl = theMatrix->numSubDiags;
            int ku = theMatrix->numSuperDiags;
            
            if ((row <= col + kl) && (col <= row + ku))
            {
                index = col * (2 * kl + ku + 1) + kl + ku + row - col;
            }
            
            break;
        }
        default:
            break;
    }
    
    if (index < 0)
    {
        return result;
    }
    
    __CLPK_doublecomplex *complexBuff = (__CLPK_doublecomplex *)theMatrix->buffer;
    
    result = complexBuff[index];
    
    return result;
}

BLAS_Matrix *MultiplyComplexMatrices(__CLPK_doublecomplex alpha, int transA, BLAS_Matrix *A, int transB, BLAS_Matrix *B, __CLPK_doublecomplex beta, BLAS_Matrix *C)
{
    if (A == NULL || B == NULL)
    {
        DLog("A and B must be non-NULL");
        return NULL;
    }
    
    if ((A->type != generalMatrix && A->type != vectorMatrix) || (B->type != generalMatrix && B->type != vectorMatrix) || (C != NULL && (C->type != generalMatrix && C->type != vectorMatrix)))
    {
        DLog("At this time, only general matrices or vectors can be passed to this routine");
        return NULL;
    }
    
    // do some error-checking first
    if (A->precision != B->precision)
    {
        DLog("Both matrices must be of the same type");
        return NULL;
    }
    
    if (A->precision != complexPrecisionMatrix)
    {
        DLog("This call is for Double Precision Matrices!");
        return NULL;
    }
    
    __CLPK_integer Arows = (transA == CblasNoTrans ? A->numRows : A->numCols);
    __CLPK_integer Acols = (transA == CblasNoTrans ? A->numCols : A->numRows);
    __CLPK_integer Brows = (transB == CblasNoTrans ? B->numRows : B->numCols);
    __CLPK_integer Bcols = (transB == CblasNoTrans ? B->numCols : B->numRows);
    
    if (Acols != Brows)
    {
        DLog("Incompatible A- & B-matrix dimensions!");
        return NULL;
    }
    
    if (C != NULL && (C->numRows != Arows || C->numCols != Bcols))
    {
        DLog("Incompatible C-matrix dimensions!");
        return NULL;
    }
    
    __CLPK_integer k = Acols;
    __CLPK_integer m = Arows;
    __CLPK_integer n = Bcols;
    
    __CLPK_integer lda = (transA == CblasNoTrans ? m : k);
    __CLPK_integer ldb = (transB == CblasNoTrans ? k : n);
    __CLPK_integer ldc = m;
    
    BLAS_Matrix *newC = (C != NULL ? CopyMatrix(C) : CreateMatrix(generalMatrix, complexPrecisionMatrix, ldc, n, 0, 0));
    
    if (newC == NULL)
    {
        DLog("Could not create newC matrix");
        return NULL;
    }
    
    cblas_zgemm(CblasColMajor, transA, transB, m, n, k, &alpha, A->buffer, lda, B->buffer, ldb, &beta, newC->buffer, ldc);
    
    return newC;
}

BLAS_Matrix *MultiplyDoubleMatrices(__CLPK_doublereal alpha, int transA, BLAS_Matrix *A, int transB, BLAS_Matrix *B, __CLPK_doublereal beta, BLAS_Matrix *C)
{
    if (A == NULL || B == NULL)
    {
        DLog("A and B must be non-NULL");
        return NULL;
    }
    
    if ((A->type != generalMatrix && A->type != vectorMatrix) || (B->type != generalMatrix && B->type != vectorMatrix) || (C != NULL && (C->type != generalMatrix && C->type != vectorMatrix)))
    {
        DLog("At this time, only general matrices or vectors can be passed to this routine");
        return NULL;
    }
    
    // do some error-checking first
    if (A->precision != B->precision)
    {
        DLog("Both matrices must be of the same type");
        return NULL;
    }
    
    if (A->precision != doublePrecisionMatrix)
    {
        DLog("This call is for Double Precision Matrices!");
        return NULL;
    }
    
    __CLPK_integer Arows = (transA == CblasNoTrans ? A->numRows : A->numCols);
    __CLPK_integer Acols = (transA == CblasNoTrans ? A->numCols : A->numRows);
    __CLPK_integer Brows = (transB == CblasNoTrans ? B->numRows : B->numCols);
    __CLPK_integer Bcols = (transB == CblasNoTrans ? B->numCols : B->numRows);
    
    if (Acols != Brows)
    {
        DLog("Incompatible A- & B-matrix dimensions!");
        return NULL;
    }
    
    if (C != NULL && (C->numRows != Arows || C->numCols != Bcols))
    {
        DLog("Incompatible C-matrix dimensions!");
        return NULL;
    }
    
    __CLPK_integer k = Acols;
    __CLPK_integer m = Arows;
    __CLPK_integer n = Bcols;
    
    __CLPK_integer lda = (transA == CblasNoTrans ? m : k);
    __CLPK_integer ldb = (transB == CblasNoTrans ? k : n);
    __CLPK_integer ldc = m;
    
    BLAS_Matrix *newC = (C != NULL ? CopyMatrix(C) : CreateMatrix(generalMatrix, doublePrecisionMatrix, ldc, n, 0, 0));
    
    if (newC == NULL)
    {
        DLog("Could not create newC matrix");
        return NULL;
    }
    
    cblas_dgemm(CblasColMajor, transA, transB, m, n, k, alpha, A->buffer, lda, B->buffer, ldb, beta, newC->buffer, ldc);
    
    return newC;
}

BLAS_Matrix *TransposeMatrix(const BLAS_Matrix *srcMatrix)
{
    // Diagonal matrices and vectors (as far as this library is concerned) are the same as their transposes.
    if (srcMatrix->type == diagonalMatrix || isVector(srcMatrix))
    {
        return CopyMatrix(srcMatrix);
    }
    
    BLAS_Matrix *newMatrix = CreateMatrix(srcMatrix->type, srcMatrix->precision, srcMatrix->numCols, srcMatrix->numRows, srcMatrix->numSuperDiags, srcMatrix->numSubDiags);
    
    for (int i=0; i<srcMatrix->numCols; i++)
    {
        for (int j=0; j<srcMatrix->numRows; j++)
        {
            if (srcMatrix->precision == doublePrecisionMatrix)
            {
                __CLPK_doublereal value = GetDoubleValue(srcMatrix, j, i);
                SetDoubleValue(newMatrix, i, j, value);
            }
            else if (srcMatrix->precision == complexPrecisionMatrix)
            {
                __CLPK_doublecomplex value = GetComplexValue(srcMatrix, j, i);
                SetComplexValue(newMatrix, i, j, value);
            }
            else
            {
                DLog("Illegal precision specified!");
                free(newMatrix);
                return NULL;
            }
        }
    }
    
    return newMatrix;
}

bool isVector(const BLAS_Matrix *matrix)
{
    return (matrix->numRows == 1 || matrix->numCols == 1);
}

BLAS_Matrix *CopyMatrix(const BLAS_Matrix *srcMatrix)
{
    if (srcMatrix == NULL)
    {
        DLog("Cannot copy NULL matrix!");
        return NULL;
    }
    
    size_t buffSize = 0;
    
    void *buff = AllocateMatrix(srcMatrix->type, srcMatrix->precision, srcMatrix->numRows, srcMatrix->numCols, srcMatrix->numSubDiags, srcMatrix->numSuperDiags, &buffSize);
    
    if (buff == NULL)
    {
        DLog("Could not allocate buffer for matrix!");
        return NULL;
    }
    
    
    memcpy(buff, srcMatrix->buffer, buffSize);
    
    BLAS_Matrix *newMatrix = malloc(sizeof(BLAS_Matrix));
    
    if (newMatrix == NULL)
    {
        DLog("Could not allocate memory for matrix struct!");
        free(buff);
        return NULL;
    }
    
    newMatrix->type = srcMatrix->type;
    newMatrix->precision = srcMatrix->precision;
    newMatrix->numRows = srcMatrix->numRows;
    newMatrix->numCols = srcMatrix->numCols;
    newMatrix->numSubDiags = srcMatrix->numSubDiags;
    newMatrix->numSuperDiags = srcMatrix->numSuperDiags;
    newMatrix->buffer = buff;
    newMatrix->bufferSize = buffSize;
    
    return newMatrix;
}

BLAS_Matrix *CreateVector(MatrixPrecision precision, unsigned int numElements)
{
    return CreateMatrix(vectorMatrix, precision, numElements, 0, 0, 0);
}

BLAS_Matrix *CreateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals)
{
    // Do a bunch of error-checking to start
    if (rows == 0 || columns == 0)
    {
        DLog("Either or both of rows and columns are 0");
        return NULL;
    }
    
    if (type != vectorMatrix && (rows == 1 || columns == 1))
    {
        DLog("Only vectors can be specified with 1 row or 1 column! Fixing!");
        
        if (rows == 1)
        {
            rows = columns;
            columns = 1;
        }
        
        type = vectorMatrix;
        
        subDiagonals = 0;
        superDiagonals = 0;
    }
    
    if (type == vectorMatrix && rows > 1 && columns > 1)
    {
        // I don't think this code can ever be executed
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
    
    size_t buffSize = 0;
    
    void *buff = AllocateMatrix(type , precision, rows, columns, subDiagonals, superDiagonals, &buffSize);
    
    if (buff == NULL)
    {
        free(matrix);
        DLog("Error creating the matrix buffer!");
        return NULL;
    }
    
    matrix->buffer = buff;
    matrix->bufferSize = buffSize;
    
    return matrix;
}

// Allocate memory for a matrix of the given type & precision
__CLPK_doublereal *AllocateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals, size_t *buffSize)
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
    
    size_t numElements = 0;
    // Take care of the different matrix types
    if (type == vectorMatrix)
    {
        // Vectors only have a single row (or column). We will take the higher of rows or columns as the dimension of the vector.
        numElements = (rows > columns ? rows : columns);
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
        *buffSize = 0;
    }
    else
    {
        *buffSize = numElements * elementSize;
    }
    
    return result;
}

char *MatrixAsString(BLAS_Matrix *theMatrix)
{
    char buffer[32567] = "";
    
    for (int j=0; j<theMatrix->numRows; j++)
    {
        strcat(buffer, "[ ");
        for (int i=0; i<theMatrix->numCols; i++)
        {
            if (theMatrix->precision == doublePrecisionMatrix)
            {
                __CLPK_doublereal value = GetDoubleValue(theMatrix, j, i);
                char tempString[11];
                snprintf(tempString, 11, "%10.3G", value);
                strcat(buffer, tempString);
            }
            else if (theMatrix->precision == complexPrecisionMatrix)
            {
                __CLPK_doublecomplex value = GetComplexValue(theMatrix, j, i);
                char tempString[21];
                snprintf(tempString, 21, "%10.3G%+10.3G", value.r, value.i);
            }
            else
            {
                DLog("Illegal precision specification");
                return NULL;
            }
        }
        strcat(buffer, " ]\n");
    }
    
    char *result = malloc(strlen(buffer) + 1);
    strcpy(result, buffer);
    return result;
}

// Private function to convert the element at the given row and column into a string. The calling function is responsible for freeing the returned pointer.
char *ElementString(BLAS_Matrix *matrix, unsigned int row, unsigned int col)
{
    // allow space for a ridiculously large string
    char tempString[256];
    
    if (matrix->precision == doublePrecisionMatrix)
    {
        __CLPK_doublereal value = GetDoubleValue(matrix, row, col);
        snprintf(tempString, 256, "%.3G", value);
    }
    else if (matrix->precision == complexPrecisionMatrix)
    {
        __CLPK_doublecomplex value = GetComplexValue(matrix, row, col);
        
        snprintf(tempString, 256, "%.3G%+.3G", value.r, value.i);
    }
    else
    {
        DLog("Illegal precision specification");
        return NULL;
    }
    
    char *result = malloc(strlen(tempString) + 1);
    strcpy(result, tempString);
    return result;
}
