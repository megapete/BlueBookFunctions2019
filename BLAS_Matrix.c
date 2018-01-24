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

// Private static variable to initialize a SolveResult struct to 0
// static const struct _tagSolveResult emptySolveResult = {0};

// Calling routines for the Set...Value should check the return value - if false, something went wrong

MatrixStatus SetDoubleValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublereal value)
{
    if (theMatrix->precision != doublePrecisionMatrix)
    {
        DLog("Matrix is not double precision");
        return PRECISION_MISMATCH_MATRIX_ERROR;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Index out of range!");
        return ILLEGAL_INDEX_MATRIX_ERROR;
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
        DLog("Unreachable index for this matrix type! Ignoring!");
        return UNREACHABLE_INDEX_MATRIX_ERROR;
    }
    
    theMatrix->buffer[index] = value;
    
    return NO_MATRIX_ERROR;
}

MatrixStatus SetComplexValue(BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, __CLPK_doublecomplex value)
{
    if (theMatrix->precision != complexPrecisionMatrix)
    {
        DLog("Matrix is not complex precision");
        return PRECISION_MISMATCH_MATRIX_ERROR;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Index out of range!");
        return ILLEGAL_INDEX_MATRIX_ERROR;
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
        return UNREACHABLE_INDEX_MATRIX_ERROR;
    }
    
    __CLPK_doublecomplex *compBuff = (__CLPK_doublecomplex *)theMatrix->buffer;
    compBuff[index] = value;
    
    return NO_MATRIX_ERROR;
}

// Calling routines of the Get...Value functions should test the real value returned for the GET_VALUE_ERROR error code.

__CLPK_doublereal GetDoubleValue(const BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, MatrixStatus *status)
{
    if (status != NULL)
    {
        *status = NO_MATRIX_ERROR;
    }
    
    if (theMatrix->precision != doublePrecisionMatrix)
    {
        DLog("Matrix is not double precision");
        
        if (status != NULL)
        {
            *status = PRECISION_MISMATCH_MATRIX_ERROR;
        }
        
        return GET_VALUE_ERROR;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Index out of range!");
        
        if (status != NULL)
        {
            *status = ILLEGAL_INDEX_MATRIX_ERROR;
        }
        
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
        case positiveDefiniteMatrix:
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

__CLPK_doublecomplex GetComplexValue(const BLAS_Matrix *theMatrix, unsigned int row, unsigned int col, MatrixStatus *status)
{
    if (status != NULL)
    {
        *status = NO_MATRIX_ERROR;
    }
    
    __CLPK_doublecomplex result = {0.0, 0.0};
    
    if (theMatrix->precision != complexPrecisionMatrix)
    {
        DLog("Matrix is not complex precision");
        
        if (status != NULL)
        {
            *status = PRECISION_MISMATCH_MATRIX_ERROR;
        }
        
        result.r = GET_VALUE_ERROR;
        return result;
    }
    
    if (row >= theMatrix->numRows || col >= theMatrix->numCols)
    {
        DLog("Index out of range!");
        
        if (status != NULL)
        {
            *status = ILLEGAL_INDEX_MATRIX_ERROR;
        }
        
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

void ReleaseSolveResult(SolveResult sResult)
{
    free(sResult.A);
    free(sResult.AF);
    free(sResult.B);
    free(sResult.bErr);
    free(sResult.C);
    free(sResult.fErr);
    free(sResult.R);
    free(sResult.X);
}


SolveResult Solve(bool useExtendedMethod, BLAS_Matrix *A, BLAS_Matrix *B, BLAS_Factoring factor, BLAS_Transposing transpose, BLAS_Equilibrating equilibrate)
{
    SolveResult result = {0};
    
    if (factor == BLAS_Factored)
    {
        DLog("The use of pre-factored matrices is NOT yet implemented!");
        result.status = ILLEGAL_FACTORING_SOLVE_ERROR;
        return result;
    }
    
    if (A->numCols != B->numRows)
    {
        DLog("Matrix dimensions are not compatible!");
        result.status = DIMENSION_MISMATCH_MATRIX_ERROR;
        return result;
    }
    
    if (B->type != generalMatrix)
    {
        DLog("B matrix must be a general matrix");
        result.status = ILLEGAL_TYPE_MATRIX_ERROR;
        return result;
    }
    
    // Set up some variables that are common to all the routines
    result.A = malloc(A->bufferSize);
    result.B = malloc(B->bufferSize);
    memcpy(result.A, A->buffer, A->bufferSize);
    memcpy(result.B, B->buffer, B->bufferSize);
    
    __CLPK_integer n = A->numRows;
    result.lda = n;
    result.ldb = n;
    __CLPK_integer nrhs = B->numCols;
    result.ipiv = malloc(n);
    
    // kl and ku are for triangular matrices but we'll create and initialize them even though other matrix types will not use them
    __CLPK_integer kl = A->numSubDiags;
    __CLPK_integer ku = A->numSuperDiags;
    
    if (!useExtendedMethod)
    {
        // Simple driver routines
        
        if (A->precision == doublePrecisionMatrix)
        {
            if (A->type == generalMatrix)
            {
                dgesv_(&n, &nrhs, result.A, &result.lda, result.ipiv, result.B, &result.ldb, &result.info);
            }
            else if (A->type == bandedMatrix || A->type == upperTriangularMatrix || A->type == lowerTriangularMatrix || A->type == diagonalMatrix)
            {
                __CLPK_integer ldab = 2 * kl + ku + 1;
                dgbsv_(&n, &kl, &ku, &nrhs,result.A, &ldab, result.ipiv, result.B, &result.ldb, &result.info);
            }
            else if (A->type == symmetricMatrix)
            {
                char uplo[2] = "U";
                __CLPK_doublereal workSize;
                __CLPK_integer lWork = -1;
                
                dsysv_(uplo, &n, &nrhs, result.A, &result.lda, result.ipiv, result.B, &result.ldb, &workSize, &lWork, &result.info);
                
                lWork = (__CLPK_integer)workSize;
                __CLPK_doublereal work[lWork];
                
                dsysv_(uplo, &n, &nrhs, result.A, &result.lda, result.ipiv, result.B, &result.ldb, work, &lWork, &result.info);
            }
            else if (A->type == positiveDefiniteMatrix)
            {
                char uplo[2] = "U";
                
                dposv_(uplo, &n, &nrhs, result.A, &result.lda, result.B, &result.ldb, &result.info);
            }
            else
            {
                DLog("Illegal matrix type!");
                result.status = ILLEGAL_TYPE_MATRIX_ERROR;
                return result;
            }
        }
        else if (A->precision == complexPrecisionMatrix)
        {
            if (A->type == generalMatrix)
            {
                zgesv_(&n, &nrhs, result.A, &result.lda, result.ipiv, result.B, &result.ldb, &result.info);
            }
            else if (A->type == bandedMatrix || A->type == upperTriangularMatrix || A->type == lowerTriangularMatrix || A->type == diagonalMatrix)
            {
                __CLPK_integer ldab = 2 * kl + ku + 1;
                zgbsv_(&n, &kl, &ku, &nrhs,result.A, &ldab, result.ipiv, result.B, &result.ldb, &result.info);
            }
            else if (A->type == symmetricMatrix)
            {
                char uplo[2] = "U";
                __CLPK_doublecomplex workSize;
                __CLPK_integer lWork = -1;
                
                zsysv_(uplo, &n, &nrhs, result.A, &result.lda, result.ipiv, result.B, &result.ldb, &workSize, &lWork, &result.info);
                
                lWork = (__CLPK_integer)workSize.r;
                __CLPK_doublecomplex work[lWork];
                
                zsysv_(uplo, &n, &nrhs, result.A, &result.lda, result.ipiv, result.B, &result.ldb, work, &lWork, &result.info);
            }
            else if (A->type == positiveDefiniteMatrix)
            {
                char uplo[2] = "U";
                
                zposv_(uplo, &n, &nrhs, result.A, &result.lda, result.B, &result.ldb, &result.info);
            }
            
            else
            {
                DLog("Illegal matrix type!");
                result.status = ILLEGAL_TYPE_MATRIX_ERROR;
                return result;
            }
        }
        else
        {
            DLog("Unimplemented precision!");
            result.status = ILLEGAL_PRECISION_MATRIX_ERROR;
            return result;
        }
    }
    else // routines with extended in/out arguments (not sure what they all do)
    {
        if (A->precision == doublePrecisionMatrix)
        {
            if (A->type == generalMatrix)
            {
                // dgesv_(&n, &nrhs, result.A, &result.lda, result.ipiv, result.B, &result.ldb, &result.info);
                char fact = 'N';
                if (factor == BLAS_Equilibrate)
                {
                    fact = 'E';
                }
                
                char trans = 'N';
                if (transpose == BLAS_Transpose)
                {
                    trans = 'T';
                }
                else if (transpose == BLAS_ConjugateTranspose)
                {
                    trans = 'C';
                }
                
                result.AF = malloc(A->bufferSize);
                result.ldaf = n;
                result.equed = 'N';
                result.R = calloc(n, sizeof(__CLPK_doublereal));
                result.C = calloc(n, sizeof(__CLPK_doublereal));
                result.X = malloc(B->bufferSize);
                result.ldx = n;
                result.fErr = calloc(nrhs, sizeof(__CLPK_doublereal));
                result.bErr = calloc(nrhs, sizeof(__CLPK_doublereal));
                result.work = calloc(4 * n, sizeof(__CLPK_doublereal));
                result.iWork = calloc(n, sizeof(__CLPK_integer));
                
                dgesvx_(&fact, &trans, &n, &nrhs, result.A, &result.lda, result.AF, &result.ldaf, result.ipiv, &result.equed, result.R, result.C, result.B, &result.ldb, result.X, &result.ldx, &result.rCond, result.fErr, result.bErr, result.work, result.iWork, &result.info);
            }
            else
            {
                DLog("Unimplemented matrix type!");
                result.status = ILLEGAL_TYPE_MATRIX_ERROR;
                return result;
            }
        }
        else if (A->precision == complexPrecisionMatrix)
        {
            if (A->type == generalMatrix)
            {
                // dgesv_(&n, &nrhs, result.A, &result.lda, result.ipiv, result.B, &result.ldb, &result.info);
                char fact = 'N';
                if (factor == BLAS_Equilibrate)
                {
                    fact = 'E';
                }
                
                char trans = 'N';
                if (transpose == BLAS_Transpose)
                {
                    trans = 'T';
                }
                else if (transpose == BLAS_ConjugateTranspose)
                {
                    trans = 'C';
                }
                
                result.AF = malloc(A->bufferSize);
                result.ldaf = n;
                result.equed = 'N';
                result.R = calloc(n, sizeof(__CLPK_doublereal));
                result.C = calloc(n, sizeof(__CLPK_doublereal));
                result.X = malloc(B->bufferSize);
                result.ldx = n;
                result.fErr = calloc(nrhs, sizeof(__CLPK_doublereal));
                result.bErr = calloc(nrhs, sizeof(__CLPK_doublereal));
                result.work = calloc(2 * n, sizeof(__CLPK_doublecomplex));
                result.rWork = calloc(2 * n, sizeof(__CLPK_doublereal));
                
                zgesvx_(&fact, &trans, &n, &nrhs, result.A, &result.lda, result.AF, &result.ldaf, result.ipiv, &result.equed, result.R, result.C, result.B, &result.ldb, result.X, &result.ldx, &result.rCond, result.fErr, result.bErr, result.work, result.rWork, &result.info);
            }
            else
            {
                DLog("Unimplemented matrix type!");
                result.status = ILLEGAL_TYPE_MATRIX_ERROR;
                return result;
            }
        }
        else
        {
            DLog("Unimplemented precision!");
            result.status = ILLEGAL_PRECISION_MATRIX_ERROR;
            return result;
        }
    }
    
    if (result.info != 0)
    {
        result.status = SOLVE_MATRIX_ERROR;
    }
    
    return result;
}

BLAS_Matrix *SolveSimpleLinearSystem(BLAS_Matrix *A, BLAS_Matrix *B)
{
    SolveResult result = Solve(false, A, B, BLAS_NoFactor, BLAS_NoTranspose, BLAS_NoEquilibration);
    
    BLAS_Matrix *newMatrix = CreateMatrix(B->type, B->precision, B->numRows, B->numCols, 0, 0);
    
    if (result.status == NO_MATRIX_ERROR)
    {
        memcpy(newMatrix->buffer, result.X, newMatrix->bufferSize);
    }
    else
    {
        newMatrix->status = result.status;
        newMatrix->lapack_info = result.info;
    }
    
    ReleaseSolveResult(result);
    
    return newMatrix;
}

BLAS_Matrix *Inverse(const BLAS_Matrix *srcMatrix)
{
    if (srcMatrix == NULL)
    {
        DLog("Attempt to invert NULL matrix!");
        return NULL;
    }
    
    BLAS_Matrix *A = CopyMatrix(srcMatrix);
    
    if (srcMatrix->numRows != srcMatrix->numCols)
    {
        DLog("Only square matrices can be inverted");
        A->status = EXPECTED_SQUARE_MATRIX_ERROR;
        return A;
    }
    
    if (A->status != NO_MATRIX_ERROR)
    {
        DLog("Could not copy source matrix");
        return A;
    }
    
    __CLPK_integer m = A->numRows;
    __CLPK_integer n = A->numCols;
    __CLPK_integer lda = m;
    
    int ipivCount = (m < n ? m : n);
    __CLPK_integer ipiv[ipivCount];
    
    __CLPK_integer info = 0;
    
    if (A->precision == doublePrecisionMatrix)
    {
        dgetrf_(&m, &n, A->buffer, &lda, ipiv, &info);
        
        if (info != 0)
        {
            DLog("An error occurred in dgetrf_(): %d", info);
            A->status = LAPACK_ROUTINE_ERROR;
            A->lapack_info = info;
            
            return A;
        }
        
        __CLPK_doublereal workSize;
        __CLPK_integer lWork = -1;
        
        dgetri_(&n, A->buffer, &lda, ipiv, &workSize, &lWork, &info);
        
        lWork = (__CLPK_integer)workSize;
        
        __CLPK_doublereal work[lWork];
        
        dgetri_(&n, A->buffer, &lda, ipiv, work, &lWork, &info);
        
        if (info != 0)
        {
            DLog("An error occurred in dgetri_(): %d", info);
            A->status = LAPACK_ROUTINE_ERROR;
            A->lapack_info = info;
            
            return A;
        }
    }
    else
    {
        zgetrf_(&m, &n, (__CLPK_doublecomplex *)(A->buffer), &lda, ipiv, &info);
        
        if (info != 0)
        {
            DLog("An error occurred in zgetrf_(): %d", info);
            A->status = LAPACK_ROUTINE_ERROR;
            A->lapack_info = info;
            
            return A;
        }
        
        __CLPK_doublecomplex workSize;
        __CLPK_integer lWork = -1;
        
        zgetri_(&n, (__CLPK_doublecomplex *)(A->buffer), &lda, ipiv, &workSize, &lWork, &info);
        
        lWork = (__CLPK_integer)workSize.r;
        
        __CLPK_doublecomplex work[lWork];
        
        zgetri_(&n, (__CLPK_doublecomplex *)(A->buffer), &lda, ipiv, work, &lWork, &info);
        
        if (info != 0)
        {
            DLog("An error occurred in zgetri_(): %d", info);
            A->status = LAPACK_ROUTINE_ERROR;
            A->lapack_info = info;
            
            return A;
        }
    }
    
    return A;
}

BLAS_Matrix *MultiplyComplexMatrices(__CLPK_doublecomplex alpha, int transA, BLAS_Matrix *A, int transB, BLAS_Matrix *B, __CLPK_doublecomplex beta, BLAS_Matrix *C)
{
    if (A == NULL || B == NULL)
    {
        DLog("A and B must be non-NULL");
        return NULL;
    }
    
    BLAS_Matrix *dummyMatrix = CopyMatrix(A);
    
    if ((A->type != generalMatrix && A->type != vectorMatrix) || (B->type != generalMatrix && B->type != vectorMatrix) || (C != NULL && (C->type != generalMatrix && C->type != vectorMatrix)))
    {
        DLog("At this time, only general matrices or vectors can be passed to this routine");
        dummyMatrix->status = ILLEGAL_TYPE_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    // do some error-checking first
    if (A->precision != B->precision)
    {
        DLog("Both matrices must be of the same type");
        dummyMatrix->status = PRECISION_MISMATCH_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    if (A->precision != complexPrecisionMatrix)
    {
        DLog("This call is for Complex Precision Matrices!");
        dummyMatrix->status = ILLEGAL_PRECISION_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    __CLPK_integer Arows = (transA == CblasNoTrans ? A->numRows : A->numCols);
    __CLPK_integer Acols = (transA == CblasNoTrans ? A->numCols : A->numRows);
    __CLPK_integer Brows = (transB == CblasNoTrans ? B->numRows : B->numCols);
    __CLPK_integer Bcols = (transB == CblasNoTrans ? B->numCols : B->numRows);
    
    if (Acols != Brows)
    {
        DLog("Incompatible A- & B-matrix dimensions!");
        dummyMatrix->status = DIMENSION_MISMATCH_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    if (C != NULL && (C->numRows != Arows || C->numCols != Bcols))
    {
        dummyMatrix->status = DIMENSION_MISMATCH_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    // if we get here, we no longer need dummyMatrix
    ReleaseMatrix(dummyMatrix);
    
    __CLPK_integer k = Acols;
    __CLPK_integer m = Arows;
    __CLPK_integer n = Bcols;
    
    __CLPK_integer lda = (transA == CblasNoTrans ? m : k);
    __CLPK_integer ldb = (transB == CblasNoTrans ? k : n);
    __CLPK_integer ldc = m;
    
    BLAS_Matrix *newC = (C != NULL ? CopyMatrix(C) : CreateMatrix(generalMatrix, complexPrecisionMatrix, ldc, n, 0, 0));
    
    if (newC->status != NO_MATRIX_ERROR)
    {
        DLog("Could not create newC matrix");
        return newC;
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
    
    BLAS_Matrix *dummyMatrix = CopyMatrix(A);
    
    if ((A->type != generalMatrix && A->type != vectorMatrix) || (B->type != generalMatrix && B->type != vectorMatrix) || (C != NULL && (C->type != generalMatrix && C->type != vectorMatrix)))
    {
        DLog("At this time, only general matrices or vectors can be passed to this routine");
        dummyMatrix->status = ILLEGAL_TYPE_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    // do some error-checking first
    if (A->precision != B->precision)
    {
        DLog("Both matrices must be of the same type");
        dummyMatrix->status = PRECISION_MISMATCH_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    if (A->precision != doublePrecisionMatrix)
    {
        DLog("This call is for Double Precision Matrices!");
        dummyMatrix->status = ILLEGAL_PRECISION_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    __CLPK_integer Arows = (transA == CblasNoTrans ? A->numRows : A->numCols);
    __CLPK_integer Acols = (transA == CblasNoTrans ? A->numCols : A->numRows);
    __CLPK_integer Brows = (transB == CblasNoTrans ? B->numRows : B->numCols);
    __CLPK_integer Bcols = (transB == CblasNoTrans ? B->numCols : B->numRows);
    
    if (Acols != Brows)
    {
        DLog("Incompatible A- & B-matrix dimensions!");
        dummyMatrix->status = DIMENSION_MISMATCH_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    if (C != NULL && (C->numRows != Arows || C->numCols != Bcols))
    {
        dummyMatrix->status = DIMENSION_MISMATCH_MATRIX_ERROR;
        return dummyMatrix;
    }
    
    // if we get here, we no longer need dummyMatrix
    ReleaseMatrix(dummyMatrix);
    
    __CLPK_integer k = Acols;
    __CLPK_integer m = Arows;
    __CLPK_integer n = Bcols;
    
    __CLPK_integer lda = (transA == CblasNoTrans ? m : k);
    __CLPK_integer ldb = (transB == CblasNoTrans ? k : n);
    __CLPK_integer ldc = m;
    
    BLAS_Matrix *newC = (C != NULL ? CopyMatrix(C) : CreateMatrix(generalMatrix, doublePrecisionMatrix, ldc, n, 0, 0));
    
    if (newC->status != NO_MATRIX_ERROR)
    {
        DLog("Could not create newC matrix");
        return newC;
    }
    
    cblas_dgemm(CblasColMajor, transA, transB, m, n, k, alpha, A->buffer, lda, B->buffer, ldb, beta, newC->buffer, ldc);
    
    return newC;
}

BLAS_Matrix *TransposeMatrix(const BLAS_Matrix *srcMatrix)
{
    // This is a pretty inefficient function. Generally, it would be better to to use the "transpose" flag in the BLAS and LAPACK routines instead of actually transposing the matrix (those routines just switch up the indexes instead of actually transposing any matrices).
    
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
                __CLPK_doublereal value = GetDoubleValue(srcMatrix, j, i, NULL);
                SetDoubleValue(newMatrix, i, j, value);
            }
            else if (srcMatrix->precision == complexPrecisionMatrix)
            {
                __CLPK_doublecomplex value = GetComplexValue(srcMatrix, j, i, NULL);
                SetComplexValue(newMatrix, i, j, value);
            }
            else
            {
                DLog("Illegal precision specified!");
                newMatrix->status = ILLEGAL_PRECISION_MATRIX_ERROR;
                return newMatrix;
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
        DLog("Attempt to copy NULL matrix!");
        return NULL;
    }
    
    BLAS_Matrix *newMatrix = CreateMatrix(srcMatrix->type, srcMatrix->precision, srcMatrix->numRows, srcMatrix->numCols, srcMatrix->numSubDiags, srcMatrix->numSuperDiags);
    
    if (newMatrix->status != NO_MATRIX_ERROR)
    {
        return newMatrix;
    }
    
    if (newMatrix->bufferSize != srcMatrix->bufferSize)
    {
        DLog("Buffer sizes do not match!");
        newMatrix->status = BUFFERSIZE_MATRIX_ERROR;
        return newMatrix;
    }
    
    memcpy(newMatrix->buffer, srcMatrix->buffer, srcMatrix->bufferSize);
    
    return newMatrix;
}

void ReleaseMatrix(BLAS_Matrix *matrix)
{
    if (matrix != NULL)
    {
        free(matrix->buffer);
        free(matrix);
    }
}

BLAS_Matrix *CreateVector(MatrixPrecision precision, unsigned int numElements)
{
    return CreateMatrix(vectorMatrix, precision, numElements, 0, 0, 0);
}

BLAS_Matrix *CreateMatrix(MatrixType type, MatrixPrecision precision,  unsigned int rows, unsigned int columns, unsigned int subDiagonals, unsigned int superDiagonals)
{
    BLAS_Matrix *matrix = malloc(sizeof(BLAS_Matrix));
    if (matrix == NULL)
    {
        ALog("Could not allocate memory for matrix struct!");
        return NULL;
    }
    
    // We set the buffer to NULL to ensure that any access crashes if there's an error
    matrix->buffer = NULL;
    matrix->status = NO_MATRIX_ERROR;
    
    // Do a bunch of error-checking to start
    if (rows == 0 || columns == 0)
    {
        DLog("Either or both of rows and columns are 0");
        matrix->status = ZERO_DIMENSION_MATRIX_ERROR;
        return matrix;
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
        matrix->status = VECTOR_DIMENSION_MATRIX_ERROR;
        return matrix;
    }
    
    if (type == symmetricMatrix || type == lowerTriangularMatrix || type == upperTriangularMatrix || type == hermitianMatrix)
    {
        if (rows != columns)
        {
            DLog("Expected square dimensions for this type of matrix!");
            matrix->status = EXPECTED_SQUARE_MATRIX_ERROR;
            return matrix;
        }
    }
    
    if (type == hermitianMatrix && precision != complexPrecisionMatrix)
    {
        DLog("Hermitian matrices must be complex precision!");
        matrix->status = NONCOMPLEX_HERMITIAN_MATRIX_ERROR;
        return matrix;
    }
    
    if (type == bandedMatrix)
    {
        // I think banded matrices need to be square (see https://en.wikipedia.org/wiki/Band_matrix, in the first sentence under the heading "Bandwidth", which says: "Formally, consider an n×n matrix...".
        
        if (rows != columns)
        {
            DLog("Expected square dimensions for this type of matrix!");
            matrix->status = EXPECTED_SQUARE_MATRIX_ERROR;
            return matrix;
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
    else if (type == upperTriangularMatrix)
    {
        subDiagonals = 0;
        superDiagonals = columns - 1;
    }
    else if (type == lowerTriangularMatrix)
    {
        subDiagonals = columns - 1;
        superDiagonals = 0;
    }
    else if (type == diagonalMatrix)
    {
        subDiagonals = 0;
        superDiagonals = 0;
    }
    
    // If we get here, we were successful
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
        DLog("Error creating the matrix buffer!");
        matrix->status = BUFFER_ALLOCATION_MATRIX_ERROR;
        return matrix;
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
    else if (type == generalMatrix || type == symmetricMatrix || type == positiveDefiniteMatrix || type == hermitianMatrix)
    {
        // the solving routines do not yield any space-saving for symmeteric (or positive definite) matrices
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
                __CLPK_doublereal value = GetDoubleValue(theMatrix, j, i, NULL);
                char tempString[11];
                snprintf(tempString, 11, "%10.3G", value);
                strcat(buffer, tempString);
            }
            else if (theMatrix->precision == complexPrecisionMatrix)
            {
                __CLPK_doublecomplex value = GetComplexValue(theMatrix, j, i, NULL);
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
        __CLPK_doublereal value = GetDoubleValue(matrix, row, col, NULL);
        snprintf(tempString, 256, "%.3G", value);
    }
    else if (matrix->precision == complexPrecisionMatrix)
    {
        __CLPK_doublecomplex value = GetComplexValue(matrix, row, col, NULL);
        
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
