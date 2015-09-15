#pragma once

namespace liegroups
{
    // ab <-- a * b
    // No aliasing of input and output permitted.
    template <int N, typename S>
    void mat_mult_square(S ab[N*N], const S a[N*N], const S b[N*N]);

    // Instantiated for N=2,3
    template <int N, typename S>
#ifndef _MSC_VER
    bool invert(S invm[N*N], const S m[N*N], S *detm=0);
#else
    bool invert(S *invm, const S *m, S *detm=0);
#endif

    // A is replaced with L and U (diagonal of U is not stored).
    // Permutation is stored in index.
    // Returns true on success.
    template <int N, typename S>
    bool LU_decompose(S A[N*N], int index[N]);

    // Solves Ax=b for x, where {LU,index} is the decomposition of A.
    // x can alias b.
    // x and b have stride Stride
    template <int N, int Stride, typename S>
    void LU_inverse_times_vec(S x[N*Stride], const S LU[N*N], const int index[N], const S b[N*Stride]);
    
    // Solves AX=B for X, where {LU,index} is the decomposition of A.
    // X can alias B.
    template <int N, int C, typename S>
    void LU_inverse_times_mat(S X[N*C], const S LU[N*N], const int index[N], const S B[N*C])
    {
        for (int j=0; j<C; ++j)
            LU_inverse_times_vec<N,C>(&X[j], LU, index, &B[j]);
    }

    // Compute inv = inverse(A) from its LU decomposition.
    template <int N, typename S>
    void LU_invert(S inv[N*N], const S LU[N*N], const int index[N]);
    
    
    // s <-- sqrtm(m)
    // Aliasing permitted.
    // Returns true on success.
    template <int N, typename S>
    bool sqrtm(S s[N*N], const S m[N*N], const S tol);
    
    // em <-- expm(m)
    // Aliasing permitted.
    // Returns true on success.
    template <int N, typename S>
    bool expm(S em[N*N], const S m[N*N]);    

    // lm <-- logm(m)
    // Aliasing permitted.
    // Returns true on success.
    template <int N, typename S>
    bool logm(S lm[N*N], const S m[N*N]);    
}

