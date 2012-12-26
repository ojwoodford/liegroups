#pragma once

namespace liegroups
{
    // Returns max(abs(x))
    template <int N, typename S>
    S max_norm(const S x[N]);

    // ab <-- a * b
    // No aliasing of input and output permitted.
    template <int M, int C, int N, typename S>
    void mat_mult(S ab[M*N], const S a[M*C], const S b[C*N]);

    // Instantiated for N=2,3
    template <int N, typename S>
    bool invert(S invm[N*N], const S m[N*N], S *detm=0);
    
    // B gets inv(A)*B
    // A is destroyed.
    // No aliasing permitted.
    // Returns false if A is singular.
    template <int N, int C, typename S>
    bool solve(S A[N*N], S B[N*C]);

    // s <-- sqrtm(m)
    // Aliasing permitted.
    // Returns true on success.
    template <int N, typename S>
    bool sqrtm(S s[N*N], const S m[N*N]);    
    
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

