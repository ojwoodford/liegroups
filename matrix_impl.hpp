#pragma once

#include <liegroups/matrix.hpp>
#include <liegroups/scalar.hpp>
#include <iostream>

template <typename T>
static void swap(T &a, T &b)
{
    T tmp = a;
    a = b;
    b = tmp;
}

template <typename S>
static void swap(S *a, S *b, int n)
{
    for (int i=0; i<n; ++i) {
        S tmp = a[i];
        a[i] = b[i];
        b[i] = tmp;
    }
}

template <int N, typename S>
static void copy(S out[N], const S in[N])
{
    for (int i=0; i<N; ++i)
        out[i] = in[i];
}

template <int N, typename S>
S liegroups::max_norm(const S x[N])
{
    S mv = 0;
    for (int i=0; i<N; ++i) {
        S axi = liegroups::abs(x[i]);
        mv = liegroups::max(axi, mv);
    }
    return mv;
}

namespace liegroups
{
    template <int M, int C, int N>
    struct MatMult
    {
        template <typename S>
        static void eval(S ab[M*N], const S a[M*C], const S b[C*N])
        {
            for (int i=0; i<M; ++i) {
                for (int j=0; j<N; ++j) {
                    S sum = (S)0;
                    for (int k=0; k<C; ++k) {
                        sum += a[i*C + k] * b[j + k*N];
                    }
                    ab[i*N + j] = sum;
                }
            }
        }        
    };

    template <>
    struct MatMult<3,3,3>
    {
        template <typename S>
        static void eval(S ab[3*3], const S a[3*3], const S b[3*3])
        {
            ab[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
            ab[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
            ab[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];

            ab[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
            ab[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
            ab[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];

            ab[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
            ab[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
            ab[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
        }        
    };

}
        
template <int M, int C, int N, typename S>
void liegroups::mat_mult(S ab[M*N], const S a[M*C], const S b[C*N])
{
    MatMult<M,C,N>::eval(ab, a, b);
}

template <typename S>
static bool invert2(S invm[2*2], const S m[2*2], S *detm)
{
    const S det = m[0]*m[3] - m[1]*m[2];
    if (detm)
        *detm = det;
    if (det == (S)0)
        return false;

    const S inv_det = (S)1 / det;
    invm[0] =  inv_det * m[3];
    invm[1] = -inv_det * m[1];
    invm[2] = -inv_det * m[2];
    invm[3] =  inv_det * m[0];
    return true;
}

template <typename S>
static bool invert3(S invm[3*3], const S m[3*3], S *detm)
{
    const S m00 = m[4]*m[8] - m[5]*m[7];
    const S m10 = m[1]*m[8] - m[2]*m[7];
    const S m20 = m[1]*m[5] - m[2]*m[4];

    const S det = m[0]*m00 - m[3]*m10 + m[6]*m20;
    if (detm)
        *detm = det;
    if (det == (S)0)
        return false;

    const S inv_det = (S)1 / det;
    const S m01 = m[3]*m[8] - m[5]*m[6];
    const S m02 = m[3]*m[7] - m[4]*m[6];
    const S m11 = m[0]*m[8] - m[2]*m[6];
    const S m12 = m[0]*m[7] - m[1]*m[6];
    const S m21 = m[0]*m[5] - m[2]*m[3];
    const S m22 = m[0]*m[4] - m[1]*m[3];

    invm[0] =  inv_det * m00;
    invm[1] = -inv_det * m10;
    invm[2] =  inv_det * m20;
    invm[3] = -inv_det * m01;
    invm[4] =  inv_det * m11;
    invm[5] = -inv_det * m21;
    invm[6] =  inv_det * m02;
    invm[7] = -inv_det * m12;
    invm[8] =  inv_det * m22;
    return true;
}

namespace liegroups {

    template <>
    bool invert<2,float>(float invm[2*2], const float m[2*2], float *detm)
    {
        return invert2(invm, m, detm);
    }

    template <>
    bool invert<2,double>(double invm[2*2], const double m[2*2], double *detm)
    {
        return invert2(invm, m, detm);
    }
    
    template <>
    bool invert<3,float>(float invm[3*3], const float m[3*3], float *detm)
    {
        return invert3(invm, m, detm);
    }

    template <>
    bool invert<3,double>(double invm[3*3], const double m[3*3], double *detm)
    {
        return invert3(invm, m, detm);
    }

}

template <int N, int C, typename S>
bool liegroups::solve(S A[N*N], S B[N*C])
{
    for (int i=0; i<N; ++i) {
        // Find pivot
        int pivot = i;
        S pv = liegroups::abs(A[i*(N+1)]);
        for (int j=i+1; j<N; ++j) {
            S aj = liegroups::abs(A[j*N+i]);
            if (aj > pv) {
                pv = aj;
                pivot = j;
            }
        }
        if (pv == (S)0)
            return false;
        
        // Swap rows
        if (pivot != i) {
            swap(&A[i*(N+1)], &A[pivot*N + i], N-i);
            swap(&B[i*C], &B[pivot*C], C);
        }

        const S inv_pivot = (S)1 / A[i*(N+1)];

        // B
        {
            // Scale this row of B
            S *pb = &B[i*C];        
            for (int j=0; j<C; ++j)
                pb[j] *= inv_pivot;

            // Subtract from other rows of B
            for (int r=0; r<N; ++r) {
                if (r == i)
                    continue;
                S* pr = &B[r*C];
                const S f = A[r*N + i];
                for (int j=0; j<C; ++j)
                    pr[j] -= f * pb[j];                
            }
        }

        // A
        {
            // Scale this row of A
            S *pa = &A[i*N];
            //pa[i] = (S)1;
            for (int j=i+1; j<N; ++j)
                pa[j] *= inv_pivot;
        
            // Subtract from other rows of A
            for (int r=0; r<N; ++r) {
                if (r == i)
                    continue;
                S *pr = &A[r*N];
                const S f = pr[i];
                //pr[i] = 0;
                for (int j=i+1; j<N; ++j)
                    pr[j] -= f * pa[j];
            }
        }
    }
    return true;
}

template <int N, typename S>
bool liegroups::expm(S em[N*N], const S m[N*N])
{
    S scale = max_norm<N*N>(m);
    S factor = (S)1;
    int s = 0;
    if (scale > (S)0.5) {
        const S inv_ln2 = (S)1.44269504089;
        S lg_scale = liegroups::ln(scale) * inv_ln2;
        s = (int)lg_scale + 1;
        factor = (S)1 / (S)(1 << s);
    }

    S sm[N*N];
    for (int i=0; i<N*N; ++i)
        sm[i] = m[i] * factor;
    
    S *const sm2 = em;
    mat_mult<N,N,N>(sm2, sm, sm);

    static const S c[7] = {S(1.0/2), S(3.0/26), S(5.0/312), S(5.0/3432),
                           S(1.0/11440), S(1.0/308880), S(1.0/17297280) };
    
    S tmp[N*N], A[N*N], B[N*N];
    // Construct A
    {
        for (int i=0; i<N*N; ++i)
            A[i] = sm2[i]*c[5];
        for (int i=0; i<N; ++i)
            A[i*(N+1)] += c[3];
        mat_mult<N,N,N>(tmp, sm2, A);
        for (int i=0; i<N; ++i)
            tmp[i*(N+1)] += c[1];
        mat_mult<N,N,N>(A, sm2, tmp);
        for (int i=0; i<N; ++i)
            A[i*(N+1)] += (S)1;
    }
    // Construct B
    {
        for (int i=0; i<N*N; ++i)
            tmp[i] = sm2[i]*c[6];
        for (int i=0; i<N; ++i)
            tmp[i*(N+1)] += c[4];
        mat_mult<N,N,N>(B, sm2, tmp);
        for (int i=0; i<N; ++i)
            B[i*(N+1)] += c[2];
        mat_mult<N,N,N>(tmp, sm2, B);
        for (int i=0; i<N; ++i)
            tmp[i*(N+1)] += c[0];
        mat_mult<N,N,N>(B, sm, tmp);
    }

    // Compute tmp = inv(A-B) * (A+B)
    for (int i=0; i<N*N; ++i) {
        tmp[i] = A[i] + B[i];
        A[i] -= B[i];
    }
    if (!solve<N,N>(A, tmp))
        return false;

    // Square s times
    S *in = tmp, *out = A;
    for (int i=0; i<s; ++i) {        
        mat_mult<N,N,N>(out, in, in);
        S *p = in;
        in = out;
        out = p;
    }

    copy<N*N>(em, in);
    return true;
}


template <int N, typename S>
bool liegroups::sqrtm(S s[N*N], const S m[N*N])
{
    S y[N*N], z[N*N], x[N*N];

    for (int i=0; i<N*N; ++i)
        y[i] = (S)0.5*m[i];

    if (!invert<N>(z, m))
        return false;

    for (int i=0; i<N*N; ++i)
        z[i] *= (S)0.5;
        
    for (int i=0; i<N; ++i) {
        y[i*(N+1)] += (S)0.5;
        z[i*(N+1)] += (S)0.5;
    }

    for (int pass=0; pass<20; ++pass) {
        if (!invert<N>(s, z) || !invert<N>(x, y))
            return false;

        for (int i=0; i<N*N; ++i) {
            y[i] = (S)0.5*(y[i] + s[i]);
            z[i] = (S)0.5*(z[i] + x[i]);
        }

        // Check for convergence
        S yz[N*N];
        mat_mult<N,N,N>(yz, y, z);
        for (int i=0; i<N; ++i)
            yz[i*(N+1)] -= (S)1;
        
        if (max_norm<N*N>(yz) < Constants<S>::epsilon() * (S)10) {
            copy<N*N>(s, y);
            return true;
        }
    }
    return false;
}

template <int N, typename S>
bool liegroups::logm(S lm[N*N], const S m[N*N])
{
    S tmp1[N*N], tmp2[N*N];
    S *x = tmp1, *y = tmp2;
    copy<N*N>(x, m);

    int s = 0;
    S factor = (S)1;
    const int MAX_PASSES = 20;
    int pass;
    for (pass=0; pass<MAX_PASSES; ++pass) {        
        S mv = 0;
        for (int i=0; i<N; ++i) {
            for (int j=0; j<N; ++j) {
                S xi = x[i*N+j];
                if (i == j)
                    xi -= (S)1;
                S axi = liegroups::abs(xi);
                mv = liegroups::max(axi, mv);
            }
        }
        std::cout << "mv = " << mv << std::endl;
        if (mv < (S)0.5)
            break;
        
        if (!sqrtm<N>(y, x))
            return false;
        swap(x, y);

        ++s;
        factor *= (S)2;
    }
    if (pass == MAX_PASSES)
        return false;

    for (int i=0; i<N; ++i)
        x[i*(N+1)] -= (S)1;

    static const S a[7] = {(S)1, (S)3, S(535.0/156),
                           S(145.0/78), S(1377.0/2860), S(223.0/4290),
                           S(11.0/7280)};
    
    static const S b[7] = {S(7.0/2), S(63.0/13), S(175.0/52),
                           S(175.0/143), S(63.0/286), S(7.0/429),
                           S(1.0/3432)};
        
    S A[N*N], B[N*N];
    {
        // A0 = x
        copy<N*N>(A, x);
        // B0 = b[0]*x + I
        for (int i=0; i<N*N; ++i)
            B[i] = x[i]*b[0];
        for (int i=0; i<N; ++i)
            B[i*(N+1)] += (S)1;
    }

    // y = x^(term+1)
    mat_mult<N,N,N>(y, x, x);
    S *z = lm;
    
    for (int term=1; term<7; ++term) {
        S fa = a[term];
        S fb = b[term];
        for (int i=0; i<N*N; ++i) {
            A[i] += y[i]*fa;
            B[i] += y[i]*fb;
        }
        if (term<6) {
            mat_mult<N,N,N>(z, y, x);
            swap(y, z);
        }
    }

    if (!solve<N,N>(B,A))
        return false;

    for (int i=0; i<N*N; ++i)
        lm[i] = factor * A[i];
    
    return true;
}
