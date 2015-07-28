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
static S max_abs(const S x[N])
{
    S mv = 0;
    for (int i=0; i<N; ++i) {
        S axi = liegroups::abs(x[i]);
        mv = liegroups::max(axi, mv);
    }
    return mv;
}

template <int N, typename S>
static S max_row_norm(const S x[N*N])
{
    S mv = 0;
    for (int i=0; i<N; ++i) {
        S sum = 0;
        for (int j=0; j<N; ++j)
            sum += liegroups::abs(x[i*N+j]);        
        mv = liegroups::max(sum, mv);
    }
    return mv;
}

namespace liegroups
{
    template <int N>
    struct MatMultSquare
    {
        template <typename S>
        static void eval(S ab[N*N], const S a[N*N], const S b[N*N])
        {
            for (int i=0; i<N; ++i) {
                for (int j=0; j<N; ++j) {
                    S sum = (S)0;
                    for (int k=0; k<N; ++k) {
                        sum += a[i*N + k] * b[j + k*N];
                    }
                    ab[i*N + j] = sum;
                }
            }
        }        
    };

    template <>
    struct MatMultSquare<3>
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
        
template <int N, typename S>
void liegroups::mat_mult_square(S ab[N*N], const S a[N*N], const S b[N*N])
{
    MatMultSquare<N>::eval(ab, a, b);
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
    const S m0 = m[0];
    invm[0] =  inv_det * m[3];
    invm[1] = -inv_det * m[1];
    invm[2] = -inv_det * m[2];
    invm[3] =  inv_det * m0;
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

template <int N, typename S>
bool liegroups::LU_decompose(S A[N*N], int index[N])
{
    for (int i=0; i<N; ++i)
        index[i] = i;
    
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

        S *pi = &A[i*N];
        
        // Swap rows
        if (pivot != i) {
            swap(index[i], index[pivot]);
            swap(pi, &A[pivot*N], N);
        }

        const S inv_pivot = (S)1 / pi[i];

        // Scale this row of A
        pi[i] = inv_pivot;
        for (int j=i+1; j<N; ++j)
            pi[j] *= inv_pivot;
        
        // Subtract from other rows of A
        for (int r=i+1; r<N; ++r) { 
            S *pr = &A[r*N];
            const S f = pr[i];
            for (int j=i+1; j<N; ++j)
                pr[j] -= f * pi[j];
        }        
    }
    return true;
}

// b and x have stride C
template <int N, int C, typename S>
void liegroups::LU_inverse_times_vec(S x[N*C], const S LU[N*N], const int index[N], const S b[N*C])
{
    S y[N];
    for (int i=0; i<N; ++i) {
        const S *Li = &LU[i*N];
        S yi = b[index[i]*C];
        for (int j=0; j<i; ++j)
            yi -= Li[j] * y[j];
        y[i] = Li[i] * yi;        
    }
    for (int i=N-1; i>=0; --i) {
        const S *Ui = &LU[i*N];
        S xi = y[i];
        for (int j=i+1; j<N; ++j)
            xi -= Ui[j]*x[j*C];
        x[i*C] = xi;
    }
}

template <int N, typename S>
void liegroups::LU_invert(S inv[N*N], const S LU[N*N], const int index[N])
{
    S b[N];
    for (int i=0; i<N; ++i) {
        b[i] = 0;
    }
    
    for (int c=0; c<N; ++c) {
        b[c] = 1;
        LU_inverse_times_vec<N,1>(&inv[c*N], LU, index, b);
        b[c] = 0;
    }

    // Transpose
    for (int i=0; i<N; ++i) {
        for (int j=i+1; j<N; ++j) {
            S tmp = inv[i*N + j];
            inv[i*N + j] = inv[j*N + i];
            inv[j*N + i] = tmp;
        }
    }
}

template <typename S>
static S lg(S x)
{
    const S inv_ln2 = (S)1.44269504088896;
    return liegroups::ln(x) * inv_ln2;
}

template <int N, typename S>
static bool balance(const S m[N*N],
                    S d[N], S inv_d[N])
{
    S bm[N*N];
    for (int i=0; i<N*N; ++i)
        bm[i] = liegroups::abs(m[i]);

    for (int i=0; i<N; ++i)
        inv_d[i] = d[i] = (S)1;

    bool changed = false;
    
    for (int pass=0; pass < N*2; ++pass) {
        S sc[N], sr[N];
        for (int i=0; i<N; ++i)
            sc[i] = sr[i] = (S)0;
        
        for (int i=0; i<N; ++i) {
            for (int j=0; j<N; ++j) {
                sc[i] += bm[j*N+i];
                sr[i] += bm[i*N+j];                
            }
        }

        S max_a = (S)2, max_b = (S)1;
        int argmax = -1;
        for (int i=0; i<N; ++i) {
            S a = sr[i], b = sc[i];
            if (a == (S)0 || b == (S)0)
                continue;
            if (a < b) {
                if (max_a < max_b) {
                    if (a * max_b > max_a * b)
                        continue;
                } else {
                    if (a * max_a > b * max_b)
                        continue;
                }
            } else {
                if (max_a < max_b) {
                    if (a * max_a < max_b * b)
                        continue;
                } else if (a * max_b < b * max_a)
                    continue;
            }
            max_a = a;
            max_b = b;
            argmax = i;
        }

        if (argmax == -1)
            break;

        S r = liegroups::sqrt(max_a / max_b);
        S inv_r = (S) 1 / r;
        
        for (int i=0; i<N; ++i) {
            bm[argmax*N + i] *= inv_r;
            bm[i*N + argmax] *= r;
        }

        d[argmax] *= r;
        inv_d[argmax] *= inv_r;
        changed = true;
    }

    if (changed) {
        S det_d = d[0];
        for (int i=1; i<N; ++i)
            det_d *= d[i];
        S scale = liegroups::pow(det_d, (S)1/(S)N);
        S inv_scale = 1 / scale; 
        for (int i=1; i<N; ++i) {
            d[i] *= inv_scale;
            inv_d[i] *= scale;
        }
    }
    // std::cerr << "d = [ ";
    // for (int i=0; i<N; ++i) {
    //     std::cerr.width(12);
    //     std::cerr << d[i] << " ";
    // }
    // std::cerr << std::endl;

    return changed;
}


template <int N, typename S>
bool liegroups::expm(S em[N*N], const S m[N*N])
{
    S sm[N*N];
    copy<N*N>(sm, m);
    
    S trace = m[0];
    for (int i=1; i<N; ++i)
        trace += m[i*(N+1)];
    S mu = trace * ((S)1 / (S)N);
    for (int i=0; i<N; ++i)
        sm[i*(N+1)] -= mu;
    
    S scale = max_row_norm<N>(sm);
    S pre_scale = scale;
    bool balanced = false;
    S d[N], inv_d[N];
    int s = 0;
    if (scale > (S)1.9) {
        balanced = balance<N>(sm, d, inv_d);        
        if (balanced) {
            for (int i=0,k=0; i<N; ++i) {
                for (int j=0; j<N; ++j, ++k) {
                    sm[k] *= inv_d[i] * d[j];
                }
            }
            scale = max_row_norm<N>(sm);
            if (scale > pre_scale) {
                copy<N*N>(sm, m);
                for (int i=0; i<N; ++i)
                    sm[i*(N+1)] -= mu;
                scale = pre_scale;
                balanced = false;
                //std::cerr << "balancing failed" << std::endl;
            }
        }
        
        if (scale > (S)1.9) {
            S lg_scale = lg(scale);
            s = (int)lg_scale;
            if (s > 0) {
                S factor = liegroups::pow((S)2, (S)-s);                
                for (int i=0; i<N*N; ++i)
                    sm[i] *= factor;
            }
        }        
        //std::cerr << "s = " << s << std::endl;
    } else {
        //std::cerr << "no scale " << std::endl;
    }
    
    S *const sm2 = em;
    mat_mult_square<N>(sm2, sm, sm);

    static const S c[7] = {S(1.0/2), S(3.0/26), S(5.0/312), S(5.0/3432),
                           S(1.0/11440), S(1.0/308880), S(1.0/17297280) };
    
    S tmp[N*N], A[N*N], B[N*N];
    // Construct A
    {
        for (int i=0; i<N*N; ++i)
            A[i] = sm2[i]*c[5];
        for (int i=0; i<N; ++i)
            A[i*(N+1)] += c[3];
        mat_mult_square<N>(tmp, sm2, A);
        for (int i=0; i<N; ++i)
            tmp[i*(N+1)] += c[1];
        mat_mult_square<N>(A, sm2, tmp);
        for (int i=0; i<N; ++i)
            A[i*(N+1)] += (S)1;
    }
    // Construct B
    {
        for (int i=0; i<N*N; ++i)
            tmp[i] = sm2[i]*c[6];
        for (int i=0; i<N; ++i)
            tmp[i*(N+1)] += c[4];
        mat_mult_square<N>(B, sm2, tmp);
        for (int i=0; i<N; ++i)
            B[i*(N+1)] += c[2];
        mat_mult_square<N>(tmp, sm2, B);
        for (int i=0; i<N; ++i)
            tmp[i*(N+1)] += c[0];
        mat_mult_square<N>(B, sm, tmp);
    }

    // Compute tmp = inv(A-B) * (A+B)
    for (int i=0; i<N*N; ++i) {
        tmp[i] = A[i] + B[i];
        A[i] -= B[i];
    }
    int index[N];
    if (!LU_decompose<N>(A, index))
        return false;

    LU_inverse_times_mat<N,N>(tmp, A, index, tmp);
    
    // Square s times
    S *in = tmp, *out = A;
    for (int i=0; i<s; ++i) {        
        mat_mult_square<N>(out, in, in);
        swap(in, out);
    }

    S exp_mu = liegroups::exp(mu);
    
    if (balanced) {
        // Undo balancing
        for (int i=0,k=0; i<N; ++i) {
            for (int j=0; j<N; ++j, ++k) {
                em[k] = exp_mu * (in[k] * inv_d[j] * d[i]);
            }
        }
        //std::cerr << pre_scale << " --> " << scale << std::endl;        
    } else {
        for (int i=0; i<N*N; ++i)
            em[i] = exp_mu * in[i];
    }
    
    return true;
}

template <int N, typename S>
static bool LU_invert(S invm[N*N], const S m[N*N], S* det = 0)
{
    S lu[N*N];
    copy<N*N>(lu, m);
    int index[N];
    if (!liegroups::LU_decompose<N>(lu, index))
        return false;

    if (det) {
        S d = lu[0];
        for (int i=1; i<N; ++i)
            d *= lu[i*(N+1)];
        
        *det = (S)1 / d;
    }

    for (int i=0,k=0; i<N-1; ++i) {
        invm[k++] = (S)1;
        for (int j=0; j<N; ++j)
            invm[k++] = (S)0;        
    }
    invm[N*N-1] = (S)1;

    liegroups::LU_inverse_times_mat<N,N>(invm, lu, index, invm);
    return true;
}

template <int N, typename S>
bool liegroups::sqrtm(S s[N*N], const S m[N*N], const S tol)
{
    S y[N*N], z[N*N], x[N*N];

    const S lambda_p = S(-1) / (S)(2*N);
    const S mag = max_abs<N*N>(m);
    //std::cerr << "mag = " << mag << std::endl;
    // y1, z1
    {
        S detY;
        if (!invert<N>(z, m, &detY))
            return false;

        const S lambda = liegroups::pow(liegroups::abs(detY), lambda_p);
        const S inv_lambda = (S)1 / lambda;
        const S a = (S)0.5*lambda;
        const S b = (S)0.5*inv_lambda;
        
        for (int i=0; i<N*N; ++i)
            y[i] = a*m[i];

        for (int i=0; i<N*N; ++i)
            z[i] *= b;
        
        for (int i=0; i<N; ++i) {
            y[i*(N+1)] += b;
            z[i*(N+1)] += a;
        }
    }
    
    S last_err = (S)-1;
    
    for (int pass=0; pass<20; ++pass) {
        S detY, detZ;
        if (!invert<N>(s, z, &detZ) || !invert<N>(x, y, &detY))
            break;

        const S lambda = liegroups::pow(liegroups::abs(detY * detZ), lambda_p);
        const S inv_lambda = (S)1 / lambda;
        const S a = (S)0.5*lambda;
        const S b = (S)0.5*inv_lambda;
        
        for (int i=0; i<N*N; ++i) {
            y[i] = a*y[i] + b*s[i];
            z[i] = a*z[i] + b*x[i];
        }

        if (pass >= 4) {            
            // Check for convergence
            S yy[N*N];
            mat_mult_square<N>(yy, y, y);
            for (int i=0; i<N*N; ++i)
                yy[i] -= m[i];
        
            S err = max_abs<N*N>(yy);
            if (last_err > 0 && err >= last_err) {
                //std::cerr << "non-decreasing error: " << last_err << " --> " << err << std::endl;
                break;
            }
            //std::cerr << pass << "\t" << err << std::endl;
            if (err < tol * mag) {
                copy<N*N>(s, y);
                return true;
            }
            last_err = err;
        }
    }

    if (0)
    {
        S yy[N*N];
        mat_mult_square<N>(yy, y, y);
        for (int i=0; i<N*N; ++i)
            yy[i] -= m[i];
        std::cerr.precision(19);
        std::cerr << "mag = " << mag << std::endl;
        std::cerr << "max err = " << tol * mag << std::endl;
        for (int i=0; i<N; ++i) {
            for (int j=0; j<N; ++j) {
                std::cerr.width(25);
                std::cerr << yy[i*N+j];
            }
            std::cerr << std::endl;
        }
    }
    //std::cerr << "sqrtm failed" << std::endl;
    
    return false;
}

template <int N, typename S>
static bool logm_newton(S lm[N*N], const S m[N*N], int iter)
{    
    S nlm[N*N];
    for (int i=0; i<N*N; ++i)
        nlm[i] = -lm[i];
    
    for (int pass=0; pass<iter; ++pass) {        
        S enlm[N*N];
        if (!liegroups::expm<N>(enlm, nlm))
            return false;

        S d[N*N];
        liegroups::mat_mult_square<N>(d, enlm, m);
        
        for (int i=0; i<N; ++i)
            d[i*(N+1)] -= (S)1;
        
        for (int i=0; i<N*N; ++i)
            nlm[i] -= d[i];
    }
    
    for (int i=0; i<N*N; ++i)
        lm[i] = -nlm[i];
    
    return true;
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
        S x_minus_id[N*N];
        copy<N*N>(x_minus_id, x);
        for (int i=0; i<N; ++i)
            x_minus_id[i*(N+1)] -= (S)1;
        
        if (max_row_norm<N>(x_minus_id) < (S)0.25)
            break;
        
        if (!sqrtm<N>(y, x, Constants<S>::epsilon()*(S)(10*N)))
            return false;
        swap(x, y);

        ++s;
        factor *= (S)2;
    }
    if (pass == MAX_PASSES)
        return false;

    for (int i=0; i<N; ++i)
        x[i*(N+1)] -= (S)1;

    // Partial fraction coefficients for Pade of log(1+x), degree 7
    static const S alpha[7] = {
        (S)0.2089795918367347,
        (S)0.1909150252525595,
        (S)0.1909150252525595,
        (S)0.1398526957446383,
        (S)0.1398526957446383,
        (S)0.0647424830844349,
        (S)0.0647424830844349
    };
    
    static const S beta[7] = {
        (S)0.5000000000000000,
        (S)0.7029225756886985,
        (S)0.2970774243113014,
        (S)0.1292344072003028,
        (S)0.8707655927996972,
        (S)0.0254460438286208,
        (S)0.9745539561713792
    };

    S A[N*N], B[N*N];
    int index[N];

    for (int i=0; i<N*N; ++i)
        lm[i] = (S)0;
    
    for (int term=0; term<7; ++term)
    {        
        for (int i=0; i<N*N; ++i)
            B[i] = beta[term] * x[i];
        for (int i=0; i<N; ++i)
            B[i*(N+1)] += (S)1;
        if (!LU_decompose<N>(B, index))
            return false;
        LU_inverse_times_mat<N,N>(A, B, index, x);
        const S f = factor * alpha[term];
        for (int i=0; i<N*N; ++i)
            lm[i] += f * A[i];
    }
    return logm_newton<N>(lm, m, 2);
}
