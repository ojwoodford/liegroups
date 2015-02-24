#include <liegroups/so3.hpp>
#include <liegroups/exp_coefs.hpp>
#include <liegroups/exp_helpers.hpp>
#include <liegroups/scalar.hpp>
#include <liegroups/matrix.hpp>
#include <cmath>

template <> const liegroups::SO3<float>
liegroups::SO3<float>::identity = { {1.f, 0.f, 0.f,
                                     0.f, 1.f, 0.f,
                                     0.f, 0.f, 1.f} };

template <> const liegroups::SO3<double>
liegroups::SO3<double>::identity = { {1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0} };

template <typename S>
void liegroups::SO3<S>::ad_multiply(S ada_b[], const S a[], const S b[])
{
    const S b0 = b[0];
    const S b1 = b[1];
    const S b2 = b[2];
    ada_b[0] = a[1]*b2 - a[2]*b1;
    ada_b[1] = a[2]*b0 - a[0]*b2;
    ada_b[2] = a[0]*b1 - a[1]*b0;
}

template void liegroups::SO3<float>::ad_multiply(float[], const float[], const float[]);
template void liegroups::SO3<double>::ad_multiply(double[], const double[], const double[]);

template <typename S>
void liegroups::SO3<S>::ad(S ada[3*3], const S a[])
{
    ada[0] = ada[4] = ada[8] = 0;
    ada[1] = -a[2];
    ada[2] = a[1];
    ada[3] = a[2];
    ada[5] = -a[0];
    ada[6] = -a[1];
    ada[7] = a[0];
}

template void liegroups::SO3<float>::ad(float[], const float[]);
template void liegroups::SO3<double>::ad(double[], const double[]);

template <class S>
void liegroups::multiply(SO3<S> &ab, const SO3<S> &a, const SO3<S> &b)
{
    mat_mult_square<3>(ab.R, a.R, b.R);
}

template void liegroups::multiply<float>(SO3<float>&, const SO3<float>&, const SO3<float> &);
template void liegroups::multiply<double>(SO3<double>&, const SO3<double>&, const SO3<double> &);


template <typename S1, typename S2>
static S2 dot3(const S1 a[3], const S2 b[3])
{
    return (S2)(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

template <class S>
void liegroups::multiply_a_binv(SO3<S> &abinv, const SO3<S> &a, const SO3<S> &b)
{
    abinv.R[0] = dot3(&a.R[0], &b.R[0]);
    abinv.R[1] = dot3(&a.R[0], &b.R[3]);
    abinv.R[2] = dot3(&a.R[0], &b.R[6]);
    abinv.R[3] = dot3(&a.R[3], &b.R[0]);
    abinv.R[4] = dot3(&a.R[3], &b.R[3]);
    abinv.R[5] = dot3(&a.R[3], &b.R[6]);
    abinv.R[6] = dot3(&a.R[6], &b.R[0]);
    abinv.R[7] = dot3(&a.R[6], &b.R[3]);
    abinv.R[8] = dot3(&a.R[6], &b.R[6]);
}

template void liegroups::multiply_a_binv<float>(SO3<float>&, const SO3<float>&, const SO3<float> &);
template void liegroups::multiply_a_binv<double>(SO3<double>&, const SO3<double>&, const SO3<double> &);

template <class S>
liegroups::SO3<S> liegroups::inverse(const SO3<S> &g)
{
    const S r1 = g.R[1], r2 = g.R[2], r5 = g.R[5];
    
    SO3<S> ginv;
    ginv.R[0] = g.R[0];
    ginv.R[1] = g.R[3];
    ginv.R[2] = g.R[6];
    ginv.R[3] = r1;
    ginv.R[4] = g.R[4];
    ginv.R[5] = g.R[7];
    ginv.R[6] = r2;
    ginv.R[7] = r5;
    ginv.R[8] = g.R[8];
    return ginv;
}

template liegroups::SO3<float> liegroups::inverse<float>(const SO3<float>&);
template liegroups::SO3<double> liegroups::inverse<double>(const SO3<double>&);

template <class S>
void liegroups::invert(SO3<S> &g)
{
    g = inverse(g);
}

template void liegroups::invert<float>(SO3<float>&);
template void liegroups::invert<double>(SO3<double>&);

template <typename S>
static void normalize3(S x[3])
{
    S xx = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    S factor = (S)1 / (S)liegroups::sqrt(xx);
    x[0] *= factor;
    x[1] *= factor;
    x[2] *= factor;
}

template <class S>
void liegroups::rectify(SO3<S> &g)
{
    normalize3(&g.R[0]);
    
    S xy = dot3(&g.R[0], &g.R[3]);
    g.R[3] -= xy*g.R[0];
    g.R[4] -= xy*g.R[1];
    g.R[5] -= xy*g.R[2];
    normalize3(&g.R[3]);

    g.R[6] = g.R[1]*g.R[5] - g.R[2]*g.R[4];
    g.R[7] = g.R[2]*g.R[3] - g.R[0]*g.R[5];
    g.R[8] = g.R[0]*g.R[4] - g.R[1]*g.R[3];
}

template void liegroups::rectify<float>(SO3<float>&);
template void liegroups::rectify<double>(SO3<double>&);

template <class S, class X>
void liegroups::transform_point(X y[3], const SO3<S> &g, const X x[3])
{
    S y0 = dot3(&g.R[0], x);
    S y1 = dot3(&g.R[3], x);
    S y2 = dot3(&g.R[6], x);
    y[0] = y0;
    y[1] = y1;
    y[2] = y2;
}

template void liegroups::transform_point<float,float>(float[3], const SO3<float> &, const float[3]);
template void liegroups::transform_point<double,double>(double[3], const SO3<double> &, const double[3]);
template void liegroups::transform_point<double,float>(float[3], const SO3<double> &, const float[3]);
template void liegroups::transform_point<float,double>(double[3], const SO3<float> &, const double[3]);

template <class S, class X>
void liegroups::transform_point_by_inverse(X y[3], const SO3<S> &g, const X x[3])
{
    S y0 = g.R[0]*x[0] + g.R[3]*x[1] + g.R[6]*x[2];
    S y1 = g.R[1]*x[0] + g.R[4]*x[1] + g.R[7]*x[2];
    S y2 = g.R[2]*x[0] + g.R[5]*x[1] + g.R[8]*x[2];
    y[0] = y0;
    y[1] = y1;
    y[2] = y2;
}

template void liegroups::transform_point_by_inverse<float,float>(float[3], const SO3<float> &, const float[3]);
template void liegroups::transform_point_by_inverse<double,double>(double[3], const SO3<double> &, const double[3]);
template void liegroups::transform_point_by_inverse<double,float>(float[3], const SO3<double> &, const float[3]);
template void liegroups::transform_point_by_inverse<float,double>(double[3], const SO3<float> &, const double[3]);

template <class S>
void liegroups::exp(SO3<S> &X, const S w[3])
{
    const S theta_sq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
    const ExpCoefs<S> coefs(theta_sq);
    compute_exp_matrix3(X.R, coefs.cos_theta, coefs.A, coefs.B, w);
}

template void liegroups::exp<float>(SO3<float> &, const float[3]);
template void liegroups::exp<double>(SO3<double> &, const double[3]);

template <class S>
void liegroups::log(S w[3], const SO3<S> &X)
{
    w[0] = (S)0.5 * (X.R[7] - X.R[5]);
    w[1] = (S)0.5 * (X.R[2] - X.R[6]);
    w[2] = (S)0.5 * (X.R[3] - X.R[1]);
    
    S tr = X.R[0] + X.R[4] + X.R[8];
    S ct = (S)0.5  * (tr - (S)1);
    S st2 = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];

    if (ct > (S)0.999856) {
        // Small angles
        // Taylor expansion of f(x) = arcsin(x) / x
        // x^2 = st2
        S f = (S)1 + st2*((S)(1/6.0) + st2*((S)(3/40.0) + st2*(S)(5/112.0)));
        w[0] *= f;
        w[1] *= f;
        w[2] *= f;
        return;
    }

    if (ct > (S)-0.99) {
        S theta = liegroups::acos(ct);
        S st = liegroups::sqrt(st2);
        S factor = theta / st;
        w[0] *= factor;
        w[1] *= factor;
        w[2] *= factor;
        return;
    }
   
    // Angles near pi    
    S st = liegroups::sqrt(st2);    
    S theta = (S)3.1415926535897932 - liegroups::asin(st);
    S invB = (theta*theta) / ((S)1 - ct);

    S w00 = invB*(X.R[0] - ct);
    S w11 = invB*(X.R[4] - ct);
    S w22 = invB*(X.R[8] - ct);

    S w01 = invB*(S)0.5*(X.R[1] + X.R[3]);
    S w02 = invB*(S)0.5*(X.R[2] + X.R[6]);
    S w12 = invB*(S)0.5*(X.R[5] + X.R[7]);

    // Take sqrt of biggest element of w
    if (w00 > w11) {
        if (w00 > w22) {
            w[0] = (S)(w[0] < 0 ? -1 : 1) * liegroups::sqrt(w00);
            S inv_w0 = (S)1/w[0];
            w[1] = w01 * inv_w0;
            w[2] = w02 * inv_w0;
        } else {
            w[2] = (S)(w[2] < 0 ? -1 : 1) * liegroups::sqrt(w22);
            S inv_w2 = (S)1/w[2];
            w[0] = w02 * inv_w2;
            w[1] = w12 * inv_w2;
        }
    } else if (w11 > w22) {
        w[1] = (S)(w[1] < 0 ? -1 : 1) * liegroups::sqrt(w11);
        S inv_w1 = (S)1/w[1];
        w[0] = w01 * inv_w1;
        w[2] = w12 * inv_w1;
    } else {
        w[2] = (S)(w[2] < 0 ? -1 : 1) * liegroups::sqrt(w22);
        S inv_w2 = (S)1/w[2];
        w[0] = w02 * inv_w2;
        w[1] = w12 * inv_w2;
    }
}

template void liegroups::log<float>(float[3], const SO3<float> &);
template void liegroups::log<double>(double[3], const SO3<double> &);

template <class S>
void liegroups::adjoint(S adj[3*3], const SO3<S> &g)
{
    for (int i=0; i<9; ++i)
        adj[i] = g.R[i];
}

template void liegroups::adjoint<float>(float[3*3], const SO3<float> &);
template void liegroups::adjoint<double>(double[3*3], const SO3<double> &);

template <class S>
void liegroups::adjoint_multiply(S y[3], const SO3<S> &g, const S x[3])
{
    transform_point(y, g, x);
}

template void liegroups::adjoint_multiply<float>(float[3], const SO3<float> &, const float[3]);
template void liegroups::adjoint_multiply<double>(double[3], const SO3<double> &, const double[3]);

template <class S>
void liegroups::adjoint_T_multiply(S y[3], const SO3<S> &g, const S x[3])
{
    transform_point_by_inverse(y, g, x);
}

template void liegroups::adjoint_T_multiply<float>(float[3], const SO3<float> &, const float[3]);
template void liegroups::adjoint_T_multiply<double>(double[3], const SO3<double> &, const double[3]);

template <class S>
static void unit_perp3(S p[3], const S v[3])
{
    const S v00 = v[0]*v[0], v11 = v[1]*v[1], v22 = v[2]*v[2];
    const S v01 = v[0]*v[1], v02 = v[0]*v[2], v12 = v[1]*v[2];
    if (v00 < v11) {
        if (v00 < v22) {
            p[0] = (S)1 - v00;
            p[1] = -v01;
            p[2] = -v02;
        } else {
            p[0] = -v02;
            p[1] = -v12;
            p[2] = (S)1 - v22;
        }
    } else if (v00 < v22) {
        if (v00 < v11) {
            p[0] = (S)1 - v00;
            p[1] = -v01;
            p[2] = -v02;
        } else {
            p[0] = -v01;
            p[1] = (S)1 - v11;
            p[2] = -v12;
        }
    } else if (v11 < v22) {
            p[0] = -v01;
            p[1] = (S)1 - v11;
            p[2] = -v12;
    } else {
        p[0] = -v02;
        p[1] = -v12;
        p[2] = (S)1 - v22;
    }
    
    normalize3(p);
}

template <class S>
bool liegroups::compute_rotation_between_unit_vectors(SO3<S> &R, const S a[3], const S b[3])
{
    const S cos_theta = dot3(a,b);

    if (cos_theta < (S)-0.9) {
        const S neg_a[3] = {-a[0], -a[1], -a[2] };
        SO3<S> neg_a_to_b;
        if (!compute_rotation_between_unit_vectors(neg_a_to_b, neg_a, b))
            return false;

        SO3<S> a_to_neg_a;
        S p[3];
        unit_perp3(p, a);
        const S C = (S)2;
        const S Cp00 = C*p[0]*p[0], Cp11 = C*p[1]*p[1], Cp22 = C*p[2]*p[2];
        const S Cp01 = C*p[0]*p[1], Cp02 = C*p[0]*p[2], Cp12 = C*p[1]*p[2];
        a_to_neg_a.R[0] = (S)1 - (Cp11 + Cp22);
        a_to_neg_a.R[4] = (S)1 - (Cp00 + Cp22);
        a_to_neg_a.R[8] = (S)1 - (Cp00 + Cp11);
        a_to_neg_a.R[1] = Cp01;
        a_to_neg_a.R[2] = Cp02;
        a_to_neg_a.R[3] = Cp01;
        a_to_neg_a.R[5] = Cp12;
        a_to_neg_a.R[6] = Cp02;
        a_to_neg_a.R[7] = Cp12;

        multiply(R, neg_a_to_b, a_to_neg_a);
        return true;
    }

    const S C = (S)1 / ((S)1 + cos_theta);
    
    const S w[3] = {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
        
    const S w00 = w[0]*w[0];
    const S w11 = w[1]*w[1];
    const S w22 = w[2]*w[2];
    R.R[0] = (S)1 - C*(w11 + w22);
    R.R[4] = (S)1 - C*(w00 + w22);
    R.R[8] = (S)1 - C*(w00 + w11);
    const S Cw01 = C*w[0]*w[1];
    const S Cw02 = C*w[0]*w[2];
    const S Cw12 = C*w[1]*w[2];
    R.R[1] = Cw01 - w[2];
    R.R[2] = Cw02 + w[1];
    R.R[3] = Cw01 + w[2];
    R.R[5] = Cw12 - w[0];
    R.R[6] = Cw02 - w[1];
    R.R[7] = Cw12 + w[0];
    return true;
}

template bool liegroups::compute_rotation_between_unit_vectors<float>(SO3<float> &, const float[3], const float[3]);
template bool liegroups::compute_rotation_between_unit_vectors<double>(SO3<double> &, const double[3], const double[3]);


template <class S>
void liegroups::exp_diff(SO3<S> &X, S dexp[3*3], const S w[3])
{
    ExpCoefs<S> coefs(dot3(w,w));

    compute_exp_matrix3(X.R, coefs.cos_theta, coefs.A, coefs.B, w);
    compute_exp_matrix3(dexp, coefs.A, coefs.B, coefs.C, w);
}

template void liegroups::exp_diff<float>(SO3<float>&, float[3*3], const float [3]);
template void liegroups::exp_diff<double>(SO3<double>&, double[3*3], const double [3]);

template <typename S>
void liegroups::log_diff(S w[3], S dlog[3*3], const SO3<S> &X)
{
    log(w, X);
    
    const S theta_sq = dot3(w,w);
    const ExpCoefs<S> coefs(theta_sq);

    S e;
    if (coefs.cos_theta < (S)0.9) {
        e = (coefs.B - (S)0.5 * coefs.A) / (1 - coefs.cos_theta);
    } else {
        e = ((S)0.5 * coefs.B - coefs.C) / coefs.A;
    }

    compute_exp_matrix3(dlog, (S)1 - theta_sq * e, (S)-0.5, e, w);    
}

template void liegroups::log_diff<float>(float[], float[], const SO3<float>&);
template void liegroups::log_diff<double>(double[], double[], const SO3<double>&);
