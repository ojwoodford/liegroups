#include <liegroups/sim3.hpp>
#include <liegroups/scalar.hpp>
#include <liegroups/matrix.hpp>
#include <liegroups/exp_coefs.hpp>
#include <liegroups/exp_helpers.hpp>
#include <cmath>
#include <iostream>

template <> const liegroups::Sim3<float>
liegroups::Sim3<float>::identity = { liegroups::SE3<float>::identity, 1.0f, 1.0f };

template <> const liegroups::Sim3<double>
liegroups::Sim3<double>::identity = { liegroups::SE3<double>::identity, 1.0, 1.0 };

template <class S>
void liegroups::multiply(Sim3<S> &ab, const Sim3<S> &a, const Sim3<S> &b)
{
    multiply(ab.rigid, a.rigid, b.rigid);
    S tf = b.inv_scale - (S)1;
    ab.rigid.t[0] += tf*a.rigid.t[0];
    ab.rigid.t[1] += tf*a.rigid.t[1];
    ab.rigid.t[2] += tf*a.rigid.t[2];
    ab.scale = a.scale * b.scale;
    ab.inv_scale = a.inv_scale * b.inv_scale;
}

template void liegroups::multiply<float>(Sim3<float>&, const Sim3<float>&, const Sim3<float> &);
template void liegroups::multiply<double>(Sim3<double>&, const Sim3<double>&, const Sim3<double> &);

template <class S>
void liegroups::multiply_a_binv(Sim3<S> &abinv, const Sim3<S> &a, const Sim3<S> &b)
{
    multiply_a_binv(abinv.rigid, a.rigid, b.rigid);
    abinv.rigid.t[0] *= b.scale;
    abinv.rigid.t[1] *= b.scale;
    abinv.rigid.t[2] *= b.scale;
    abinv.scale = a.scale * b.inv_scale;
    abinv.inv_scale = a.inv_scale * b.scale;
}

template void liegroups::multiply_a_binv<float>(Sim3<float>&, const Sim3<float>&, const Sim3<float> &);
template void liegroups::multiply_a_binv<double>(Sim3<double>&, const Sim3<double>&, const Sim3<double> &);

template <class S>
void liegroups::invert(Sim3<S> &g)
{
    invert(g.rigid);
    S s = g.scale;
    g.rigid.t[0] *= s;
    g.rigid.t[1] *= s;
    g.rigid.t[2] *= s;
    g.scale = g.inv_scale;
    g.inv_scale = s;
}

template void liegroups::invert<float>(Sim3<float>&);
template void liegroups::invert<double>(Sim3<double>&);

template <class S>
void liegroups::rectify(Sim3<S> &g)
{
    rectify(g.rigid);
    g.inv_scale = (S)1 / g.scale;
}

template void liegroups::rectify<float>(Sim3<float>&);
template void liegroups::rectify<double>(Sim3<double>&);

template <class S, class X>
void liegroups::transform_point(X y[3], const Sim3<S> &g, const X x[3])
{
    transform_point(y, g.rigid, x);
    y[0] *= g.scale;
    y[1] *= g.scale;
    y[2] *= g.scale;
}

template void liegroups::transform_point<float,float>(float[3], const Sim3<float> &, const float[3]);
template void liegroups::transform_point<double,double>(double[3], const Sim3<double> &, const double[3]);
template void liegroups::transform_point<double,float>(float[3], const Sim3<double> &, const float[3]);
template void liegroups::transform_point<float,double>(double[3], const Sim3<float> &, const double[3]);

template <class S, class X>
void liegroups::transform_point_by_inverse(X y[3], const Sim3<S> &g, const X x[3])
{
    y[0] = g.inv_scale * x[0];
    y[1] = g.inv_scale * x[1];
    y[2] = g.inv_scale * x[2];
    transform_point_by_inverse(y, g.rigid, y);
}

template void liegroups::transform_point_by_inverse<float,float>(float[3], const Sim3<float> &, const float[3]);
template void liegroups::transform_point_by_inverse<double,double>(double[3], const Sim3<double> &, const double[3]);
template void liegroups::transform_point_by_inverse<double,float>(float[3], const Sim3<double> &, const float[3]);
template void liegroups::transform_point_by_inverse<float,double>(double[3], const Sim3<float> &, const double[3]);

template <class S>
static void compute_V(S &a, S &b, S theta, S lambda, S ct, S st, S exp_lambda)
{
    S theta_sq = theta * theta;
    S lambda_sq = lambda * lambda;
    S z = theta_sq + lambda_sq;

    if (z < liegroups::Constants<S>::sqrt_epsilon()) {
        a = (S)1 - (S)0.5*lambda + (lambda_sq - theta_sq)*(S)(1.0/6);
        b = theta*((S)0.5 - lambda*(S)(1.0/6));
    } else {
        S inv_s = (S)1 / exp_lambda;
        S inv_z = (S)1 / z;
        a = inv_z*(lambda*ct + theta*st - lambda*inv_s);
        b = inv_z*(lambda*st - theta*ct + theta*inv_s);
    }
}

template <typename S>
static void cross(S axb[3], const S a[3], const S b[3])
{
    const S b0 = b[0];
    const S b1 = b[1];
    const S b2 = b[2];
    axb[0] = a[1]*b2 - a[2]*b1;
    axb[1] = a[2]*b0 - a[0]*b2;
    axb[2] = a[0]*b1 - a[1]*b0;
}

template <typename S>
static void compute_exp_coefs(S &Xc, S &Yc, S &Zc,
                              S theta_sq, const liegroups::ExpCoefs<S> &coefs,
                              S lambda, S exp_neg_lambda)
{
    using namespace liegroups;
    
    const S lambda_sq = lambda * lambda;
    S D, alpha;
    if (theta_sq < Constants<S>::sqrt_epsilon()) {
        D = (S)(1.0/24) * (1 - theta_sq * (S)(1.0/30) * (1 - theta_sq * (S)(1.0/56)));
        if (lambda_sq > theta_sq) {
            alpha = lambda_sq / (lambda_sq + theta_sq);
        } else {
            alpha = 0;
        }
    } else {
        D = ((S)0.5 - coefs.B) / theta_sq;
        alpha = lambda_sq / (lambda_sq + theta_sq);
    }

    S beta, gamma;
    if (lambda_sq < Constants<S>::sqrt_epsilon()) {
        gamma = (S)(1.0/6) * (1 - lambda * (S)(0.25) * (1 - lambda * (S)(0.2) * (1 - lambda*(S)(1.0/6))));
        beta = (S)0.5 - lambda * gamma;
        Xc = 1 - lambda * beta;
    } else {
        const S inv_lambda = (S)1 / lambda;
        Xc = ((S)1 - exp_neg_lambda) * inv_lambda;
        beta = ((S)1 - Xc) * inv_lambda;
        gamma = ((S)0.5 - beta) * inv_lambda;
    }

    Yc = alpha * beta + (1 - alpha) * (coefs.B - lambda * coefs.C);
    Zc = alpha * gamma + (1 - alpha) * (coefs.C - lambda * D);
}

template <class S>
void liegroups::exp(Sim3<S> &X, const S uwl[7])
{
    const S *w = &uwl[3];
    const S lambda = uwl[6];
    const S theta_sq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
    const ExpCoefs<S> coefs(theta_sq);
    
    compute_exp_matrix3(X.rigid.R.R, coefs.cos_theta, coefs.A, coefs.B, w);
    X.scale = liegroups::exp(lambda);
    X.inv_scale = (S)1 / X.scale;

    S Xc, Yc, Zc;
    compute_exp_coefs(Xc, Yc, Zc, theta_sq, coefs, lambda, X.inv_scale);
    
    S wxu[3], wxwxu[3];
    cross(wxu, w, uwl);
    cross(wxwxu, w, wxu);

    for (int i=0; i<3; ++i)
        X.rigid.t[i] = Xc * uwl[i] + Yc * wxu[i] + Zc * wxwxu[i];
}

template void liegroups::exp<float>(Sim3<float> &, const float[7]);
template void liegroups::exp<double>(Sim3<double> &, const double[7]);

template <class S>
void liegroups::log(S uwl[7], const Sim3<S> &X)
{    
    uwl[6] = liegroups::ln(X.scale);

    S *w = &uwl[3];
    log(w, X.rigid.R);

    const S theta_sq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
    const ExpCoefs<S> coefs(theta_sq);
    
    S Xc, Yc, Zc;
    compute_exp_coefs(Xc, Yc, Zc, theta_sq, coefs, uwl[6], X.inv_scale);

    S V[3*3];
    compute_exp_matrix3(V, Xc - theta_sq*Zc, Yc, Zc, w);

    const S *const t = X.rigid.t;
    invert<3>(V, V);
    for (int i=0; i<3; ++i)
        uwl[i] = V[i*3]*t[0] + V[i*3+1]*t[1] + V[i*3+2]*t[2];
}

template void liegroups::log<float>(float[7], const Sim3<float> &);
template void liegroups::log<double>(double[7], const Sim3<double> &);

template <class S>
void liegroups::adjoint(S adj[7*7], const Sim3<S> &g)
{
    const S *const R = g.rigid.R.R;
    const S s = g.scale;
    const S sx = s * g.rigid.t[0];
    const S sy = s * g.rigid.t[1];
    const S sz = s * g.rigid.t[2];
    
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            adj[i*7 + j] = s * R[i*3 + j];
            adj[(i+3)*7 + j] = 0;
            adj[(i+3)*7 + j+3] = R[i*3 + j];
        }
        adj[3 + i] = sy * R[6 + i] - sz * R[3 + i];
        adj[10 + i] = sz * R[i] - sx * R[6 + i];
        adj[17 + i] = sx * R[3 + i] - sy * R[i];
        adj[(i+3)*7 + 6] = 0;
        adj[42 + i] = 0;
        adj[45 + i] = 0;
    }
    adj[6] = -sx;
    adj[13] = -sy;
    adj[20] = -sz;
    adj[48] = 1;
}

template void liegroups::adjoint<float>(float[7*7], const Sim3<float> &);
template void liegroups::adjoint<double>(double[7*7], const Sim3<double> &);

template <class S>
void liegroups::adjoint_multiply(S y[7], const Sim3<S> &g, const S x[7])
{
    const S *const R = g.rigid.R.R;
    const S *const t = g.rigid.t;
    const S s = g.scale;

    const S x3 = x[3], x4 = x[4], x5 = x[5];
    y[3] = R[0]*x3 + R[1]*x4 + R[2]*x5;
    y[4] = R[3]*x3 + R[4]*x4 + R[5]*x5;
    y[5] = R[6]*x3 + R[7]*x4 + R[8]*x5;
    const S x0 = x[0], x1 = x[1], x2 = x[2];
    y[0] = s * (R[0]*x0 + R[1]*x1 + R[2]*x2 + t[1]*y[5] - t[2]*y[4] - t[0]*x[6]);
    y[1] = s * (R[3]*x0 + R[4]*x1 + R[5]*x2 + t[2]*y[3] - t[0]*y[5] - t[1]*x[6]);
    y[2] = s * (R[6]*x0 + R[7]*x1 + R[8]*x2 + t[0]*y[4] - t[1]*y[3] - t[2]*x[6]);
    y[6] = x[6];
}

template void liegroups::adjoint_multiply<float>(float[7], const Sim3<float> &, const float[7]);
template void liegroups::adjoint_multiply<double>(double[7], const Sim3<double> &, const double[7]);

template <class S>
void liegroups::adjoint_T_multiply(S y[7], const Sim3<S> &g, const S x[7])
{
    const S *const R = g.rigid.R.R;
    const S *const t = g.rigid.t;
    const S s = g.scale;

    const S sx0 = s*x[0], sx1 = s*x[1], sx2 = s*x[2];
    y[0] = R[0]*sx0 + R[3]*sx1 + R[6]*sx2;
    y[1] = R[1]*sx0 + R[4]*sx1 + R[7]*sx2;
    y[2] = R[2]*sx0 + R[5]*sx1 + R[8]*sx2;
    const S a = x[3] - (t[1]*sx2 - t[2]*sx1);
    const S b = x[4] - (t[2]*sx0 - t[0]*sx2);
    const S c = x[5] - (t[0]*sx1 - t[1]*sx0);
    y[3] = R[0]*a + R[3]*b + R[6]*c;
    y[4] = R[1]*a + R[4]*b + R[7]*c;
    y[5] = R[2]*a + R[5]*b + R[8]*c;
    y[6] = x[6] - (t[0]*sx0 + t[1]*sx1 + t[2]*sx2);
}

template void liegroups::adjoint_T_multiply<float>(float[7], const Sim3<float> &, const float[7]);
template void liegroups::adjoint_T_multiply<double>(double[7], const Sim3<double> &, const double[7]);
