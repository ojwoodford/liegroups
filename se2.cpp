#include <liegroups/se2.hpp>
#include <liegroups/scalar.hpp>
#include <liegroups/exp_coefs.hpp>
#include <cmath>

template <> const liegroups::SE2<float>
liegroups::SE2<float>::identity = { {1.f, 0.f}, {0.f, 0.f} };

template <> const liegroups::SE2<double>
liegroups::SE2<double>::identity = { {1.0, 0.0}, {0.0, 0.0} };

template <typename S>
void liegroups::SE2<S>::ad_multiply(S ada_b[3], const S a[3], const S b[3])
{
    const S b0 = b[0];
    ada_b[0] = a[1]*b[2] - a[2]*b[1];
    ada_b[1] = a[2]*b0   - a[0]*b[2];
    ada_b[2] = 0;
}

template void liegroups::SE2<float>::ad_multiply(float[3], const float[3], const float[3]);
template void liegroups::SE2<double>::ad_multiply(double[3], const double[3], const double[3]);

template <typename S>
void liegroups::SE2<S>::ad(S ada[3*3], const S a[3])
{
    ada[0] = ada[4] = 0;
    ada[1] = -a[2];
    ada[2] = a[1];
    ada[3] = a[2];
    ada[5] = -a[0];
    ada[6] = ada[7] = ada[8] = 0;    
}

template void liegroups::SE2<float>::ad(float[3*3], const float a[3]);
template void liegroups::SE2<double>::ad(double[3*3], const double a[3]);

template <class S>
void liegroups::multiply(SE2<S> &ab, const SE2<S> &a, const SE2<S> &b)
{
    S r0 = a.r[0]*b.r[0] - a.r[1]*b.r[1];
    S r1 = a.r[0]*b.r[1] + a.r[1]*b.r[0];
    S t0 = a.r[0]*b.t[0] + a.r[1]*b.t[1] + a.t[0];
    S t1 = a.r[0]*b.t[1] - a.r[1]*b.t[0] + a.t[1];
    ab.r[0] = r0;
    ab.r[1] = r1;
    ab.t[0] = t0;
    ab.t[1] = t1;
}

template void liegroups::multiply<float>(SE2<float>&, const SE2<float>&, const SE2<float> &);
template void liegroups::multiply<double>(SE2<double>&, const SE2<double>&, const SE2<double> &);

template <class S>
void liegroups::multiply_a_binv(SE2<S> &abinv, const SE2<S> &a, const SE2<S> &b)
{    
    S r0 = a.r[0]*b.r[0] + a.r[1]*b.r[1];
    S r1 = -a.r[0]*b.r[1] + a.r[1]*b.r[0];
    S t0 = a.t[0] - r0*b.t[0] - r1*b.t[1];
    S t1 = a.t[1] - r0*b.t[1] + r1*b.t[0];
    abinv.r[0] = r0;
    abinv.r[1] = r1;
    abinv.t[0] = t0;
    abinv.t[1] = t1;
}

template void liegroups::multiply_a_binv<float>(SE2<float>&, const SE2<float>&, const SE2<float> &);
template void liegroups::multiply_a_binv<double>(SE2<double>&, const SE2<double>&, const SE2<double> &);

template <class S>
void liegroups::invert(SE2<S> &g)
{
    g.r[1] = -g.r[1];
    S t0 = -(g.r[0]*g.t[0] + g.r[1]*g.t[1]);
    g.t[1] = g.r[1]*g.t[0] - g.r[0]*g.t[1];
    g.t[0] = t0;
}

template void liegroups::invert<float>(SE2<float>&);
template void liegroups::invert<double>(SE2<double>&);

template <class S>
void liegroups::rectify(SE2<S> &g)
{
    S rr = g.r[0]*g.r[0] + g.r[1]*g.r[1];
    S factor = (S)1 / (S)::sqrt(rr);
    g.r[0] *= factor;
    g.r[1] *= factor;
}

template void liegroups::rectify<float>(SE2<float>&);
template void liegroups::rectify<double>(SE2<double>&);

template <class S, class X>
void liegroups::transform_point(X y[2], const SE2<S> &g, const X x[2])
{
    X y0 = (X)(g.r[0]*x[0] + g.r[1]*x[1] + g.t[0]);
    X y1 = (X)(-g.r[1]*x[0] + g.r[0]*x[1] + g.t[1]);
    y[0] = y0;
    y[1] = y1;
}

template void liegroups::transform_point<float,float>(float[2], const SE2<float> &, const float[2]);
template void liegroups::transform_point<double,double>(double[2], const SE2<double> &, const double[2]);
template void liegroups::transform_point<double,float>(float[2], const SE2<double> &, const float[2]);
template void liegroups::transform_point<float,double>(double[2], const SE2<float> &, const double[2]);

template <class S, class X>
void liegroups::transform_point_by_inverse(X y[2], const SE2<S> &g, const X x[2])
{
    S x0 = (S)x[0] - g.t[0];
    S x1 = (S)x[1] - g.t[1];
    y[0] = (X)(g.r[0]*x0 - g.r[1]*x1);
    y[1] = (X)(g.r[1]*x0 + g.r[0]*x1);
}

template void liegroups::transform_point_by_inverse<float,float>(float[2], const SE2<float> &, const float[2]);
template void liegroups::transform_point_by_inverse<double,double>(double[2], const SE2<double> &, const double[2]);
template void liegroups::transform_point_by_inverse<double,float>(float[2], const SE2<double> &, const float[2]);
template void liegroups::transform_point_by_inverse<float,double>(double[2], const SE2<float> &, const double[2]);

template <class S>
void liegroups::exp(SE2<S> &X, const S x[3])
{
    const S theta_sq = x[2]*x[2];
    const ExpCoefs<S> coefs(theta_sq);
    const S Bt = coefs.B*x[2];
    X.r[0] = coefs.cos_theta;
    X.r[1] = -x[2]*coefs.A;
    X.t[0] = coefs.A*x[0] - Bt*x[1];
    X.t[1] = Bt*x[0] + coefs.A*x[1];
}

template void liegroups::exp<float>(SE2<float> &, const float[3]);
template void liegroups::exp<double>(SE2<double> &, const double[3]);

template <class S>
void liegroups::exp_diff(SE2<S> &X, S dexp[3*3], const S x[3])
{
    const S theta_sq = x[2]*x[2];
    const ExpCoefs<S> coefs(theta_sq);
    const S Bt = coefs.B*x[2];
    X.r[0] = coefs.cos_theta;
    X.r[1] = -x[2]*coefs.A;
    X.t[0] = coefs.A*x[0] - Bt*x[1];
    X.t[1] = Bt*x[0] + coefs.A*x[1];

    dexp[0] = coefs.A;
    dexp[4] = coefs.A;
    dexp[6] = 0;
    dexp[7] = 0;
    dexp[8] = (S)1;
    
    dexp[1] = -Bt;
    dexp[3] = Bt;
    
    const S Ct = coefs.C*x[2];
    dexp[2] = coefs.B * x[1] + Ct * x[0];
    dexp[5] = Ct * x[1] - coefs.B * x[0];
}

template void liegroups::exp_diff<float>(SE2<float>&, float[3*3], const float[3]);
template void liegroups::exp_diff<double>(SE2<double>&, double[3*3], const double[3]);

template <class S>
S liegroups::SO2_log(S r00, S r01)
{
    return liegroups::atan2(-r01, r00);
}

template float liegroups::SO2_log(float r00, float r01);
template double liegroups::SO2_log(double r00, double r01);

template <class S>
void liegroups::log(S x[3], const SE2<S> &X)
{
    const S theta = SO2_log(X.r[0], X.r[1]);
    const S theta_sq = theta * theta;
    const ExpCoefs<S> coefs(theta_sq);
    const S f = (S)0.5 * coefs.A / coefs.B;
    const S ht = (S)0.5 * theta;

    x[0] = f * X.t[0] + ht * X.t[1];
    x[1] = f * X.t[1] - ht * X.t[0];
    x[2] = theta;
}

template void liegroups::log<float>(float[3], const SE2<float> &);
template void liegroups::log<double>(double[3], const SE2<double> &);

template <class S>
void liegroups::log_diff(S x[3], S dlog[3*3], const SE2<S> &X)
{
    const S theta = SO2_log(X.r[0], X.r[1]);
    const S theta_sq = theta * theta;
    const ExpCoefs<S> coefs(theta_sq);
    const S f = (S)0.5 * coefs.A / coefs.B;    
    const S ht = (S)0.5 * theta;
    
    x[0] = f * X.t[0] + ht * X.t[1];
    x[1] = f * X.t[1] - ht * X.t[0];
    x[2] = theta;

    const S Ct = coefs.C * theta;
    const S ex = Ct * x[0] + coefs.B * x[1];
    const S ey = Ct * x[1] - coefs.B * x[0];
    dlog[0] = dlog[4] = f;
    dlog[1] = ht;
    dlog[3] = -ht;
    dlog[2] = -(f * ex + ht * ey);
    dlog[5] = ht * ex - f * ey;
    dlog[6] = dlog[7] = 0;
    dlog[8] = 1;
}

template void liegroups::log_diff<float>(float[3], float[3*3], const SE2<float> &);
template void liegroups::log_diff<double>(double[3], double[3*3], const SE2<double> &);

template <class S>
void liegroups::adjoint(S adj[3*3], const SE2<S> &g)
{
    adj[0] = adj[4] = g.r[0];
    adj[1] = g.r[1];
    adj[3] = -g.r[1];
    adj[2] = g.t[1];
    adj[5] = -g.t[0];
    adj[6] = adj[7] = (S)0;
    adj[8] = (S)1;
}

template void liegroups::adjoint<float>(float[3*3], const SE2<float> &);
template void liegroups::adjoint<double>(double[3*3], const SE2<double> &);

template <class S>
void liegroups::adjoint_multiply(S y[3], const SE2<S> &g, const S x[3])
{
    S y0 = g.r[0]*x[0] + g.r[1]*x[1] + g.t[1]*x[2];
    S y1 = -g.r[1]*x[0] + g.r[0]*x[1] - g.t[0]*x[2];
    y[0] = y0;
    y[1] = y1;
    y[2] = x[2];
}

template void liegroups::adjoint_multiply<float>(float[3], const SE2<float> &, const float[3]);
template void liegroups::adjoint_multiply<double>(double[3], const SE2<double> &, const double[3]);

template <class S>
void liegroups::adjoint_T_multiply(S y[3], const SE2<S> &g, const S x[3])
{
    S y0 = g.r[0]*x[0] - g.r[1]*x[1];
    S y1 = g.r[0]*x[1] + g.r[1]*x[0];
    S y2 = g.t[1]*x[0] - g.t[0]*x[1] + x[2];
    y[0] = y0;
    y[1] = y1;
    y[2] = y2;    
}

template void liegroups::adjoint_T_multiply<float>(float[3], const SE2<float> &, const float[3]);
template void liegroups::adjoint_T_multiply<double>(double[3], const SE2<double> &, const double[3]);
