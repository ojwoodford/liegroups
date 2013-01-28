#include <liegroups/sim2.hpp>
#include <liegroups/scalar.hpp>
#include <cmath>
#include <iostream>

template <> const liegroups::Sim2<float>
liegroups::Sim2<float>::identity = { liegroups::SE2<float>::identity, 1.0f, 1.0f };

template <> const liegroups::Sim2<double>
liegroups::Sim2<double>::identity = { liegroups::SE2<double>::identity, 1.0, 1.0 };

template <class S>
void liegroups::multiply(Sim2<S> &ab, const Sim2<S> &a, const Sim2<S> &b)
{
    multiply(ab.rigid, a.rigid, b.rigid);
    S tf = b.inv_scale - (S)1;
    ab.rigid.t[0] += tf*a.rigid.t[0];
    ab.rigid.t[1] += tf*a.rigid.t[1];
    ab.scale = a.scale * b.scale;
    ab.inv_scale = a.inv_scale * b.inv_scale;
}

template void liegroups::multiply<float>(Sim2<float>&, const Sim2<float>&, const Sim2<float> &);
template void liegroups::multiply<double>(Sim2<double>&, const Sim2<double>&, const Sim2<double> &);

template <class S>
void liegroups::multiply_a_binv(Sim2<S> &abinv, const Sim2<S> &a, const Sim2<S> &b)
{
    multiply_a_binv(abinv.rigid, a.rigid, b.rigid);
    abinv.rigid.t[0] *= b.scale;
    abinv.rigid.t[1] *= b.scale;
    abinv.scale = a.scale * b.inv_scale;
    abinv.inv_scale = a.inv_scale * b.scale;
}

template void liegroups::multiply_a_binv<float>(Sim2<float>&, const Sim2<float>&, const Sim2<float> &);
template void liegroups::multiply_a_binv<double>(Sim2<double>&, const Sim2<double>&, const Sim2<double> &);

template <class S>
void liegroups::invert(Sim2<S> &g)
{
    invert(g.rigid);
    S s = g.scale;
    g.rigid.t[0] *= s;
    g.rigid.t[1] *= s;
    g.scale = g.inv_scale;
    g.inv_scale = s;
}

template void liegroups::invert<float>(Sim2<float>&);
template void liegroups::invert<double>(Sim2<double>&);

template <class S>
void liegroups::rectify(Sim2<S> &g)
{
    rectify(g.rigid);
    g.inv_scale = (S)1 / g.scale;
}

template void liegroups::rectify<float>(Sim2<float>&);
template void liegroups::rectify<double>(Sim2<double>&);

template <class S, class X>
void liegroups::transform_point(X y[2], const Sim2<S> &g, const X x[2])
{
    transform_point(y, g.rigid, x);
    y[0] *= g.scale;
    y[1] *= g.scale;
}

template void liegroups::transform_point<float,float>(float[2], const Sim2<float> &, const float[2]);
template void liegroups::transform_point<double,double>(double[2], const Sim2<double> &, const double[2]);
template void liegroups::transform_point<double,float>(float[2], const Sim2<double> &, const float[2]);
template void liegroups::transform_point<float,double>(double[2], const Sim2<float> &, const double[2]);

template <class S, class X>
void liegroups::transform_point_by_inverse(X y[2], const Sim2<S> &g, const X x[2])
{
    y[0] = g.inv_scale * x[0];
    y[1] = g.inv_scale * x[1];
    transform_point_by_inverse(y, g.rigid, y);
}

template void liegroups::transform_point_by_inverse<float,float>(float[2], const Sim2<float> &, const float[2]);
template void liegroups::transform_point_by_inverse<double,double>(double[2], const Sim2<double> &, const double[2]);
template void liegroups::transform_point_by_inverse<double,float>(float[2], const Sim2<double> &, const float[2]);
template void liegroups::transform_point_by_inverse<float,double>(double[2], const Sim2<float> &, const double[2]);

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

template <class S>
void liegroups::exp(Sim2<S> &X, const S x[4])
{
    S theta_sq = x[2]*x[2];
    S ct, st;
    if (theta_sq < Constants<S>::sqrt_epsilon()) {
        ct = (S)1 - theta_sq*(S)0.5*((S)1 - theta_sq*(S)(1/12.0));
        st = x[2]*((S)1 - theta_sq*(S)(1/6.0)*((S)1 - theta_sq*(S)(1/20.0)));
    } else {
        ct = liegroups::cos(x[2]);
        st = liegroups::sin(x[2]);
    }
    
    const S s = liegroups::exp(x[3]);    

    S a, b;
    compute_V(a, b, x[2], x[3], ct, st, s);

    X.rigid.r[0] = ct;
    X.rigid.r[1] = -st;
    X.rigid.t[0] = a*x[0] - b*x[1];
    X.rigid.t[1] = b*x[0] + a*x[1];
    X.scale = s;
    X.inv_scale = (S)1 / s;
}

template void liegroups::exp<float>(Sim2<float> &, const float[4]);
template void liegroups::exp<double>(Sim2<double> &, const double[4]);

template <class S>
void liegroups::log(S x[4], const Sim2<S> &X)
{
    x[3] = liegroups::ln(X.scale);
    x[2] = liegroups::SO2_log(X.rigid.r[0], X.rigid.r[1]);
    
    S a, b;    
    compute_V(a, b, x[2], x[3], X.rigid.r[0], -X.rigid.r[1], X.scale);
    
    S inv_det = (S)1 / (a*a + b*b);
    x[0] = inv_det*(a*X.rigid.t[0] + b*X.rigid.t[1]);
    x[1] = inv_det*(a*X.rigid.t[1] - b*X.rigid.t[0]);
}

template void liegroups::log<float>(float[4], const Sim2<float> &);
template void liegroups::log<double>(double[4], const Sim2<double> &);

template <class S>
void liegroups::adjoint(S adj[4*4], const Sim2<S> &g)
{
    adj[0] = g.scale * g.rigid.r[0];
    adj[1] = g.scale * g.rigid.r[1];
    adj[4] = -adj[1];
    adj[5] = adj[0];
    adj[2] = g.scale * g.rigid.t[1];
    adj[3] = g.scale * -g.rigid.t[0];
    adj[6] = adj[3];
    adj[7] = -adj[2];
    adj[8] = adj[9] = adj[12] = adj[13] = adj[11] = adj[14] = 0;
    adj[10] = adj[15] = (S)1;
}

template void liegroups::adjoint<float>(float[4*4], const Sim2<float> &);
template void liegroups::adjoint<double>(double[4*4], const Sim2<double> &);

template <class S>
void liegroups::adjoint_multiply(S y[4], const Sim2<S> &g, const S x[4])
{
    const S *r = g.rigid.r;
    const S *t = g.rigid.t;
    S y0 = g.scale * (r[0]*x[0] + r[1]*x[1] + t[1]*x[2] - t[0]*x[3]);
    S y1 = g.scale * (r[0]*x[1] - r[1]*x[0] - t[0]*x[2] - t[1]*x[3]);
    y[0] = y0;
    y[1] = y1;
    y[2] = x[2];
    y[3] = x[3];
}

template void liegroups::adjoint_multiply<float>(float[4], const Sim2<float> &, const float[4]);
template void liegroups::adjoint_multiply<double>(double[4], const Sim2<double> &, const double[4]);

template <class S>
void liegroups::adjoint_T_multiply(S y[4], const Sim2<S> &g, const S x[4])
{
    const S *r = g.rigid.r;
    const S *t = g.rigid.t;
    S y0 = g.scale * (r[0]*x[0] - r[1]*x[1]);
    S y1 = g.scale * (r[1]*x[0] + r[0]*x[1]);
    y[2] = x[2] + g.scale * (t[1]*x[0] - t[0]*x[1]);
    y[3] = x[3] - g.scale * (t[0]*x[0] + t[1]*x[1]);
    y[0] = y0;
    y[1] = y1;
}

template void liegroups::adjoint_T_multiply<float>(float[3], const Sim2<float> &, const float[3]);
template void liegroups::adjoint_T_multiply<double>(double[3], const Sim2<double> &, const double[3]);
