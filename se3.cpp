#include <liegroups/se3.hpp>
#include <liegroups/scalar.hpp>
#include <liegroups/matrix.hpp>
#include <cmath>
#include <liegroups/exp_coefs.hpp>
#include <liegroups/exp_helpers.hpp>

using namespace liegroups;

template <> const SE3<float>
SE3<float>::identity = { SO3<float>::identity, {0.f, 0.f, 0.f} };

template <> const SE3<double>
SE3<double>::identity = { SO3<double>::identity, {0.0, 0.0, 0.0} };

template <typename S1, typename S2>
static S2 dot3(const S1 a[3], const S2 b[3])
{
    return (S2)(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

template <class S>
void liegroups::multiply(SE3<S> &ab, const SE3<S> &a, const SE3<S> &b)
{
    ab.t[0] = a.t[0] + dot3(&a.R.R[0], b.t);
    ab.t[1] = a.t[1] + dot3(&a.R.R[3], b.t);
    ab.t[2] = a.t[2] + dot3(&a.R.R[6], b.t);
    multiply(ab.R, a.R, b.R);
}

template void liegroups::multiply<float>(SE3<float>&, const SE3<float>&, const SE3<float> &);
template void liegroups::multiply<double>(SE3<double>&, const SE3<double>&, const SE3<double> &);

template <class S>
void liegroups::multiply_a_binv(SE3<S> &abinv, const SE3<S> &a, const SE3<S> &b)
{
    multiply_a_binv(abinv.R, a.R, b.R);    
    abinv.t[0] = a.t[0] - dot3(&abinv.R.R[0], b.t);
    abinv.t[1] = a.t[1] - dot3(&abinv.R.R[3], b.t);
    abinv.t[2] = a.t[2] - dot3(&abinv.R.R[6], b.t);
}

template void liegroups::multiply_a_binv<float>(SE3<float>&, const SE3<float>&, const SE3<float> &);
template void liegroups::multiply_a_binv<double>(SE3<double>&, const SE3<double>&, const SE3<double> &);

template <class S>
liegroups::SE3<S> liegroups::inverse(const SE3<S> &g)
{
    SE3<S> ginv;
    ginv.R = inverse(g.R);
    ginv.t[0] = -dot3(&ginv.R.R[0], g.t);
    ginv.t[1] = -dot3(&ginv.R.R[3], g.t);
    ginv.t[2] = -dot3(&ginv.R.R[6], g.t);
    return ginv;
}

template liegroups::SE3<float> liegroups::inverse<float>(const SE3<float>&);
template liegroups::SE3<double> liegroups::inverse<double>(const SE3<double>&);

template <class S>
void liegroups::invert(SE3<S> &g)
{
    g = inverse(g);
}

template void liegroups::invert<float>(SE3<float>&);
template void liegroups::invert<double>(SE3<double>&);

template <class S>
void liegroups::rectify(SE3<S> &g)
{
    rectify(g.R);
}

template void liegroups::rectify<float>(SE3<float>&);
template void liegroups::rectify<double>(SE3<double>&);

template <class S, class X>
void liegroups::transform_point(X y[3], const SE3<S> &g, const X x[3])
{
    transform_point(y, g.R, x);
    y[0] += (X)g.t[0];
    y[1] += (X)g.t[1];
    y[2] += (X)g.t[2];
}

template void liegroups::transform_point<float,float>(float[3], const SE3<float> &, const float[3]);
template void liegroups::transform_point<double,double>(double[3], const SE3<double> &, const double[3]);
template void liegroups::transform_point<double,float>(float[3], const SE3<double> &, const float[3]);
template void liegroups::transform_point<float,double>(double[3], const SE3<float> &, const double[3]);

template <class S, class X>
void liegroups::transform_point_by_inverse(X y[3], const SE3<S> &g, const X x[3])
{
    y[0] = x[0] - g.t[0];
    y[1] = x[1] - g.t[1];
    y[2] = x[2] - g.t[2];
    transform_point_by_inverse(y, g.R, y);
}

template void liegroups::transform_point_by_inverse<float,float>(float[3], const SE3<float> &, const float[3]);
template void liegroups::transform_point_by_inverse<double,double>(double[3], const SE3<double> &, const double[3]);
template void liegroups::transform_point_by_inverse<double,float>(float[3], const SE3<double> &, const float[3]);
template void liegroups::transform_point_by_inverse<float,double>(double[3], const SE3<float> &, const double[3]);

template <typename S>
static void cross(S axb[3], const S a[3], const S b[3])
{
    axb[0] = a[1]*b[2] - a[2]*b[1];
    axb[1] = a[2]*b[0] - a[0]*b[2];
    axb[2] = a[0]*b[1] - a[1]*b[0];
}

    
template <class S>
void liegroups::exp(SE3<S> &X, const S uw[6])
{
    const S *w = &uw[3];
    const S theta_sq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
    const ExpCoefs<S> coefs(theta_sq);
    
    compute_exp_matrix3(X.R.R, coefs.cos_theta, coefs.A, coefs.B, w);

    S wxu[3], wxwxu[3];
    cross(wxu, w, uw);
    cross(wxwxu, w, wxu);

    for (int i=0; i<3; ++i)
        X.t[i] = uw[i] + coefs.B*wxu[i] + coefs.C*wxwxu[i];
}

template void liegroups::exp<float>(SE3<float> &, const float[3]);
template void liegroups::exp<double>(SE3<double> &, const double[3]);

template <class S>
void liegroups::log(S uw[6], const SE3<S> &X)
{
    S *w = &uw[3];
    log(w, X.R);

    const S theta_sq = dot3(w,w);
    const ExpCoefs<S> coefs(theta_sq);
    const S a = coefs.A;
    const S b = coefs.B;
    const S c = coefs.C;

    S d;
    if (theta_sq < Constants<S>::epsilon()*25) {
        d = (S)(1.0/12) + theta_sq*((S)(1.0/720) + theta_sq*(S)(1.0/30240));
    } else if (theta_sq > (S)9) {
        d = (b - (S)0.5*a) / (b*theta_sq);
    } else {
        d = (b*(S)0.5 - c) / a;
    }

    S wxt[3], wxwxt[3];
    cross(wxt, w, X.t);
    cross(wxwxt, w, wxt);

    for (int i=0; i<3; ++i) {
        uw[i] = X.t[i] - (S)0.5*wxt[i] + d*wxwxt[i];
    }
}

template void liegroups::log<float>(float[6], const SE3<float> &);
template void liegroups::log<double>(double[6], const SE3<double> &);

template <class S>
void liegroups::adjoint(S adj[6*6], const SE3<S> &g)
{    
    for (int i=0; i<3; ++i) {
        S a = g.R.R[i], b = g.R.R[i+3], c = g.R.R[i+6];
        adj[i] = adj[i+21] = a;
        adj[i+6] = adj[i+27] = b;
        adj[i+12] = adj[i+33] = c;

        adj[i+3] = g.t[1]*c - g.t[2]*b;
        adj[i+9] = g.t[2]*a - g.t[0]*c;
        adj[i+15] = g.t[0]*b - g.t[1]*a;
        
        adj[i+18] = (S)0;
        adj[i+24] = (S)0;
        adj[i+30] = (S)0;
    }
}

template void liegroups::adjoint<float>(float[6*6], const SE3<float> &);
template void liegroups::adjoint<double>(double[6*6], const SE3<double> &);

template <class S>
void liegroups::adjoint_multiply(S y[6], const SE3<S> &g, const S x[6])
{
    transform_point(&y[0], g.R, &x[0]);
    transform_point(&y[3], g.R, &x[3]);
    
    y[0] += g.t[1]*y[5] - g.t[2]*y[4];
    y[1] += g.t[2]*y[3] - g.t[0]*y[5];    
    y[2] += g.t[0]*y[4] - g.t[1]*y[3];
}

template void liegroups::adjoint_multiply<float>(float[6], const SE3<float> &, const float[6]);
template void liegroups::adjoint_multiply<double>(double[6], const SE3<double> &, const double[6]);

template <class S>
void liegroups::adjoint_T_multiply(S y[6], const SE3<S> &g, const S x[6])
{
    y[3] = x[3] - (g.t[1]*x[2] - g.t[2]*x[1]);
    y[4] = x[4] - (g.t[2]*x[0] - g.t[0]*x[2]);
    y[5] = x[5] - (g.t[0]*x[1] - g.t[1]*x[0]);
    transform_point_by_inverse(&y[3], g.R, &y[3]);
    transform_point_by_inverse(&y[0], g.R, &x[0]);
}

template void liegroups::adjoint_T_multiply<float>(float[6], const SE3<float> &, const float[6]);
template void liegroups::adjoint_T_multiply<double>(double[6], const SE3<double> &, const double[6]);


template <class S>
void liegroups::exp_diff(SE3<S> &exp_uw, S dexp[6*6], const S uw[6])
{
    const S *const u = &uw[0];
    const S *const w = &uw[3];
    
    const S w00 = w[0]*w[0];
    const S w11 = w[1]*w[1];
    const S w22 = w[2]*w[2];
    const S theta_sq = w00 + w11 + w22;
    DiffExpCoefs<S> coefs(theta_sq);

    // exp_uw.R = I + Awx + Bwx^2
    compute_exp_matrix3(exp_uw.R.R, coefs.cos_theta, coefs.A, coefs.B, w);
    
    // dexp = [ V  X ]
    //        [ 0  V ]
    //
    // V = A*I + B*w_x + C*ww'
    //
    compute_exp_matrix3(&dexp[0], &dexp[6], &dexp[12], coefs.A, coefs.B, coefs.C, w);
    
    // exp_uw.t = V * u
    exp_uw.t[0] = dot3(&dexp[ 0], u);
    exp_uw.t[1] = dot3(&dexp[ 6], u);
    exp_uw.t[2] = dot3(&dexp[12], u);
    
    // X = B*u_x + C*(uw' + wu') + (w'u) * (A1*I + B1*w_x + C1*ww')
    {
        const S uw0 = u[0]*w[0];
        const S uw1 = u[1]*w[1];
        const S uw2 = u[2]*w[2];
        const S wtu = uw0 + uw1 + uw2;

        // Diagonal
        dexp[ 3] = coefs.C*2*uw0 + wtu*(coefs.A1 + coefs.C1 * w00);
        dexp[10] = coefs.C*2*uw1 + wtu*(coefs.A1 + coefs.C1 * w11);
        dexp[17] = coefs.C*2*uw2 + wtu*(coefs.A1 + coefs.C1 * w22);

        // Symmetric
        const S s01 = coefs.C*(u[0]*w[1] + u[1]*w[0]);
        const S s02 = coefs.C*(u[0]*w[2] + u[2]*w[0]);
        const S s12 = coefs.C*(u[1]*w[2] + u[2]*w[1]);
        
        // Anti-symmetric
        const S a10 = coefs.B*u[2] + wtu*coefs.B1*w[2];
        const S a02 = coefs.B*u[1] + wtu*coefs.B1*w[1];
        const S a21 = coefs.B*u[0] + wtu*coefs.B1*w[0];
        
        dexp[ 4] = s01 - a10;
        dexp[ 9] = s01 + a10;
        dexp[ 5] = s02 + a02;
        dexp[15] = s02 - a02;
        dexp[11] = s12 - a21;
        dexp[16] = s12 + a21;
    }

    // Zero lower left, and copy V to lower right
    for (int i=0; i<3; ++i) {
        dexp[18+i] = 0;
        dexp[24+i] = 0;
        dexp[30+i] = 0;
        dexp[21+i] = dexp[i];
        dexp[27+i] = dexp[6+i];
        dexp[33+i] = dexp[12+i];
    }
}

template void liegroups::exp_diff<float>(SE3<float>&, float[6*6], const float [6]);
template void liegroups::exp_diff<double>(SE3<double>&, double[6*6], const double [6]);
