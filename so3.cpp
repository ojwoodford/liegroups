#include <liegroups/so3.hpp>
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

template <class S>
void liegroups::multiply(SO3<S> &ab, const SO3<S> &a, const SO3<S> &b)
{
    mat_mult<3,3,3>(ab.R, a.R, b.R);
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
    normalize3(&g.R[1]);

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
    S w00 = w[0]*w[0], w11 = w[1]*w[1], w22 = w[2]*w[2];
    S theta_sq = w00 + w11 + w22;
    S a, b;
    if (theta_sq < Constants<S>::sqrt_epsilon()) {
        a = (S)1 - theta_sq*((S)(1.0/6) + theta_sq*(S)(1.0/120));
        b = (S)0.5 - theta_sq*((S)(1.0/24) - theta_sq*(S)(1.0/720));
    } else {
        S theta = liegroups::sqrt(theta_sq);
        S inv_theta = (S)1 / theta;
        S c = liegroups::cos(theta);
        S s = liegroups::sin(theta);
        a = inv_theta * s;
        b = inv_theta * inv_theta * ((S)1 - c);
    }

    X.R[0] = (S)1 - b * (w11 + w22);
    X.R[4] = (S)1 - b * (w00 + w22);
    X.R[8] = (S)1 - b * (w00 + w11);

    S Bab = b*w[0]*w[1];
    S Bac = b*w[0]*w[2];
    S Bbc = b*w[1]*w[2];
    X.R[1] = Bab - a*w[2];
    X.R[3] = Bab + a*w[2];
    X.R[2] = Bac + a*w[1];
    X.R[6] = Bac - a*w[1];
    X.R[5] = Bbc - a*w[0];
    X.R[7] = Bbc + a*w[0];    
}

template void liegroups::exp<float>(SO3<float> &, const float[3]);
template void liegroups::exp<double>(SO3<double> &, const double[3]);

template <class S>
void liegroups::log(S w[3], const SO3<S> &X)
{
    w[0] = X.R[7] - X.R[5];
    w[1] = X.R[2] - X.R[6];
    w[2] = X.R[3] - X.R[1];    
    
    S tr = X.R[0] + X.R[4] + X.R[8];
    S ct = (S)0.5  * (tr - (S)1);

    if (ct > (S)0.99999) {
        S rsq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
        S f = (S)0.5 + rsq*((S)(1/48.0) + rsq*(S)(3/1280.0));
        w[0] *= f;
        w[1] *= f;
        w[2] *= f;
        return;
    }

    ct = liegroups::max((S)-1, ct);
    S theta = liegroups::acos(ct);
    
    if (ct > (S)-0.99) {
        S st = liegroups::sqrt((S)1 - ct*ct);
        S factor = (S)0.5 * theta / st;
        w[0] *= factor;
        w[1] *= factor;
        w[2] *= factor;
        return;
    }

    S theta_sq = theta*theta;
    S inv_B = theta_sq / ((S)1 - ct);

    S a = liegroups::sqrt(liegroups::max((S)0, inv_B*(X.R[0] - ct)));
    S b = liegroups::sqrt(liegroups::max((S)0, inv_B*(X.R[4] - ct)));
    S c = liegroups::sqrt(liegroups::max((S)0, inv_B*(X.R[8] - ct)));
    
    w[0] = w[0] < (S)0 ? -a : a;
    w[1] = w[1] < (S)0 ? -b : b;
    w[2] = w[2] < (S)0 ? -c : c;    
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
