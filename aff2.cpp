#include <liegroups/aff2.hpp>
#include <liegroups/scalar.hpp>
#include <liegroups/matrix.hpp>
#include <liegroups/matrix_impl.hpp>
#include <cmath>
#include <iostream>

template <> const liegroups::Aff2<float>
liegroups::Aff2<float>::identity = { {1.f, 0.f, 0.f, 1.f}, {0.f, 0.f} };

template <> const liegroups::Aff2<double>
liegroups::Aff2<double>::identity = { {1.0, 0.0, 0.0, 1.0}, {0.0, 0.0} };

template <class S>
void liegroups::multiply(Aff2<S> &ab, const Aff2<S> &a, const Aff2<S> &b)
{
    ab.A[0] = a.A[0] * b.A[0] + a.A[1] * b.A[2];
    ab.A[1] = a.A[0] * b.A[1] + a.A[1] * b.A[3];
    ab.A[2] = a.A[2] * b.A[0] + a.A[3] * b.A[2];
    ab.A[3] = a.A[2] * b.A[1] + a.A[3] * b.A[3];
    ab.t[0] = a.A[0] * b.t[0] + a.A[1] * b.t[1] + a.t[0];
    ab.t[1] = a.A[2] * b.t[0] + a.A[3] * b.t[1] + a.t[1];
}

template void liegroups::multiply<float>(Aff2<float>&, const Aff2<float>&, const Aff2<float> &);
template void liegroups::multiply<double>(Aff2<double>&, const Aff2<double>&, const Aff2<double> &);

template <class S>
void liegroups::multiply_a_binv(Aff2<S> &abinv, const Aff2<S> &a, const Aff2<S> &b)
{
    multiply(abinv, a, inverse(b));
}

template void liegroups::multiply_a_binv<float>(Aff2<float>&, const Aff2<float>&, const Aff2<float> &);
template void liegroups::multiply_a_binv<double>(Aff2<double>&, const Aff2<double>&, const Aff2<double> &);

template <class S>
void liegroups::invert(Aff2<S> &g)
{
    invert<2>(g.A, g.A);
    const S t0 = g.t[0];
    const S t1 = g.t[1];
    g.t[0] = -(g.A[0] * t0 + g.A[1] * t1);
    g.t[1] = -(g.A[2] * t0 + g.A[3] * t1);
}

template void liegroups::invert<float>(Aff2<float>&);
template void liegroups::invert<double>(Aff2<double>&);

template <class S, class X>
void liegroups::transform_point(X y[2], const Aff2<S> &g, const X x[2])
{
    const X x0 = x[0];
    const X x1 = x[1];
    y[0] = (X)(g.A[0]*x0 + g.A[1]*x1 + g.t[0]);
    y[1] = (X)(g.A[2]*x0 + g.A[3]*x1 + g.t[1]);
}

template void liegroups::transform_point<float,float>(float[2], const Aff2<float> &, const float[2]);
template void liegroups::transform_point<double,double>(double[2], const Aff2<double> &, const double[2]);
template void liegroups::transform_point<double,float>(float[2], const Aff2<double> &, const float[2]);
template void liegroups::transform_point<float,double>(double[2], const Aff2<float> &, const double[2]);

template <class S, class X>
void liegroups::transform_point_by_inverse(X y[2], const Aff2<S> &g, const X x[2])
{
    transform_point(y, inverse(g), x);
}

template void liegroups::transform_point_by_inverse<float,float>(float[2], const Aff2<float> &, const float[2]);
template void liegroups::transform_point_by_inverse<double,double>(double[2], const Aff2<double> &, const double[2]);
template void liegroups::transform_point_by_inverse<double,float>(float[2], const Aff2<double> &, const float[2]);
template void liegroups::transform_point_by_inverse<float,double>(double[2], const Aff2<float> &, const double[2]);

template <typename S>
static void to_alg(S h[3*3], const S v[6])
{
    h[0] = v[3] + v[4];
    h[1] = v[5] - v[2];
    h[2] = v[0];
    h[3] = v[2] + v[5];
    h[4] = v[3] - v[4];
    h[5] = v[1];
    h[6] = 0;
    h[7] = 0;
    h[8] = 0;
}

template <typename S>
static void from_alg(S v[6], const S h[3*3])
{
    v[0] = h[2];
    v[1] = h[5];
    v[2] = (S)0.5 * (h[3] - h[1]);
    v[3] = (S)0.5 * (h[0] + h[4]);
    v[4] = (S)0.5 * (h[0] - h[4]);
    v[5] = (S)0.5 * (h[1] + h[3]);
}

template <class S>
void liegroups::exp(Aff2<S> &X, const S x[6])
{
    S m[3*3];
    to_alg(m, x);

    S em[3*3];
    if (liegroups::expm<3>(em, m)) {
        X.A[0] = em[0];
        X.A[1] = em[1];
        X.t[0] = em[2];
        X.A[2] = em[3];
        X.A[3] = em[4];
        X.t[1] = em[5];
    } else {
        X = Aff2<S>::identity;
    }

}

template void liegroups::exp<float>(Aff2<float> &, const float[6]);
template void liegroups::exp<double>(Aff2<double> &, const double[6]);

template <class S>
bool liegroups::log(S x[6], const Aff2<S> &X)
{
    S m[3*3];
    const S em[3*3] = {X.A[0], X.A[1], X.t[0],
                       X.A[2], X.A[3], X.t[1],
                       0, 0, (S)1};
    if  (!liegroups::logm<3>(m, em)) {
        //std::cerr << "logm failed" << std::endl;
        return false;
    }

    from_alg(x, m);
    return true;
}

template bool liegroups::log<float>(float[6], const Aff2<float> &);
template bool liegroups::log<double>(double[6], const Aff2<double> &);

template <class S>
void liegroups::adjoint(S adj[6*6], const Aff2<S> &g)
{
    const S x = g.t[0], y = g.t[1];
    const S a = g.A[0], b = g.A[1], c = g.A[2], d = g.A[3];
    const S aa = a*a, ab = a*b, ac = a*c, ad = a*d;
    const S bb = b*b, bc = b*c, bd = b*d;
    const S cc = c*c, cd = c*d;
    const S dd = d*d;
    const S f = (S)1 / (ad - bc);
    const S fx = f*x, fy = f*y;
    const S half_f = (S)0.5 * f;
    
    adj[0]  = a;
    adj[1]  = b;
    adj[2]  = fy*(aa + bb) - fx*(ac + bd);
    adj[3]  = -x;
    adj[4]  = fy*2*ab - fx*(ad + bc);
    adj[5]  = fx*(ac - bd) - fy*(aa - bb);

    adj[6]  = c;
    adj[7]  = d;
    adj[8]  = fy*(ac + bd) - fx*(cc + dd);
    adj[9]  = -y;
    adj[10] = fy*(ad + bc) - fx*2*cd;
    adj[11] = fx*(cc - dd) - fy*(ac - bd);

    adj[12] = 0;
    adj[13] = 0;
    adj[14] = half_f * (aa + bb + cc + dd);
    adj[15] = 0;
    adj[16] = f * (ab + cd);
    adj[17] = half_f * (bb + dd - aa - cc);

    adj[18] = 0;
    adj[19] = 0;
    adj[20] = 0;
    adj[21] = 1;
    adj[22] = 0;
    adj[23] = 0;

    adj[24] = 0;
    adj[25] = 0;
    adj[26] = f * (ac + bd);
    adj[27] = 0;
    adj[28] = f * (ad + bc);
    adj[29] = f * (bd - ac);

    adj[30] = 0;
    adj[31] = 0;
    adj[32] = half_f * (cc + dd - aa - bb);
    adj[33] = 0;
    adj[34] = f * (cd - ab);
    adj[35] = half_f * (aa + dd - bb - cc);
}

template void liegroups::adjoint<float>(float[6*6], const Aff2<float> &);
template void liegroups::adjoint<double>(double[6*6], const Aff2<double> &);

template <class S>
void liegroups::adjoint_multiply(S y[6], const Aff2<S> &g, const S x[6])
{
    S adj[6*6];
    adjoint(adj, g);
    const S v[6] = {x[0], x[1], x[2], x[3], x[4], x[5]};
    const S *a = adj;
    for (int i=0; i<6; ++i) {
        y[i] = a[0]*v[0] + a[1]*v[1] + a[2]*v[2] + a[3]*v[3] + a[4]*v[4] + a[5]*v[5];
        a += 6;
    }
}

template void liegroups::adjoint_multiply<float>(float[4], const Aff2<float> &, const float[4]);
template void liegroups::adjoint_multiply<double>(double[4], const Aff2<double> &, const double[4]);

template <class S>
void liegroups::adjoint_T_multiply(S y[4], const Aff2<S> &g, const S x[4])
{
    S adj[6*6];
    adjoint(adj, g);
    const S v[6] = {x[0], x[1], x[2], x[3], x[4], x[5]};
    const S *a = adj;
    for (int i=0; i<6; ++i) {
        y[i] = a[0]*v[0] + a[6]*v[1] + a[12]*v[2] + a[18]*v[3] + a[24]*v[4] + a[30]*v[5];
        ++a;
    }
}

template void liegroups::adjoint_T_multiply<float>(float[3], const Aff2<float> &, const float[3]);
template void liegroups::adjoint_T_multiply<double>(double[3], const Aff2<double> &, const double[3]);
