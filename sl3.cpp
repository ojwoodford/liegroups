#include <liegroups/sl3.hpp>
#include <liegroups/scalar.hpp>
#include <liegroups/matrix.hpp>
#include <cmath>
#include <iostream>

template <> const liegroups::SL3<float>
liegroups::SL3<float>::identity = { {1.f, 0.f, 0.f,
                                     0.f, 1.f, 0.f,
                                     0.f, 0.f, 1.f} };

template <> const liegroups::SL3<double>
liegroups::SL3<double>::identity = { {1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0} };

template <class S>
void liegroups::multiply(SL3<S> &ab, const SL3<S> &a, const SL3<S> &b)
{
    mat_mult_square<3>(ab.H, a.H, b.H);
}

template void liegroups::multiply<float>(SL3<float>&, const SL3<float>&, const SL3<float> &);
template void liegroups::multiply<double>(SL3<double>&, const SL3<double>&, const SL3<double> &);


template <typename S1, typename S2>
static S2 dot3(const S1 a[3], const S2 b[3])
{
    return (S2)(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

template <class S>
void liegroups::multiply_a_binv(SL3<S> &abinv, const SL3<S> &a, const SL3<S> &b)
{
    abinv = a * inverse(b);
}

template void liegroups::multiply_a_binv<float>(SL3<float>&, const SL3<float>&, const SL3<float> &);
template void liegroups::multiply_a_binv<double>(SL3<double>&, const SL3<double>&, const SL3<double> &);

template <class S>
liegroups::SL3<S> liegroups::inverse(const SL3<S> &g)
{
    SL3<S> ginv = g;
    invert(ginv);
    return ginv;
}

template liegroups::SL3<float> liegroups::inverse<float>(const SL3<float>&);
template liegroups::SL3<double> liegroups::inverse<double>(const SL3<double>&);

template <class S>
void liegroups::invert(SL3<S> &g)
{
    const S m00 = g.H[4]*g.H[8] - g.H[5]*g.H[7];
    const S m10 = g.H[1]*g.H[8] - g.H[2]*g.H[7];
    const S m20 = g.H[1]*g.H[5] - g.H[2]*g.H[4];
    const S m01 = g.H[3]*g.H[8] - g.H[5]*g.H[6];
    const S m02 = g.H[3]*g.H[7] - g.H[4]*g.H[6];
    const S m11 = g.H[0]*g.H[8] - g.H[2]*g.H[6];
    const S m12 = g.H[0]*g.H[7] - g.H[1]*g.H[6];
    const S m21 = g.H[0]*g.H[5] - g.H[2]*g.H[3];
    const S m22 = g.H[0]*g.H[4] - g.H[1]*g.H[3];

    g.H[0] =  m00;
    g.H[1] = -m10;
    g.H[2] =  m20;
    g.H[3] = -m01;
    g.H[4] =  m11;
    g.H[5] = -m21;
    g.H[6] =  m02;
    g.H[7] = -m12;
    g.H[8] =  m22;
}

template void liegroups::invert<float>(SL3<float>&);
template void liegroups::invert<double>(SL3<double>&);

template <class S>
void liegroups::rectify(SL3<S> &g)
{
    const S m00 = g.H[4]*g.H[8] - g.H[5]*g.H[7];
    const S m10 = g.H[1]*g.H[8] - g.H[2]*g.H[7];
    const S m20 = g.H[1]*g.H[5] - g.H[2]*g.H[4];

    const S det = g.H[0]*m00 - g.H[3]*m10 + g.H[6]*m20;
    const S factor = liegroups::pow(det, (S)(-1.0/3));
    for (int i=0; i<9; ++i)
        g.H[i] *= factor;
}

template void liegroups::rectify<float>(SL3<float>&);
template void liegroups::rectify<double>(SL3<double>&);

template <class S, class X>
void liegroups::transform_point(X y[3], const SL3<S> &g, const X x[3])
{
    S y0 = dot3(&g.H[0], x);
    S y1 = dot3(&g.H[3], x);
    S y2 = dot3(&g.H[6], x);
    y[0] = y0;
    y[1] = y1;
    y[2] = y2;
}

template void liegroups::transform_point<float,float>(float[3], const SL3<float> &, const float[3]);
template void liegroups::transform_point<double,double>(double[3], const SL3<double> &, const double[3]);
template void liegroups::transform_point<double,float>(float[3], const SL3<double> &, const float[3]);
template void liegroups::transform_point<float,double>(double[3], const SL3<float> &, const double[3]);

template <class S, class X>
void liegroups::transform_point_by_inverse(X y[3], const SL3<S> &g, const X x[3])
{
    SL3<S> ginv = inverse(g);
    transform_point(y, ginv, x);
}

template void liegroups::transform_point_by_inverse<float,float>(float[3], const SL3<float> &, const float[3]);
template void liegroups::transform_point_by_inverse<double,double>(double[3], const SL3<double> &, const double[3]);
template void liegroups::transform_point_by_inverse<double,float>(float[3], const SL3<double> &, const float[3]);
template void liegroups::transform_point_by_inverse<float,double>(double[3], const SL3<float> &, const double[3]);

template <typename S>
static void to_alg(S h[3*3], const S v[8])
{
    h[0] = v[3] + v[4];
    h[1] = v[5] - v[2];
    h[2] = v[0];
    h[3] = v[2] + v[5];
    h[4] = v[3] - v[4];
    h[5] = v[1];
    h[6] = v[6];
    h[7] = v[7];
    h[8] = (S)-2 * v[3];
}

template <typename S>
static void from_alg(S v[8], const S h[3*3])
{
    v[0] = h[2];
    v[1] = h[5];
    v[2] = (S)0.5 * (h[3] - h[1]);
    v[3] = (S)0.5 * (h[0] + h[4]);
    v[4] = (S)0.5 * (h[0] - h[4]);
    v[5] = (S)0.5 * (h[1] + h[3]);
    v[6] = h[6];
    v[7] = h[7];
}

template <class S>
void liegroups::exp(SL3<S> &X, const S v[8])
{
    S h[3*3];
    to_alg(h, v);

    if (!liegroups::expm<3>(X.H, h))
        X = SL3<S>::identity;
    rectify(X);
}

template void liegroups::exp<float>(SL3<float> &, const float[8]);
template void liegroups::exp<double>(SL3<double> &, const double[8]);

template <class S>
bool liegroups::log(S v[8], const SL3<S> &X)
{
    S h[3*3];
    if  (!liegroups::logm<3>(h, X.H)) {
        //std::cerr << "logm failed" << std::endl;
        return false;
    }

    from_alg(v, h);
    return true;
}

template bool liegroups::log<float>(float[8], const SL3<float> &);
template bool liegroups::log<double>(double[8], const SL3<double> &);

template <class S>
void liegroups::adjoint(S adj[8*8], const SL3<S> &g)
{    
    const SL3<S> ginv = inverse(g);
    const S* L = g.H;
    const S* R = ginv.H;
    const S m[8*8] = {
        L[0]*R[6], L[1]*R[6], L[1]*R[0]-L[0]*R[3], L[0]*R[0]+L[1]*R[3]-(S)2*L[2]*R[6], L[0]*R[0]-L[1]*R[3], L[0]*R[3]+L[1]*R[0], L[2]*R[0], L[2]*R[3],
        L[0]*R[7], L[1]*R[7], L[1]*R[1]-L[0]*R[4], L[0]*R[1]+L[1]*R[4]-(S)2*L[2]*R[7], L[0]*R[1]-L[1]*R[4], L[0]*R[4]+L[1]*R[1], L[2]*R[1], L[2]*R[4],
        L[0]*R[8], L[1]*R[8], L[1]*R[2]-L[0]*R[5], L[0]*R[2]+L[1]*R[5]-(S)2*L[2]*R[8], L[0]*R[2]-L[1]*R[5], L[0]*R[5]+L[1]*R[2], L[2]*R[2], L[2]*R[5],

        L[3]*R[6], L[4]*R[6], L[4]*R[0]-L[3]*R[3], L[3]*R[0]+L[4]*R[3]-(S)2*L[5]*R[6], L[3]*R[0]-L[4]*R[3], L[3]*R[3]+L[4]*R[0], L[5]*R[0], L[5]*R[3],
        L[3]*R[7], L[4]*R[7], L[4]*R[1]-L[3]*R[4], L[3]*R[1]+L[4]*R[4]-(S)2*L[5]*R[7], L[3]*R[1]-L[4]*R[4], L[3]*R[4]+L[4]*R[1], L[5]*R[1], L[5]*R[4],
        L[3]*R[8], L[4]*R[8], L[4]*R[2]-L[3]*R[5], L[3]*R[2]+L[4]*R[5]-(S)2*L[5]*R[8], L[3]*R[2]-L[4]*R[5], L[3]*R[5]+L[4]*R[2], L[5]*R[2], L[5]*R[5],

        L[6]*R[6], L[7]*R[6], L[7]*R[0]-L[6]*R[3], L[6]*R[0]+L[7]*R[3]-(S)2*L[8]*R[6], L[6]*R[0]-L[7]*R[3], L[6]*R[3]+L[7]*R[0], L[8]*R[0], L[8]*R[3],
        L[6]*R[7], L[7]*R[7], L[7]*R[1]-L[6]*R[4], L[6]*R[1]+L[7]*R[4]-(S)2*L[8]*R[7], L[6]*R[1]-L[7]*R[4], L[6]*R[4]+L[7]*R[1], L[8]*R[1], L[8]*R[4],
    };        

    for (int i=0; i<8; ++i) {
        adj[ 0+i] = m[16+i];
        adj[ 8+i] = m[40+i];
        adj[16+i] = (S)0.5*(m[24+i] - m[ 8+i]);
        adj[24+i] = (S)0.5*(m[   i] + m[32+i]);
        adj[32+i] = (S)0.5*(m[   i] - m[32+i]);
        adj[40+i] = (S)0.5*(m[ 8+i] + m[24+i]);
        adj[48+i] = m[48+i];
        adj[56+i] = m[56+i];
    }
}

template void liegroups::adjoint<float>(float[8*8], const SL3<float> &);
template void liegroups::adjoint<double>(double[8*8], const SL3<double> &);

template <class S>
void liegroups::adjoint_multiply(S y[8], const SL3<S> &g, const S x[8])
{
    S h[3*3];
    to_alg(h, x);

    const SL3<S> ginv = inverse(g);

    // Compute H*alg(h)*H^-1
    S gh[3*3];
    mat_mult_square<3>(gh, g.H, h);
    mat_mult_square<3>(h, gh, ginv.H);

    from_alg(y, h);
}

template void liegroups::adjoint_multiply<float>(float[8], const SL3<float> &, const float[8]);
template void liegroups::adjoint_multiply<double>(double[8], const SL3<double> &, const double[8]);

template <class S>
void liegroups::adjoint_T_multiply(S y[8], const SL3<S> &g, const S x[8])
{
    S adj[8*8];
    adjoint(adj, g);

    const S v[8] = {x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]};
    for (int i=0; i<8; ++i) {
        S sum = 0;
        for (int j=0; j<8; ++j) {
            sum += adj[i + j*8] * v[j];
        }
        y[i] = sum;
    }
}

template void liegroups::adjoint_T_multiply<float>(float[8], const SL3<float> &, const float[8]);
template void liegroups::adjoint_T_multiply<double>(double[8], const SL3<double> &, const double[8]);
