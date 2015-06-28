#pragma once

#include <liegroups/se3.hpp>

namespace liegroups {

    template <class S>
    struct Sim3
    {
        SE3<S> rigid;
        S scale, inv_scale;
        static const Sim3<S> identity;
        enum { DoF = 7, Dim = 3 };
        typedef S Scalar;
    };

    template <class S>
    void multiply(Sim3<S> &ab, const Sim3<S> &a, const Sim3<S> &b);

    template <class S>
    Sim3<S> operator*(const Sim3<S> &a, const Sim3<S> &b)
    {
        Sim3<S> ab;
        multiply(ab, a, b);
        return ab;
    }
    
    template <class S>
    void multiply_a_binv(Sim3<S> &abinv, const Sim3<S> &a, const Sim3<S> &b);    
    
    template <class S>
    void invert(Sim3<S> &g);

    template <class S>
    Sim3<S> inverse(const Sim3<S> &g)
    {
        Sim3<S> inv = g;
        invert(inv);
        return inv;
    }

    template <class S>
    void rectify(Sim3<S> &g);
    
    template <class S, class X>
    void transform_point(X y[3], const Sim3<S> &g, const X x[3]);

    template <class S, class X>
    void transform_point_by_inverse(X y[3], const Sim3<S> &g, const X x[3]);
    
    template <class S>
    void exp(Sim3<S> &X, const S x[7]);

    template <class S>
    void log(S x[7], const Sim3<S> &X);

    template <class S>
    void adjoint(S adj[7*7], const Sim3<S> &g);
    
    template <class S>
    void adjoint_multiply(S y[7], const Sim3<S> &g, const S x[7]);

    template <class S>
    void adjoint_T_multiply(S y[7], const Sim3<S> &g, const S x[7]);

}
