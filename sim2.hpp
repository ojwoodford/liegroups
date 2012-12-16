#pragma once

#include <liegroups/se2.hpp>

namespace liegroups {

    template <class S>
    struct Sim2
    {
        SE2<S> rigid;
        S scale, inv_scale;
        static const Sim2<S> identity;
        enum { DoF = 4, Dim = 2 };
        typedef S Scalar;
    };

    template <class S>
    void multiply(Sim2<S> &ab, const Sim2<S> &a, const Sim2<S> &b);

    template <class S>
    Sim2<S> operator*(const Sim2<S> &a, const Sim2<S> &b)
    {
        Sim2<S> ab;
        multiply(ab, a, b);
        return ab;
    }
    
    template <class S>
    void multiply_a_binv(Sim2<S> &abinv, const Sim2<S> &a, const Sim2<S> &b);    
    
    template <class S>
    void invert(Sim2<S> &g);

    template <class S>
    Sim2<S> inverse(const Sim2<S> &g)
    {
        Sim2<S> inv = g;
        invert(inv);
        return inv;
    }

    template <class S>
    void rectify(Sim2<S> &g);
    
    template <class S, class X>
    void transform_point(X y[2], const Sim2<S> &g, const X x[2]);

    template <class S, class X>
    void transform_point_by_inverse(X y[2], const Sim2<S> &g, const X x[2]);
    
    template <class S>
    void exp(Sim2<S> &X, const S x[4]);

    template <class S>
    void log(S x[4], const Sim2<S> &X);

    template <class S>
    void adjoint(S adj[4*4], const Sim2<S> &g);
    
    template <class S>
    void adjoint_multiply(S y[4], const Sim2<S> &g, const S x[4]);

    template <class S>
    void adjoint_T_multiply(S y[4], const Sim2<S> &g, const S x[4]);

}
