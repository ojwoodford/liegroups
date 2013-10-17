#pragma once

namespace liegroups {

    template <class S>
    struct Aff2
    {
        S A[2*2];
        S t[2];
        static const Aff2<S> identity;
        enum { DoF = 6, Dim = 2 };
        typedef S Scalar;
    };

    template <class S>
    void multiply(Aff2<S> &ab, const Aff2<S> &a, const Aff2<S> &b);

    template <class S>
    Aff2<S> operator*(const Aff2<S> &a, const Aff2<S> &b)
    {
        Aff2<S> ab;
        multiply(ab, a, b);
        return ab;
    }
    
    template <class S>
    void multiply_a_binv(Aff2<S> &abinv, const Aff2<S> &a, const Aff2<S> &b);    
    
    template <class S>
    void invert(Aff2<S> &g);

    template <class S>
    Aff2<S> inverse(const Aff2<S> &g)
    {
        Aff2<S> inv = g;
        invert(inv);
        return inv;
    }

    template <class S> void rectify(Aff2<S> &g) {}
    
    template <class S, class X>
    void transform_point(X y[2], const Aff2<S> &g, const X x[2]);

    template <class S, class X>
    void transform_point_by_inverse(X y[2], const Aff2<S> &g, const X x[2]);
    
    template <class S>
    void exp(Aff2<S> &X, const S x[6]);

    template <class S>
    bool log(S x[6], const Aff2<S> &X);

    template <class S>
    void adjoint(S adj[6*6], const Aff2<S> &g);
    
    template <class S>
    void adjoint_multiply(S y[6], const Aff2<S> &g, const S x[6]);

    template <class S>
    void adjoint_T_multiply(S y[6], const Aff2<S> &g, const S x[6]);

}
