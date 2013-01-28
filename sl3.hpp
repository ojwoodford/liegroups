#pragma once

namespace liegroups {

    template <class S>
    struct SL3
    {
        S H[3*3];

        static const SL3<S> identity;
        enum { DoF = 8, Dim = 3 };
        typedef S Scalar;
    };

    template <class S>
    void multiply(SL3<S> &ab, const SL3<S> &a, const SL3<S> &b);

    template <class S>
    SL3<S> operator*(const SL3<S> &a, const SL3<S> &b)
    {
        SL3<S> ab;
        multiply(ab, a, b);
        return ab;
    }
    
    template <class S>
    void multiply_a_binv(SL3<S> &abinv, const SL3<S> &a, const SL3<S> &b);    
    
    template <class S>
    void invert(SL3<S> &g);

    template <class S>
    SL3<S> inverse(const SL3<S> &g);

    template <class S>
    void rectify(SL3<S> &g);
    
    template <class S, class X>
    void transform_point(X y[3], const SL3<S> &g, const X x[3]);

    template <class S, class X>
    void transform_point_by_inverse(X y[3], const SL3<S> &g, const X x[3]);
    
    template <class S>
    void exp(SL3<S> &X, const S x[8]);

    template <class S>
    bool log(S x[8], const SL3<S> &X);
    
    template <class S>
    void adjoint(S adj[8*8], const SL3<S> &g);
    
    template <class S>
    void adjoint_multiply(S y[8], const SL3<S> &g, const S x[8]);

    template <class S>
    void adjoint_T_multiply(S y[8], const SL3<S> &g, const S x[8]);

}
