#pragma once

#include <liegroups/so3.hpp>

namespace liegroups {

    template <class S>
    struct SE3
    {
        SO3<S> R;
        S t[3];
        
        static const SE3<S> identity;
        enum { DoF = 6, Dim = 3 };
        typedef S Scalar;

        // Output can alias input
        static void ad_multiply(S ada_b[6], const S a[6], const S b[6]);
        
        static void ad(S ada[6*6], const S a[]);        
    };

    template <class S>
    void multiply(SE3<S> &ab, const SE3<S> &a, const SE3<S> &b);

    template <class S>
    SE3<S> operator*(const SE3<S> &a, const SE3<S> &b)
    {
        SE3<S> ab;
        multiply(ab, a, b);
        return ab;
    }
    
    template <class S>
    void multiply_a_binv(SE3<S> &abinv, const SE3<S> &a, const SE3<S> &b);    
    
    template <class S>
    void invert(SE3<S> &g);

    template <class S>
    SE3<S> inverse(const SE3<S> &g);

    template <class S>
    void rectify(SE3<S> &g);
    
    template <class S, class X>
    void transform_point(X y[3], const SE3<S> &g, const X x[3]);

    template <class S, class X>
    void transform_point_by_inverse(X y[3], const SE3<S> &g, const X x[3]);
    
    template <class S>
    void exp(SE3<S> &X, const S x[6]);

    // Compute X = exp(x),
    //      dexp = diff(log(exp(x + d)*exp(-x)), d) at d = 0
    template <class S>
    void exp_diff(SE3<S> &X, S dexp[6*6], const S x[6]);

    template <class S>
    void log(S x[6], const SE3<S> &X);

    // Compute x = log(X),
    //      dlog = diff(log(exp(d) * X), d) at d = 0
    template <class S>
    void log_diff(S x[6], S dlog[6*6], const SE3<S> &X);
    
    template <class S>
    void adjoint(S adj[6*6], const SE3<S> &g);
    
    template <class S>
    void adjoint_multiply(S y[6], const SE3<S> &g, const S x[6]);

    template <class S>
    void adjoint_T_multiply(S y[6], const SE3<S> &g, const S x[6]);
    
}
