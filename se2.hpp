#pragma once

namespace liegroups {

    template <class S>
    struct SE2
    {
        S r[2];
        S t[2];

        static const SE2<S> identity;
        enum { DoF = 3, Dim = 2 };
        typedef S Scalar;

        // Output can alias input
        static void ad_multiply(S ada_b[3], const S a[3], const S b[3]);

        static void ad(S ada[3*3], const S a[3]);        
    };
   
    template <class S>
    void multiply(SE2<S> &ab, const SE2<S> &a, const SE2<S> &b);

    template <class S>
    SE2<S> operator*(const SE2<S> &a, const SE2<S> &b)
    {
        SE2<S> ab;
        multiply(ab, a, b);
        return ab;
    }
    
    template <class S>
    void multiply_a_binv(SE2<S> &abinv, const SE2<S> &a, const SE2<S> &b);    
    
    template <class S>
    void invert(SE2<S> &g);

    template <class S>
    SE2<S> inverse(const SE2<S> &g)
    {
        SE2<S> inv = g;
        invert(inv);
        return inv;
    }

    template <class S>
    void rectify(SE2<S> &g);
    
    template <class S, class X>
    void transform_point(X y[2], const SE2<S> &g, const X x[2]);

    template <class S, class X>
    void transform_point_by_inverse(X y[2], const SE2<S> &g, const X x[2]);
    
    template <class S>
    void exp(SE2<S> &X, const S x[3]);

    // Compute X = exp(x),
    //      dexp = diff(log(exp(x + d)*exp(-x)), d) at d = 0
    template <class S>
    void exp_diff(SE2<S> &X, S dexp[3*3], const S x[3]);
    
    template <class S>
    void log(S x[3], const SE2<S> &X);

    // Compute x = log(X),
    //      dlog = diff(log(exp(d) * X), d) at d = 0
    template <class S>
    void log_diff(S x[3], S dlog[3*3], const SE2<S> &X);

    template <class S>
    S SO2_log(S r00, S r01);
    
    template <class S>
    void adjoint(S adj[3*3], const SE2<S> &g);
    
    template <class S>
    void adjoint_multiply(S y[3], const SE2<S> &g, const S x[3]);

    template <class S>
    void adjoint_T_multiply(S y[3], const SE2<S> &g, const S x[3]);

}
