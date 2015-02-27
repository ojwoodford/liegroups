#pragma once

namespace liegroups {

    template <class S>
    struct SO3
    {
        S R[3*3];

        static const SO3<S> identity;
        enum { DoF = 3, Dim = 3 };
        typedef S Scalar;

        // Output can alias input
        static void ad_multiply(S ada_b[3], const S a[3], const S b[3]);
        
        static void ad(S ada[3*3], const S a[]);
    };

    template <class S>
    void multiply(SO3<S> &ab, const SO3<S> &a, const SO3<S> &b);

    template <class S>
    SO3<S> operator*(const SO3<S> &a, const SO3<S> &b)
    {
        SO3<S> ab;
        multiply(ab, a, b);
        return ab;
    }
    
    template <class S>
    void multiply_a_binv(SO3<S> &abinv, const SO3<S> &a, const SO3<S> &b);    
    
    template <class S>
    void invert(SO3<S> &g);

    template <class S>
    SO3<S> inverse(const SO3<S> &g);

    template <class S>
    void rectify(SO3<S> &g);
    
    template <class S, class X>
    void transform_point(X y[3], const SO3<S> &g, const X x[3]);

    template <class S, class X>
    void transform_point_by_inverse(X y[3], const SO3<S> &g, const X x[3]);
    
    template <class S>
    void exp(SO3<S> &X, const S x[3]);

    // Compute X = exp(x),
    //      dexp = diff(log(exp(x + d)*exp(-x)), d) at d = 0
    template <class S>
    void exp_diff(SO3<S> &X, S dexp[3*3], const S x[3]);

    template <class S>
    void log(S x[3], const SO3<S> &X);

    // Compute x = log(X),
    //      dlog = diff(log(exp(d) * X), d) at d = 0
    template <class S>
    void log_diff(S x[3], S dlog[3*3], const SO3<S> &X);
    
    template <class S>
    void adjoint(S adj[3*3], const SO3<S> &g);
    
    template <class S>
    void adjoint_multiply(S y[3], const SO3<S> &g, const S x[3]);

    template <class S>
    void adjoint_T_multiply(S y[3], const SO3<S> &g, const S x[3]);

    // Expects |a| == |b| == 1
    // Returns true on success
    template <class S>
    bool compute_rotation_between_unit_vectors(SO3<S> &R, const S a[3], const S b[3]);
}
