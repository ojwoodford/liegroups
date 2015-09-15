#pragma once

#include <liegroups/scalar.hpp>
#include <liegroups/matrix.hpp>

#include <liegroups/se2.hpp>
#include <liegroups/so3.hpp>
#include <liegroups/se3.hpp>
#include <iostream>
#include <cstdlib>

#ifndef __PRETTY_FUNCTION__
#define __PRETTY_FUNCTION__ __FUNCTION__
#endif

namespace {
    
    inline double rand_uniform()
    {
        return (double)std::rand() / (double)RAND_MAX;
    }

    template <int N, class S>
    void random_vec(S x[N], S scale = (S)1)
    {
        for (int i=0; i<N; ++i)
            x[i] = scale * (S)(rand_uniform()*2.0 - 1.0);
    }

    template <class S>
    std::ostream &print_vec(std::ostream &out, const S x[], int n)
    {
        const int w = out.precision() + 8;
        for (int i=0; i<n; ++i) {
            out.width(w);
            out << x[i];
        }
        return out;
    }

    template <class S>
    std::ostream &print_mat(std::ostream &out, const S x[], int m, int n)
    {
        const int w = out.precision() + 6;
        for (int i=0; i<m; ++i) {
            for (int j=0; j<n; ++j) {
                out.width(w);
                out << x[i*n + j];
            }
            out << std::endl;
        }
        return out;
    }

    template <class S>
    S dist_sq(const S x[], const S y[], int n)
    {
        S sum = 0;
        for (int i=0; i<n; ++i) {
            S d = x[i] - y[i];
            sum += d*d;
        }
        return sum;
    }

    template <class S>
    S max_abs(const S x[], int n)
    {
        S mv = 0;
        for (int i=0; i<n; ++i) {
            S d = liegroups::abs(x[i]);
            mv = liegroups::max(mv, d);
        }
        return mv;
    }
    
    template <class S>
    S max_abs_diff(const S x[], const S y[], int n)
    {
        S mv = 0;
        for (int i=0; i<n; ++i) {
            S d = liegroups::abs(x[i] - y[i]);
            mv = liegroups::max(mv, d);
        }
        return mv;
    }

    static int check_count = 0;

    template <class G, class S>
    void CHECK_ERROR(S err, S max_err, const char *name)
    {
        if (liegroups::abs(err) > max_err) {
            std::cerr << __PRETTY_FUNCTION__ << std::endl;
            std::cerr << "!!!!! " << check_count << " Failed "
                      << name
                      << ": err = " << err << ", max is " << max_err
                      << std::endl;
            std::exit(1);
        }
        ++check_count;
    }
    
    template <class G, typename S>
    void CHECK_ERROR(const S est[], const S ref[], int n, S max_err, const char *name)
    {
        S err = max_abs_diff(est, ref, n);
        if (err > max_err) {
            std::cerr << __PRETTY_FUNCTION__ << std::endl;
            std::cerr << "!!!!! " << check_count << " Failed "
                      << name                
                      << std::endl << "est = ";
            print_vec(std::cerr, est, n) << std::endl;
            std::cerr << "ref = ";
            print_vec(std::cerr, ref, n) << std::endl;
            std::cerr << "err is " << err << ", max is " << max_err << std::endl;
            std::exit(1);
        }
        ++check_count;
    }

    template <class G, typename S>
    void CHECK_ERROR(const S est[], const S ref[], int m, int n, S max_err, const char *name)
    {
        S err = max_abs_diff(est, ref, m*n);
        if (err > max_err) {
            std::cerr << __PRETTY_FUNCTION__ << std::endl;
            std::cerr << "!!!!! " << check_count << " Failed "
                      << name                
                      << std::endl << "est = " << std::endl;
            print_mat(std::cerr, est, m, n) << std::endl;
            std::cerr << "ref = " << std::endl;
            print_mat(std::cerr, ref, m, n) << std::endl;
            std::cerr << "err is " << err << ", max is " << max_err << std::endl;
            std::exit(1);
        }
        ++check_count;
    }
    
    template <class G>
    struct VecGen
    {
        typedef typename G::Scalar S;
        static void gen_vec(S x[]) {
            random_vec<G::DoF>(x);
        }
    };

    template <class S>
    struct VecGen<liegroups::SE2<S> >
    {
        static void gen_vec(S x[]) {
            const S pi = (S)3.1415926535897932384626433;
            random_vec<3>(x, pi);
        }
    };

    template <class S>
    struct VecGen<liegroups::SO3<S> >
    {
        static void gen_vec(S x[]) {
            const S pi = (S)3.1415926535897932384626433;
            random_vec<3>(x, pi);
            S xx = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
            S f = (S)0.999999;
            S max_theta = f * pi;
            if (xx >= max_theta * max_theta) {
                S f = max_theta / liegroups::sqrt(xx);
                x[0] *= f;
                x[1] *= f;
                x[2] *= f;
            }
        }
    };

    template <class S>
    struct VecGen<liegroups::SE3<S> >
    {
        static void gen_vec(S x[]) {
            random_vec<3>(&x[0]);
            VecGen<liegroups::SO3<S> >::gen_vec(&x[3]);
        }
    };


    template <class X, typename Out> struct ChangeScalar;

    template <template <typename S> class G, typename S, typename Out>
    struct ChangeScalar<G<S>, Out>
    {
        typedef G<Out> type;
    };

    
}
