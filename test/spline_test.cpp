#include "test_helpers.hpp"
#include <liegroups/spline.hpp>
#include <liegroups/se2_io.hpp>
#include <liegroups/so3_io.hpp>
#include <liegroups/se3_io.hpp>
#include <cassert>

using namespace std;
using namespace liegroups;

template <typename S>
void copy(SE2<double> &out, const SE2<S> &in)
{
    for (int i=0; i<2; ++i) {
        out.r[i] = in.r[i];
        out.t[i] = in.t[i];
    }
    rectify(out);
}

template <typename S>
void copy(SO3<double> &out, const SO3<S> &in)
{
    for (int i=0; i<9; ++i) {
        out.R[i] = in.R[i];
    }
    rectify(out);
}

template <typename S>
void copy(SE3<double> &out, const SE3<S> &in)
{
    copy(out.R, in.R);
    for (int i=0; i<3; ++i) {
        out.t[i] = in.t[i];
    }
}

template <class G>
void diff(typename G::Scalar d[G::DoF], const G &hi, const G &lo, typename G::Scalar eps)
{
    typedef typename ChangeScalar<G,double>::type Gd;
    Gd hi_d, lo_d;
    copy(hi_d, hi);
    copy(lo_d, lo);

    Gd delta;
    multiply_a_binv(delta, hi_d, lo_d);
    double log_delta[G::DoF];
    log(log_delta, delta);

    for (int i=0; i<G::DoF; ++i) {
        d[i] = (typename G::Scalar)(log_delta[i] / (2*eps));
    }    
}

const char *const dy_dp_name[6][3] = {
    { "dy_dy0", "ddy_dy0", "dd2y_dy0" },
    { "dy_ddy0", "ddy_ddy0", "dd2y_ddy0" },
    { "dy_dd2y0", "ddy_dd2y0", "dd2y_dd2y0" },
    { "dy_dy1", "ddy_dy1", "dd2y_dy1" },
    { "dy_ddy1", "ddy_ddy1", "dd2y_ddy1" },
    { "dy_dd2y1", "ddy_dd2y1", "dd2y_dd2y1" }
};

template <class G>
void test_spline_segment()
{
    typedef typename G::Scalar S;

    const int N = G::DoF;
    const S max_err = (S)10 * Constants<S>::epsilon() * (S)N;
    const S max_big_err = liegroups::sqrt(Constants<S>::sqrt_epsilon());
    
    const S eps = liegroups::sqrt(liegroups::Constants<S>::sqrt_epsilon());
    
    // Boundary conditions
    const S t[2] = {-1, 2};
    G y[2];
    S dy[2][N];
    S d2y[2][N];

    for (int i=0; i<2; ++i) {
        S log_y[N];
        VecGen<G>::gen_vec(log_y);
        if (i == 1) {
            for (int j=0; j<N; ++j) {
                log_y[j] *= (S)0.75;
            }
        }
        exp(y[i], log_y);

        random_vec<N>(dy[i], (S)1);
        random_vec<N>(d2y[i], (S)2);
    }

    y[1] = y[1] * y[0];

    
    QuinticSplineSegment<G> seg;
    bool init_success = seg.init(t[0], t[1], y[0], y[1], dy[0], dy[1], d2y[0], d2y[1]);
    assert(init_success);

    G yt;
    S dyt[N];
    S d2yt[N];
    
    // Check boundaries
    for (int i=0; i<2; ++i) {
        seg.eval(yt, dyt, d2yt, 0, t[i]);

        G delta;
        multiply_a_binv(delta, yt, y[i]);
            
        S log_delta[N];
        log(log_delta, delta);

        CHECK_ERROR<G>(max_abs(log_delta, N), max_big_err, "eval_boundary_y");
        CHECK_ERROR<G>(max_abs_diff(dyt, dy[i], N), max_big_err, "eval_boundary_dy");
        CHECK_ERROR<G>(max_abs_diff(d2yt, d2y[i], N), max_big_err, "eval_boundary_d2y");
    }

    // Check differentials
    S dy_dp[6][3][N*N];
    const S dtt = (t[1] - t[0]) * (S)0.1;
    for (S tt=t[0] - dtt; tt<=t[1] + dtt; tt += dtt)
    {
        seg.eval(yt, dyt, d2yt, dy_dp, tt);

        G y_lo, y_hi;
        S dy_lo[N], dy_hi[N];
        S d2y_lo[N], d2y_hi[N];
        S dummy[N];
        
        {
            const S dt = eps;
            seg.eval(y_lo, dy_lo, dummy, 0, tt - dt);
            seg.eval(y_hi, dy_hi, dummy, 0, tt + dt);

            S dy_num[N];            
            diff(dy_num, y_hi, y_lo, dt);
            CHECK_ERROR<G>(dyt, dy_num, N, max_big_err, "dy");
            
            S d2y_num[N];
            for (int i=0; i<N; ++i) {
                d2y_num[i] = (S)(((double)dy_hi[i] - (double)dy_lo[i]) / (2 * dt));
            }
            CHECK_ERROR<G>(d2yt, d2y_num, N, max_big_err, "d2y");
        }

        QuinticSplineSegment<G> segmod;
        
        for (int i=0; i<2; ++i) {
            const G oldy = y[i];

            S dy_num[N*N];
            S ddy_num[N*N];
            S dd2y_num[N*N];
            S mod[N] = {0};
            
            for (int j=0; j<3; ++j) {
                for (int k=0; k<N; ++k) {
                    G delta;
                    S old;
                    if (j == 0) {
                        mod[k] = eps;
                        exp(delta, mod);
                        multiply(y[i], delta, oldy);
                    } else if (j == 1) {
                        old = dy[i][k];
                        dy[i][k] = old + eps;
                    } else {
                        old = d2y[i][k];
                        d2y[i][k] = old + eps;
                    }
                    bool init_success = segmod.init(t[0], t[1], y[0], y[1], dy[0], dy[1], d2y[0], d2y[1]);
                    assert(init_success);
                    segmod.eval(y_hi, dy_hi, d2y_hi, 0, tt);

                    if (j == 0) {
                        multiply(y[i], inverse(delta), oldy);
                    } else if (j == 1) {
                        dy[i][k] = old - eps;
                    } else {
                        d2y[i][k] = old - eps;
                    }
                    
                    init_success = segmod.init(t[0], t[1], y[0], y[1], dy[0], dy[1], d2y[0], d2y[1]);
                    assert(init_success);
                    segmod.eval(y_lo, dy_lo, d2y_lo, 0, tt);

                    if (j == 0) {
                        mod[k] = 0;
                    } else if (j == 1) {
                        dy[i][k] = old;
                    } else {
                        d2y[i][k] = old;
                    }

                    S dy_num_k[N];
                    diff(dy_num_k, y_hi, y_lo, eps);

                    for (int l=0; l<N; ++l) {
                        dy_num[k + l*N] = dy_num_k[l];
                        ddy_num[k + l*N] = (S)(((double)dy_hi[l] - (double)dy_lo[l]) / (2 * eps));
                        dd2y_num[k + l*N] = (S)(((double)d2y_hi[l] - (double)d2y_lo[l]) / (2 * eps));
                    }
                }
                
                y[i] = oldy;
                                
                CHECK_ERROR<G>(dy_num, dy_dp[i*3 + j][0], N, N, max_big_err, dy_dp_name[i*3+j][0]);
                CHECK_ERROR<G>(ddy_num, dy_dp[i*3 + j][1], N, N, max_big_err, dy_dp_name[i*3+j][1]);
                CHECK_ERROR<G>(dd2y_num, dy_dp[i*3 + j][2], N, N, max_big_err, dy_dp_name[i*3+j][2]);
            }
        }
    }        
}


int main()
{
    srand(47);
    const int outer_passes = 100;
    for (int op = 0; op < outer_passes; ++op) {
        const int passes = 100;
        
        for (int pass=0; pass<passes; ++pass) test_spline_segment<SE2<float> >();
        for (int pass=0; pass<passes; ++pass) test_spline_segment<SE2<double> >();

        for (int pass=0; pass<passes; ++pass) test_spline_segment<SO3<float> >();
        for (int pass=0; pass<passes; ++pass) test_spline_segment<SO3<double> >();

        for (int pass=0; pass<passes; ++pass) test_spline_segment<SE3<float> >();
        for (int pass=0; pass<passes; ++pass) test_spline_segment<SE3<double> >();
    }
    
    return 0;
}
