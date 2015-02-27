#include <liegroups/se2.hpp>
#include <liegroups/so3.hpp>
#include <liegroups/se3.hpp>

#include <liegroups/matrix_impl.hpp>

#include "spline.hpp"

namespace {
    template <int N, typename S>
    void scale(S out[N], const S in[N], S scale)
    {
        for (int i=0; i<N; ++i)
            out[i] = in[i] * scale;
    }

    template <int N, typename S>
    void zero(S out[]) {
        for (int i=0; i<N; ++i)
            out[i] = 0;
    }

    template <int N, typename S>
    void add(S out[N], const S a[N], const S b[N])
    {
        for (int i=0; i<N; ++i)
            out[i] = a[i] + b[i];
    }

    template <int N, typename S>
    void add(S out[N], const S b[N])
    {
        add<N,S>(out, out, b);
    }

    template <int N, typename S>
    void add_scaled(S out[N], const S a[N], const S b[N], S bscale)
    {
        for (int i=0; i<N; ++i)
            out[i] = a[i] + b[i] * bscale;
    }

    template <int N, typename S>
    void set_identity(S m[], S d=1)
    {
        zero<N*N>(m);
        for (int i=0; i<N; ++i)
            m[i*(N+1)] = d;
    }

    template <int N, typename S>
    void add_identity(S m[], S d=1)
    {
        for (int i=0; i<N; ++i)
            m[i*(N+1)] += d;
    }

    template <int N, typename S, class Op>
    void mat_mult_square(S ab[], const S a[], const S b[], const Op &op)
    {
        for (int i=0; i<N; ++i) {
            const S* const ai = &a[i*N];
            for (int j=0; j<N; ++j) {
                S sum = 0;
                for (int k=0; k<N; ++k)
                    sum += ai[k] * b[j + k*N];
                op(ab[i*N + j], sum);
            }
        }
    }

    namespace ops
    {
        struct Add {
            template <class L, class R>
            void operator()(L &lhs, const R &rhs) const {
                lhs += rhs;
            }
        };

        struct Subtract {
            template <class L, class R>
            void operator()(L &lhs, const R &rhs) const {
                lhs -= rhs;
            }
        };
    }
}

template <class G, typename S>
void liegroups::eval_product_chain(G &y,
                                   S dy[],
                                   S d2y[],
                                   S dy_dp[][3][G::DoF * G::DoF],
                                   int n,
                                   const S pp[],
                                   const S bb[],
                                   G scratch0[],
                                   S scratch1[])
{
    const int N = G::DoF;
    
    if (n == 0) {
        y = G::identity;
        zero<N>(dy);
        zero<N>(d2y);
        
        if (dy_dp) {
            for (int i=0; i<3; ++i) {
                for (int j=0; j<n; ++j) {
                    zero<N*N>(dy_dp[j][i]);
                }
            }
        }
        
        return;
    }

    G *const As = scratch0;
    S *const A1s = scratch1;

    const S *pi = pp;
    const S *bi = bb;
    for (int i=0; i<n; ++i, pi += N, bi += 3) {
        const S b = bi[0];
        const S db = bi[1];
        const S d2b = bi[2];

        S x[N];
        scale<N>(x, pi, b);

        G A;
        S D[N*N];
        
        if (dy_dp) {
            exp_diff(A, D, x);
        } else {
            exp(A, x);
        }

        if (i == 0) {
            y = A;
            scale<N>(dy, pi, db);
            scale<N>(d2y, pi, d2b);

            if (dy_dp) {
                scale<N*N>(dy_dp[0][0], D, b);
                set_identity<N>(dy_dp[0][1], db);
                set_identity<N>(dy_dp[0][2], d2b);
            }            
        } else {
            S A1[N];
            scale<N>(A1, pi, db);

            S AdA_dy[N];
            S AdA_d2y[N];
            adjoint_multiply(AdA_dy, A, dy);
            adjoint_multiply(AdA_d2y, A, d2y);

            S adA1_AdA_dy[N];
            G::ad_multiply(adA1_AdA_dy, A1, AdA_dy);
            
            y = A * y;
            
            add<N>(dy, A1, AdA_dy);
            
            scale<N>(d2y, pi, d2b);
            add<N>(d2y, AdA_d2y);
            add<N>(d2y, adA1_AdA_dy);

            if (dy_dp) {
                As[i] = A;
                copy<N>(&A1s[i*N], A1);
                
                S *const dy_dpi = dy_dp[i][0];
                S *const ddy_dpi = dy_dp[i][1];
                S *const dd2y_dpi = dy_dp[i][2];
                
                scale<N*N>(dy_dpi, D, b);

                S tmp[N*N];
                
                S neg_AdA_dy[N];
                scale<N>(neg_AdA_dy, AdA_dy, (S)-1);
                G::ad(tmp, neg_AdA_dy);
                mat_mult_square<N>(ddy_dpi, tmp, dy_dpi);

                G::ad(tmp, A1);                
                mat_mult_square<N>(dd2y_dpi, tmp, ddy_dpi);
                
                add_identity<N>(ddy_dpi, db);
                
                G::ad(tmp, AdA_d2y);
                ::mat_mult_square<N>(dd2y_dpi, tmp, dy_dpi, ops::Subtract());
                
                G::ad(tmp, AdA_dy);
                add_scaled<N*N>(dd2y_dpi, dd2y_dpi, tmp, -db);
                add_identity<N>(dd2y_dpi, d2b);
            }
        }        
    }

    if (dy_dp) {
        G Li = As[n-1];
        S *const li = &A1s[(n-1)*N];
        for (int i=n-2; i>=0; --i) {            
            S *const dy_dpi = dy_dp[i][0];
            S *const ddy_dpi = dy_dp[i][1];
            S *const dd2y_dpi = dy_dp[i][2];

            S AdLi[N*N];
            adjoint(AdLi, Li);

            S tmp[N*N];
            
            mat_mult_square<N>(tmp, AdLi, dy_dpi);
            copy<N*N>(dy_dpi, tmp);

            mat_mult_square<N>(tmp, AdLi, ddy_dpi);
            copy<N*N>(ddy_dpi, tmp);

            mat_mult_square<N>(tmp, AdLi, dd2y_dpi);
            copy<N*N>(dd2y_dpi, tmp);

            G::ad(tmp, li);
            ::mat_mult_square<N>(dd2y_dpi, tmp, ddy_dpi, ops::Add());

            if (i > 0) {
                adjoint_multiply(tmp, Li, &A1s[i*N]);
                add<N>(li, tmp);
                
                Li = Li * As[i];
            }
        }
    }
}

template <class G>
liegroups::QuinticSplineSegment<G>::QuinticSplineSegment()
{
    t0 = 0;
    t1 = 1;

    y0 = G::identity;
    delta = G::identity;
    set_identity<N>(dlog_delta_ddelta);

    zero<5*N>(pp);
}

template <class G>
bool liegroups::QuinticSplineSegment<G>::init(S t0_, S t1_,
                                              const G& y0_, const G& y1,
                                              const S dy0[G::DoF], const S dy1[G::DoF],
                                              const S d2y0[G::DoF], const S d2y1[G::DoF])
{
    t0 = t0_;
    t1 = t1_;

    const S dt = t1 - t0;
        
    y0 = y0_;
    multiply_a_binv(delta, y1, y0);
        
    scale<N>(&pp[0], dy0, dt);
    scale<N>(&pp[N], d2y0, dt*dt);

    scale<N>(&pp[3*N], dy1, dt);
    scale<N>(&pp[4*N], d2y1, dt*dt);
        
    S *const log_delta = &pp[2*N];
    log_diff(log_delta, dlog_delta_ddelta, delta);
        
    return true;
}

template <class G>
void liegroups::QuinticSplineSegment<G>::eval(G &y,
                                              S dy[G::DoF],
                                              S d2y[G::DoF],
                                              S dy_dp[6][3][G::DoF * G::DoF],
                                              S t) const
{
    const S dt = t1 - t0;
    const S inv_dt = 1 / dt;
    const S s = (t - t0) * inv_dt;
     
    // All constant coefficients are zero;
    // c[5*i] is the i'th linear coefficient.
    const S c[5*5] = {
        1, 0, -6, 8, -3,
        0, (S)0.5, (S)-1.5, (S)1.5, (S)-0.5,
        0, 0, 10, -15, 6,
        0, 0, -4, 7, -3,
        0, 0, (S)0.5, -1, (S)0.5
    };

    const S ds_dt = inv_dt;
    S bb[5*3];
    for (int i=0; i<5; ++i) {
        const S* const ci = &c[5*i];
        bb[3*i] = s*(ci[0] + s*(ci[1] + s*(ci[2] + s*(ci[3] + s*ci[4]))));
        bb[3*i+1] = ds_dt * (ci[0] + s*(2*ci[1] + s*(3*ci[2] + s*(4*ci[3] + s*5*ci[4]))));
        bb[3*i+2] = ds_dt * ds_dt * (2*ci[1] + s*(6*ci[2] + s*(12*ci[3] + s*20*ci[4])));
    }

    G scratch0[5];
    S scratch1[5*N];
    eval_product_chain(y, dy, d2y,
                       (dy_dp ? &dy_dp[1] : 0),
                       5, pp, bb,
                       scratch0, scratch1);
        
    if (dy_dp) {
        S neg_Ad_delta[N*N];
        adjoint(neg_Ad_delta, delta);
        scale<N*N>(neg_Ad_delta, neg_Ad_delta, (S)-1);

        adjoint(dy_dp[0][0], y);
        for (int i=0; i<3; ++i) {
            S ab[N*N];
            mat_mult_square<N>(ab, dy_dp[3][i], dlog_delta_ddelta);
            copy<N*N>(dy_dp[3][i], ab);

            if (i==0) {
                ::mat_mult_square<N>(dy_dp[0][i], ab, neg_Ad_delta, ops::Add());
            } else {
                mat_mult_square<N>(dy_dp[0][i], ab, neg_Ad_delta);
            }

            scale<N*N>(dy_dp[1][i], dy_dp[1][i], dt);
            scale<N*N>(dy_dp[2][i], dy_dp[2][i], dt*dt);
            scale<N*N>(dy_dp[4][i], dy_dp[4][i], dt);
            scale<N*N>(dy_dp[5][i], dy_dp[5][i], dt*dt);
        }
    }
        
    y = y * y0;
}

template class liegroups::QuinticSplineSegment<liegroups::SE2<float> >;
template class liegroups::QuinticSplineSegment<liegroups::SE2<double> >;

template class liegroups::QuinticSplineSegment<liegroups::SO3<float> >;
template class liegroups::QuinticSplineSegment<liegroups::SO3<double> >;

template class liegroups::QuinticSplineSegment<liegroups::SE3<float> >;
template class liegroups::QuinticSplineSegment<liegroups::SE3<double> >;
