#pragma once

namespace liegroups {

    // Compute the group product of geodesics, and its differentials:
    //
    //    y <--- exp(a_{n-1}) * ... * exp(a_0)
    //    dy <--- dy/dt
    //    d2y <--- d/dt(dy/dt)
    //
    // where
    //
    //    a_i(t) = b_i(t) * p_i
    //    b_i(t) = bb[i*3]
    //    b_i'   = bb[i*3+1]
    //    b_i''  = bb[i*3+2]
    //    p_i    = pp[i*DoF ... i*DoF + (DoF-1)]
    //
    // (DoF = number of degrees of freedom in group G)
    //
    // If dy_dp is not null, it gets the differentials of (y,dy,d2y) by pp:
    //    dy_dp[i][0] <-- dy/dpi
    //    dy_dp[i][1] <-- ddy/dpi
    //    dy_dp[i][2] <-- dd2y/dpi   
    //
    // Differentials in the group by parameters 'q' are defined as:
    //    dy/dq = diff(log(y(q + eps) * inv(y(q))), eps) at eps = 0
    //
    // Outputs:
    // *  dy and d2y must have DoF slots
    // *  If dy_dp is not null, it must have size [n][3][DoF*DoF].
    //
    // Inputs:
    // *  n indicates how many elements are in the product chain
    // *  pp must have n*DoF slots (n vectors in the algebra)
    // *  bb must have n*3 slots (every triplet is b_i, b_i' b_i'')
    //
    // Scratch:
    // *  scratch0 must have at least n slots
    // *  scratch1 must have at least n*DoF slots
    //
    template <class G, typename S>
    void eval_product_chain(G &y,
                            S dy[],
                            S d2y[],
                            S dy_dp[][3][G::DoF * G::DoF],
                            int n,
                            const S pp[],
                            const S bb[],
                            G scratch0[],
                            S scratch1[]);


    // A time-parametrized curve in the group G specified by boundary conditions:
    // *  Time at boundaries
    // *  Value in the group at boundaries
    // *  First time derivative at boundaries
    // *  Second time derivative at boundaries
    //
    // Invalid until init() is called.
    //
    template <class G>
    class QuinticSplineSegment
    {
    public:
        
        typedef typename G::Scalar S;
        static const int N = G::DoF;

        QuinticSplineSegment();
        
        bool init(S t0, S t1,
                  const G& y0, const G& y1,
                  const S dy0[G::DoF], const S dy1[G::DoF],
                  const S d2y0[G::DoF], const S d2y1[G::DoF]);
        
        void eval(G &y,
                  S dy[G::DoF],
                  S d2y[G::DoF],
                  S dy_dp[6][3][G::DoF * G::DoF], // [ y0, dy0, d2y0, y1, dy1, d2y1 ]
                  S t) const;

    private:
        G y0;
        S t0, t1;
        
        S pp[5 * N];
        G delta;
        S dlog_delta_ddelta[N*N];
    };
}
