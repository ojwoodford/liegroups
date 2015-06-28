#pragma once

#include <cstddef>

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

    template <class G>
    struct QuinticSplineControlPoint
    {
        typedef typename G::Scalar S;

        S t;
        G y;
        S dy[G::DoF];
        S d2y[G::DoF];
    };

    
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

        // Initialize the segment from boundary conditions.
        // After this call,
        //   y(t0) =   y0,   y(t1) =   y1
        //  dy(t0) =  dy0,  dy(t1) =  dy1
        // d2y(t0) = d2y0, d2y(t1) = d2y1
        //
        // Returns true on success.
        //
        bool init(S t0, S t1,
                  const G& y0, const G& y1,
                  const S dy0[G::DoF], const S dy1[G::DoF],
                  const S d2y0[G::DoF], const S d2y1[G::DoF]);

        // Initialize the segment from boundary conditions.
        // Convenience form that calls the above.
        //
        // Returns true on success.
        //
        bool init(const QuinticSplineControlPoint<G>& c0,
                  const QuinticSplineControlPoint<G>& c1)
        {
            return init(c0.t, c1.t, c0.y, c1.y, c0.dy, c1.dy, c0.d2y, c1.d2y);
        }
            
        

        // Compute y(t), dy = y'(t), d2y = y''(t)
        // If dy_dp is not null, it gets the differentials of [y, dy, d2y] by the boundary conditions.
        // dy_dp is indexed first by column (y0, dy0, dy1, y1, dy1, d2y1),
        // then by row (y, dy, d2y).
        //
        // The time t need not be within the boundary interval.
        //
        void eval(G &y,
                  S dy[G::DoF],
                  S d2y[G::DoF],
                  S dy_dp[6][3][G::DoF * G::DoF], // [ y0, dy0, d2y0, y1, dy1, d2y1 ]
                  S t) const;

        S get_t0() const { return t0; }
        S get_t1() const { return t1; }
        
    public:
        G y0;
        S t0, t1;
        
        S pp[5 * N];
        G delta;
        S dlog_delta_ddelta[N*N];
    };


    // This class manages a buffer of spline segments,
    // matching up boundary conditions to form a C2 curve in G.
    //
    // The buffer is owned outside the instance.
    // No memory management is performed internally.
    //
    template <class G>
    class QuinticSpline
    {
    public:
        typedef typename G::Scalar S;
        
        QuinticSpline();

        // Storage lifetime must contain instance lifetime
        // Clears any initialization; init() must be called after this.
        //
        void set_storage(QuinticSplineSegment<G> *segments, size_t max_segments);

        // Returns storage capacity
        size_t get_max_segments() const { return m_max_segments; }
                
        // Returns true on succes.
        //
        // Returns false if storage doesn't have num_controls - 1 slots,
        // or if control_points[] doesn't have monotonically ascending times.
        //
        // After success, get_num_segments() == num_controls - 1.
        //
        bool init(const QuinticSplineControlPoint<G> control_points[],
                  size_t num_controls);

        bool empty() const { return m_num_segments == 0; }

        // Remove all segments.
        // Does not alter storage.
        void clear();
        
        // Add control point to back, creating an additional segment on
        // all but the first invocation on an empty spline.
        //
        // Returns false if storage is insufficient or if control_point.t is not increasing.
        //
        bool push_back(const QuinticSplineControlPoint<G> &control_point);

        // Returns number of initialized segments
        size_t get_num_segments() const { return m_num_segments; }

        // Find the index of the segment containing time t (clamped to front and back)
        size_t get_segment_index_for_time(S t) const;

        // Get the i'th segment
        // No bounds checking here!
        //
        const QuinticSplineSegment<G>& segment(size_t i) const
        {
            return m_segments[i];
        }

        // Get the segment containing time t (clamped to front and back)
        const QuinticSplineSegment<G>& get_segment_for_time(S t) const
        {
            return m_segments[get_segment_index_for_time(t)];
        }

        // Evaluate the spline at time t,
        // returning the index of the segment for t.
        //
        size_t eval(G& y,
                    S dy[G::DoF],
                    S d2y[G::DoF],
                    S dy_dp[6][3][G::DoF * G::DoF], // can be null
                    S t) const
        {
            size_t seg = get_segment_index_for_time(t);
            segment(seg).eval(y, dy, d2y, dy_dp, t);
            return seg;
        }

        // Copy segments from other.
        // Returns false if storage is insufficent (and alters nothing).
        //
        bool copy_from(const QuinticSpline<G> &other);
        
        // Swaps storage.
        // Constant time.
        void swap(QuinticSpline<G> &other);
        
    private:
        // No copying or assignment
        QuinticSpline(const QuinticSpline<G> &other);
        void operator=(const QuinticSpline<G> &other);
        
        QuinticSplineSegment<G> *m_segments;
        size_t m_max_segments;
        size_t m_num_segments;

        QuinticSplineControlPoint<G> m_last_ctrl;
        bool m_has_ctrl;
    };
}
