#pragma once

namespace liegroups {

    template <typename S>
    struct ExpCoefs
    {
        S cos_theta;
        S A, B, C;
        void compute(S theta_sq);        
        ExpCoefs() {}
        ExpCoefs(S theta_sq) { compute(theta_sq); }
    };

    template <typename S>
    struct DiffExpCoefs : public ExpCoefs<S>
    {
        S A1, B1, C1;
        void compute(S theta_sq);        
        DiffExpCoefs() {}
        DiffExpCoefs(S theta_sq) { compute(theta_sq); }
    };

    template <typename S>
    struct Diff2ExpCoefs : public DiffExpCoefs<S>
    {
        S A2, B2, C2;
        void compute(S theta_sq);        
        Diff2ExpCoefs() {}
        Diff2ExpCoefs(S theta_sq) { compute(theta_sq); }
    };    
}
