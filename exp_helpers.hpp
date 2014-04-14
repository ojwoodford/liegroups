#pragma once

namespace liegroups {

    // Compute m = A*I + B*wx + C*ww'
    template <typename S>
    void compute_exp_matrix3(S m0[3], S m1[3], S m2[3],
                             S a, S b, S c, const S w[3])
    {
        const S Cw0 = c * w[0];
        const S Cw1 = c * w[1];
        const S Cw2 = c * w[2];
        m0[0] = a + Cw0 * w[0];
        m1[1] = a + Cw1 * w[1];
        m2[2] = a + Cw2 * w[2];
    
        const S Cw01 = Cw0 * w[1];
        const S Cw02 = Cw0 * w[2];
        const S Cw12 = Cw1 * w[2];
        const S Bw0 = b * w[0];
        const S Bw1 = b * w[1];
        const S Bw2 = b * w[2];
    
        m0[1] = Cw01 - Bw2;
        m0[2] = Cw02 + Bw1;
        m1[0] = Cw01 + Bw2;
        m1[2] = Cw12 - Bw0;
        m2[0] = Cw02 - Bw1;
        m2[1] = Cw12 + Bw0;
    }

    // Convenience wrapper for contiguous m
    template <typename S>
    void compute_exp_matrix3(S m[3*3],
                             S a, S b, S c, const S w[3])
    {
        compute_exp_matrix3(&m[0], &m[3], &m[6], a, b, c, w);
    }
    
    
}
