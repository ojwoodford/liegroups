#include "test_helpers.hpp"

#include <liegroups/se2_io.hpp>
#include <liegroups/so3_io.hpp>
#include <liegroups/se3_io.hpp>

#include <liegroups/sim2.hpp>
#include <liegroups/sim2_io.hpp>
#include <liegroups/aff2.hpp>
#include <liegroups/aff2_io.hpp>
#include <liegroups/sl3.hpp>
#include <liegroups/sl3_io.hpp>
#include <liegroups/sim3.hpp>
#include <liegroups/sim3_io.hpp>

#include <cassert>

using namespace std;
using namespace liegroups;

template <class G>
void test_group()
{
    typedef typename G::Scalar S;
    const int N = G::DoF;
    const int D = G::Dim;
    const S max_err = (S)10 * Constants<S>::epsilon() * (S)N;
    const S max_big_err = Constants<S>::sqrt_epsilon();

    G g = G::identity;

    {
        S log_g[N];
        VecGen<G>::gen_vec(log_g);

        exp(g, log_g);
        
        S log_exp[N] = { (S)-9999999 };
        log(log_exp, g);
        CHECK_ERROR<G>(log_exp, log_g, N, max_big_err, "log");
    }

    {
        S log_ginvg[N];
        S log_identity[N] = {0};
        log(log_ginvg, inverse(g) * g);
        CHECK_ERROR<G>(log_ginvg, log_identity, N, max_err*10, "log(inv(a)*a)");

        G gginv;
        multiply_a_binv(gginv, g, g);
        log(log_ginvg, gginv);
        CHECK_ERROR<G>(log_ginvg, log_identity, N, max_err*10, "log(a*inv(a))");
    }
    
    {
        S adj[N*N];
        adjoint(adj, g);

        S adj_v[N], adj_v_ref[N], adj_T_v[N], adj_T_v_ref[N], v[N];
        random_vec<N>(v);

        for (int i=0; i<N; ++i) {
            S ai = 0;
            S aTi = 0;
            for (int j=0; j<N; ++j) {
                ai += adj[i*N + j] * v[j];
                aTi += adj[j*N + i] * v[j];
            }
            adj_v_ref[i] = ai;
            adj_T_v_ref[i] = aTi;
        }
        
        adjoint_multiply(adj_v, g, v);
        CHECK_ERROR<G>(adj_v, adj_v_ref, N, max_err*10, "adj");

        adjoint_T_multiply(adj_T_v, g, v);
        CHECK_ERROR<G>(adj_T_v, adj_T_v_ref, N, max_err*10, "adjT");
        
        G exp_v;
        exp(exp_v, v);
        
        G conj = g * exp_v * inverse(g);

        G exp_adj_v;        
        exp(exp_adj_v, adj_v);

        G delta;
        multiply_a_binv(delta, conj, exp_adj_v);

        S log_delta[N] = {(S)-1};
        S log_identity[N] = {(S)0};
        log(log_delta, delta);
        
        CHECK_ERROR<G>(log_delta, log_identity, N, max_big_err*10, "adj mult");
    }
    
    {     
        S x[D];
        random_vec<D>(x);

        S y[D], z[D];
        transform_point(y, g, x);
        transform_point_by_inverse(z, g, y);        
        CHECK_ERROR<G>(z, x, D, max_err*10, "transform");
        transform_point(z, inverse(g), y);
        CHECK_ERROR<G>(z, x, D, max_err*10, "transform inv");
    }
}

template <typename S>
void test_rotv2v()
{
    const S eps = Constants<S>::epsilon();
    const S max_err = (S)100 * eps;
    
    S a[3], b[3];
    S aa, bb;
    while (true) {
        random_vec<3>(a);
        random_vec<3>(b);
        
        aa = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
        bb = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
        if (aa > eps && bb > eps)
            break;
    }

    S fa = (S)1 / liegroups::sqrt(aa);
    S fb = (S)1 / liegroups::sqrt(bb);
    for (int i=0; i<3; ++i) {
        a[i] *= fa;
        b[i] *= fb;
    }

    SO3<S> R;
    compute_rotation_between_unit_vectors(R, a, b);

    S c[3];
    transform_point(c, R, a);
    CHECK_ERROR<SO3<S> >(c, b, 3, max_err, "rotv2v");
}

template <class G>
void test_exp_diff()
{
    typedef typename G::Scalar S;
    typedef typename ChangeScalar<G,double>::type Gd;
    
    const int N = G::DoF;
    const S max_err = (S)10 * Constants<S>::epsilon() * (S)N;
    const S max_big_err = Constants<S>::sqrt_epsilon();

    S err;
    
    G g;
    S log_g[N];
    VecGen<G>::gen_vec(log_g);

    S dexp[N*N];
    exp_diff(g, dexp, log_g);

    S num_dexp[N*N];
    {
        // Compute numerical Jacobian w/ double precision
        double log_g_mod[N];
        for (int i=0; i<N; ++i) {
            log_g_mod[i] = log_g[i];
        }
        for (int i=0; i<N; ++i) {
            const double eps = (S)1e-4;        
            const double old = log_g_mod[i];
            Gd hi, lo;
            log_g_mod[i] = old + eps;
            exp(hi, log_g_mod);
            log_g_mod[i] = old - eps;
            exp(lo, log_g_mod);
            log_g_mod[i] = old;

            Gd delta;
            multiply_a_binv(delta, hi, lo);

            double log_d[N];
            log(log_d, delta);

            for (int j=0; j<N; ++j) {
                num_dexp[i + j*N] = static_cast<S>(log_d[j] / (2 * eps));
            }
        }
    }
    CHECK_ERROR<G>(dexp, num_dexp, N, N, max_big_err, "exp_diff_dexp");

    S dlog[N*N];
    S log_exp[N] = { (S)-9999999 };
    log_diff(log_exp, dlog, g);

    CHECK_ERROR<G>(log_exp, log_g, N, max_big_err, "log_diff_log");

    S mat[N*N];
    mat_mult_square<N>(mat, dlog, dexp);

    S id[N*N] = {0};

    for (int i=0; i<N; ++i) {
        id[i*(N+1)] = 1;
    }

    CHECK_ERROR<G>(mat, id, N, N, max_big_err, "log_diff_diffprod");
}


int main()
{
    srand(47);
    const int passes = 100000;
    for (int pass=0; pass<passes; ++pass)
    {
        test_rotv2v<float>();
        test_rotv2v<double>();
    }
    check_count = 0;


    for (int pass=0; pass<passes; ++pass) test_exp_diff<SE2<float> >();
    for (int pass=0; pass<passes; ++pass) test_exp_diff<SE2<double> >();

    for (int pass=0; pass<passes; ++pass) test_exp_diff<SO3<float> >();
    for (int pass=0; pass<passes; ++pass) test_exp_diff<SO3<double> >();

    for (int pass=0; pass<passes; ++pass) test_exp_diff<SE3<float> >();
    for (int pass=0; pass<passes; ++pass) test_exp_diff<SE3<double> >();

    
    for (int pass=0; pass<passes; ++pass) test_group<SE2<float> >();
    for (int pass=0; pass<passes; ++pass) test_group<SE2<double> >();

    for (int pass=0; pass<passes; ++pass) test_group<Sim2<float> >();
    for (int pass=0; pass<passes; ++pass) test_group<Sim2<double> >();

    for (int pass=0; pass<passes; ++pass) test_group<Aff2<float> >();
    for (int pass=0; pass<passes; ++pass) test_group<Aff2<double> >();
    
    for (int pass=0; pass<passes; ++pass) test_group<SO3<float> >();
    for (int pass=0; pass<passes; ++pass) test_group<SO3<double> >();

    for (int pass=0; pass<passes; ++pass) test_group<SE3<float> >();
    for (int pass=0; pass<passes; ++pass) test_group<SE3<double> >();

    for (int pass=0; pass<passes; ++pass) test_group<Sim3<float> >();
    for (int pass=0; pass<passes; ++pass) test_group<Sim3<double> >();
    
    //for (int pass=0; pass<passes; ++pass) test_group<SL3<float> >();
    for (int pass=0; pass<passes; ++pass) test_group<SL3<double> >();

    
    return 0;
}
