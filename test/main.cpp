#include <liegroups/scalar.hpp>

#include <liegroups/se2.hpp>
#include <liegroups/se2_io.hpp>
#include <liegroups/sim2.hpp>
#include <liegroups/sim2_io.hpp>
#include <liegroups/aff2.hpp>
#include <liegroups/aff2_io.hpp>
#include <liegroups/sl3.hpp>
#include <liegroups/sl3_io.hpp>
#include <liegroups/so3.hpp>
#include <liegroups/so3_io.hpp>
#include <liegroups/se3.hpp>
#include <liegroups/se3_io.hpp>
#include <iostream>
#include <cstdlib>
#include <cassert>

using namespace std;
using namespace liegroups;

double rand_uniform()
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
ostream &print_vec(ostream &out, const S x[], int n)
{
    const int w = out.precision() + 8;
    for (int i=0; i<n; ++i) {
        out.width(w);
        out << x[i];
    }
    return out;
}

template <class S>
ostream &print_mat(ostream &out, const S x[], int m, int n)
{
    const int w = out.precision() + 6;
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            out.width(w);
            out << x[i*n + j];
        }
        out << endl;
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
S max_abs_diff(const S x[], const S y[], int n)
{
    S mv = 0;
    for (int i=0; i<n; ++i) {
        S d = liegroups::abs(x[i] - y[i]);
        mv = liegroups::max(mv, d);
    }
    return mv;
}

int check_count = 0;

template <class G, class S>
void CHECK_ERROR(S err, S max_err, const char *name)
{
    if (liegroups::abs(err) > max_err) {
        cerr << __PRETTY_FUNCTION__ << endl;
        cerr << "!!!!! " << check_count << " Failed "
             << name
             << ": err = " << err << ", max is " << max_err << endl;
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
struct VecGen<SE2<S> >
{
    static void gen_vec(S x[]) {
        const S pi = (S)3.1415926535897932384626433;
        random_vec<3>(x, pi);
    }
};

template <class S>
struct VecGen<SO3<S> >
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
struct VecGen<SE3<S> >
{
    static void gen_vec(S x[]) {
        random_vec<3>(&x[0]);
        VecGen<SO3<S> >::gen_vec(&x[3]);
    }
};

template <class G>
void test_group()
{
    typedef typename G::Scalar S;
    const int N = G::DoF;
    const int D = G::Dim;
    const S max_err = (S)10 * Constants<S>::epsilon() * (S)N;
    const S max_big_err = Constants<S>::sqrt_epsilon();
    //const S max_err_sq = max_err*max_err;

    G g = G::identity;

    {
        S log_g[N];
        VecGen<G>::gen_vec(log_g);

        exp(g, log_g);
        
        S log_exp[N] = { (S)-9999999 };
        log(log_exp, g);
        S err = max_abs_diff(log_exp, log_g, N);
        if (err >= max_big_err) {
            cerr.precision(19);
            cerr << g << endl;
            // for (int i=0; i<N; ++i)
            //     log_exp[i] -= log_g[i];
            
            print_vec(cerr, log_g, N) << endl;
            print_vec(cerr, log_exp, N) << endl;
        }
        CHECK_ERROR<G>(err, max_big_err, "log");
    }

    {
        S log_ginvg[N];
        S log_identity[N] = {0};
        log(log_ginvg, inverse(g) * g);
        CHECK_ERROR<G>(max_abs_diff(log_ginvg, log_identity, N), max_err*10, "log(inv(a)*a)");

        G gginv;
        multiply_a_binv(gginv, g, g);
        log(log_ginvg, gginv);
        CHECK_ERROR<G>(max_abs_diff(log_ginvg, log_identity, N), max_err*10, "log(a*inv(a))");
    }
    
    {
        S adj[N*N];
        adjoint(adj, g);
        //print_mat(cerr, adj, N, N) << endl;

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
        S err = max_abs_diff(adj_v, adj_v_ref, N);
        CHECK_ERROR<G>(err, max_err*10, "adj");

        adjoint_T_multiply(adj_T_v, g, v);
        err = max_abs_diff(adj_T_v, adj_T_v_ref, N);
        CHECK_ERROR<G>(err, max_err*10, "adjT");
        
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
        
        err = max_abs_diff(log_delta, log_identity, N);
        if (err > max_big_err*10)
        {
            cerr.precision(19);            
            cerr << g << endl;
            print_vec(cerr,v,N) << endl << endl;
            cerr << exp_v << endl;
            cerr << delta << endl;
            cerr << conj << endl;
            cerr << exp_adj_v << endl;
        }
        CHECK_ERROR<G>(err, max_big_err*10, "adj mult");
    }
    
    {     
        S x[D];
        random_vec<D>(x);

        S y[D], z[D];
        transform_point(y, g, x);
        transform_point_by_inverse(z, g, y);        
        S err = max_abs_diff(z, x, D);
        CHECK_ERROR<G>(err, max_err*10, "transform");
        transform_point(z, inverse(g), y);
        err = max_abs_diff(z, x, D);
        CHECK_ERROR<G>(err, max_err*10, "transform inv");
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
    S err = max_abs_diff(c, b, 3);
    if (err > max_err) {
        cerr << "a = ";
        print_vec(cerr, a, 3) << endl;
        cerr << "b = ";
        print_vec(cerr, b, 3) << endl;
        cerr << "R =\n" << R << endl;
        cerr << "R*a = ";
        print_vec(cerr, c, 3) << endl;
    }
    CHECK_ERROR<SO3<S> >(err, max_err, "rotv2v");
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
    
    //for (int pass=0; pass<passes; ++pass) test_group<SL3<float> >();
    for (int pass=0; pass<passes; ++pass) test_group<SL3<double> >();
    
    return 0;
}
