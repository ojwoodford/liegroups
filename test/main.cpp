#include <liegroups/scalar.hpp>

#include <liegroups/se2.hpp>
#include <liegroups/se2_io.hpp>
#include <liegroups/sim2.hpp>
#include <liegroups/sim2_io.hpp>
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
    const int w = out.precision() + 6;
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
}

template <class G>
void test_group()
{
    typedef typename G::Scalar S;
    const int N = G::DoF;
    const int D = G::Dim;
    const S max_err = (S)10 * Constants<S>::epsilon();
    const S max_big_err = Constants<S>::sqrt_epsilon();
    //const S max_err_sq = max_err*max_err;
    
    G g = G::identity;

    {
        S log_g[N];
        random_vec<N>(log_g);

        exp(g, log_g);
        
        S log_exp[N] = { (S)-9999999 };
        log(log_exp, g);
        S err = max_abs_diff(log_exp, log_g, N);
        if (err >= max_big_err) {
            cerr.precision(19);
            cerr << g << endl;
            for (int i=0; i<N; ++i)
                log_exp[i] -= log_g[i];
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

        S adj_v[N], adj_v_ref[N], v[N];
        random_vec<N>(v);

        for (int i=0; i<N; ++i) {
            S ai = 0;
            for (int j=0; j<N; ++j) {
                ai += adj[i*N + j] * v[j];
            }
            adj_v_ref[i] = ai;
        }
        
        adjoint_multiply(adj_v, g, v);
        S err = max_abs_diff(adj_v, adj_v_ref, N);
        CHECK_ERROR<G>(err, max_err*10, "adj");

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
        if (err > max_big_err)
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
    const int passes = 1000000;
    for (int pass=0; pass<passes; ++pass)
    {
        test_rotv2v<float>();
        test_rotv2v<double>();
    }
    check_count = 0;
    for (int pass=0; pass<passes; ++pass)
    {
        test_group<SE2<float> >();
        test_group<SE2<double> >();
        test_group<Sim2<float> >();
        test_group<Sim2<double> >();
        test_group<SL3<float> >();
        test_group<SL3<double> >();
        test_group<SO3<float> >();
        test_group<SO3<double> >();
        test_group<SE3<float> >();
        test_group<SE3<double> >();
        ++check_count;
    }


    
}
