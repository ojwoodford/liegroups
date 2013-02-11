#include <liegroups/so3.hpp>
#include <liegroups/so3_io.hpp>
#include <liegroups/matrix_impl.hpp>
#include <liegroups/scalar.hpp>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <fstream>

using namespace liegroups;

template <class G>
struct Measurement
{
    G z;
    G x;
    typename G::Scalar Rinv;
};

template <int M, int N, typename S, class Op>
void outer_product_T_upper(S out[M*M], const S J[N*M], S scale, const Op &op)
{
    for (int i=0; i<M; ++i) {
        for (int j=i; j<M; ++j) {
            S sum = (S)0;
            for (int k=0; k<N; ++k) {
                sum += J[i+k*M] * J[j + k*M];
            }
            op(out[i*M+j], sum * scale);
        }
    }
}

struct Add
{
    template <typename A, typename B>
    void operator()(A &a, const B &b) const { a += b; }
};

template <class G>
bool refine_extrinsics(const std::vector<Measurement<G> > &meas,
                       G &E,
                       typename G::Scalar &rms,
                       bool print_errors = false)
{
    typedef typename G::Scalar S;
    const int N = G::DoF;
    S A[N*N];
    S b[N];

    for (int i=0; i<N*N; ++i)
        A[i] = (S)0;    
    for (int i=0; i<N; ++i) {
        A[i*(N+1)] = (S)1;
        b[i] = (S)0;
    }

    S residual = (S)0;

    G Einv = inverse(E);
    for (size_t i=0; i<meas.size(); ++i) {
        G pred = Einv * meas[i].x * E;
        G delta = meas[i].z * inverse(pred);
        S v[N];
        log(v, delta);

        S J[N*N];
        adjoint(J, pred);
        for (int j=0; j<N*N; ++j)
            J[j] *= (S)-1;
        for (int j=0; j<N; ++j)
            J[j*(N+1)] += (S)1;

        outer_product_T_upper<N,N>(A, J, meas[i].Rinv, Add());

        S err_sq = (S)0;
        for (int j=0; j<N; ++j) {
            for (int k=0; k<N; ++k) {
                b[k] += J[j*N+k] * (v[j]*meas[i].Rinv);
            }
            err_sq += v[j]*v[j] * meas[i].Rinv;            
        }
        if (print_errors) {            
            std::cout << i << "\t ";
            for (int j=0; j<N; ++j) {
                std::cout.width(20);
                std::cout << v[j] << " ";
            }
            std::cout << std::endl;
        }
        residual += err_sq;
    }

    for (int i=0; i<N; ++i) {
        for (int j=i+1; j<N; ++j) {
            A[j*N+i] = A[i*N+j];
        }
    }
    
    if (print_errors) {
        std::cout << "information =" << std::endl;
        for (int i=0; i<N; ++i) {
            for (int j=0; j<N; ++j) {
                std::cout.width(20);
                std::cout << A[i*N+j] << " ";
            }
            std::cout << std::endl;
        }
    }
    int index[N];
    if (!LU_decompose<N>(A, index))
        return false;

    S update[N];
    LU_inverse_times_vec<N,1>(update, A, index, b);

    G delta;
    exp(delta, update);

    E = inverse(delta * Einv);

    rms = liegroups::sqrt(residual / (S)meas.size());
    
    return true;
}
                      
int main(int argc, char* argv[])
{
    typedef double S;
    typedef SO3<S> G;

    const int N = G::DoF;
    std::vector<Measurement<G> > meas;

    assert(argc == 3);
    std::cout.precision(12);

    {
        std::cout << "Reading refs from " << argv[1] << std::endl;
        std::ifstream in(argv[1]);
        if (!in.good()) {
            std::cerr << "Error reading " << argv[1] << std::endl;
            exit(1);
        }
        Measurement<G> m;
        while (in >> m.z) {
            m.Rinv = (S)1e4;
            meas.push_back(m);
        }
        std::cout << meas.size() << " refs" << std::endl;
    }

    {
        std::cout << "Reading meas from " << argv[2] << std::endl;
        std::ifstream in(argv[2]);
        if (!in.good()) {
            std::cerr << "Error reading " << argv[2] << std::endl;
            exit(1);
        }
        for (size_t i=0; i<meas.size(); ++i) {
            if (!(in >> meas[i].x)) {
                std::cerr << "Error reading meas " << i << std::endl;
                exit(1);
            }
        }
    }
    
    G E = G::identity;
    S rms;

    const int PASSES=16;
    for (int pass=0; pass<PASSES; ++pass) {
        if (!refine_extrinsics(meas, E, rms, pass+1==PASSES)) {
            std::cerr << "refine failed" << std::endl;
            exit(1);
        }
        std::cout << "rms = " << rms << std::endl;
    }

    std::cout << "E =" << std::endl
              << E << std::endl;

    S eps[N];
    log(eps, E);

    std::cout << "log(E) = " << std::endl;
    for (int i=0; i<N; ++i) {
        std::cout.width(25);
        std::cout << eps[i] << " ";
    }
    std::cout << std::endl;
    
    return 0;    
}
