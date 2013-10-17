#include <liegroups/se2.hpp>
#include <liegroups/se2_io.hpp>
#include <liegroups/so3.hpp>
#include <liegroups/so3_io.hpp>
#include <liegroups/se3.hpp>
#include <liegroups/se3_io.hpp>

#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>

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


template <class G>
void update_mean(const std::vector<G> &samples, G &mu)
{
    const int N = G::DoF;

    double sum[N];
    for (int k=0; k<N; ++k)
        sum[k] = 0.0;
    
    for (size_t i=0; i<samples.size(); ++i) {
        double eps[N];
        log(eps, samples[i] * inverse(mu));

        for (int k=0; k<N; ++k)
            sum[k] += eps[k];
    }

    const double inv_n = 1.0 / (double)samples.size();
    for (int k=0; k<N; ++k)
        sum[k] *= inv_n;

    G delta;
    exp(delta, sum);
    mu = delta * mu;
    rectify(mu);
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

int main(int argc, char *argv[])
{
    typedef SO3<double> G;
    const int N = G::DoF;
    
    double eps[N];
    random_vec<N>(eps);
    
    G mu;    
    exp(mu, eps);

    std::vector<G> samples(100);
    for (size_t i=0; i<samples.size(); ++i) {
        random_vec<N>(eps);
        G delta;
        exp(delta, eps);
        samples[i] = delta * mu;        
    }

    mu = samples[0];
    for (int pass=0; pass<10; ++pass) {
        update_mean(samples, mu);
    }

    double cov[N*N];
    for (int k=0; k<N*N; ++k)
        cov[k] = 0.0;

    for (size_t i=0; i<samples.size(); ++i) {
        double eps[N];
        log(eps, samples[i] * inverse(mu));

        for (int j=0; j<N; ++j) {
            for (int k=j; k<N; ++k) {
                cov[j*N+k] += eps[j] * eps[k];
            }
        }
    }

    const double inv_n = 1.0 / (double)samples.size();
    for (int j=0; j<N; ++j) {
        for (int k=j; k<N; ++k) {
            cov[j*N+k] *= inv_n;
            cov[k*N+j] = cov[j*N+k];
        }
    }

    cout.precision(15);
    cout << mu << endl;
    
    print_mat(cout, cov, N, N);
    
    return 0;
}
