#include <liegroups/se3.hpp>
#include <liegroups/se3_io.hpp>

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace liegroups;

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
    const int w = out.precision() + 8;
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
    std::istream *in = &cin;
    std::ifstream fin;
    bool read_timestamps = true;
    if (argc > 1) {
        fin.open(argv[1]);
        if (!fin.good()) {
            cerr << "Error: Could not open " << argv[1] << endl;
            exit(1);
        }
        in = &fin;
    }

    if (argc > 2) {
        int rt = atoi(argv[2]);
        read_timestamps = (rt != 0);
    }
    
    typedef SE3<double> G;
    const int N = G::DoF;

    G g;    
    std::vector<G> samples;
    if (read_timestamps) {
        double t;
        while (*in >> t >> g) {
            samples.push_back(g);
        }
    } else {
        while (*in >> g) {
            samples.push_back(g);
        }
    }

    if (samples.empty()) {
        cerr << "Error: No samples!" << endl;
        exit(1);
    }
        
    G mu = samples[0];
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

    double root_diag[N];
    for (int i=0; i<N; ++i)
        root_diag[i] = sqrt(cov[i*N+i]);

    cout << endl;
    print_mat(cout, root_diag, 1, N);
    
    return 0;
}
