#include <liegroups/se3.hpp>
#include <liegroups/se3_io.hpp>

#include <cstdlib>
#include <cassert>
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

int main(int argc, char *argv[])
{
    int n = 10;
    if (argc > 1) {
        n = atoi(argv[1]);
        if (n < 1) {
            cerr << "Argument must be a positive integer" << endl;
            exit(1);
        }
    }
    
    typedef SE3<double> G;
    const int N = G::DoF;
    
    double eps[N];
    random_vec<N>(eps);
    
    G mu;    
    exp(mu, eps);

    cout.precision(16);
    
    for (int i=0; i<n; ++i) {
        random_vec<N>(eps);
        G delta;
        exp(delta, eps);
        G g = delta * mu;
        cout << g << endl;
    }
    return 0;
}
