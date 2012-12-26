#include <liegroups/matrix.hpp>
#include <cstdlib>
#include <iostream>

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
std::ostream &print_mat(std::ostream &out, const S x[], int m, int n)
{
    const int w = out.precision() + 6;
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            out.width(w);
            out << x[i*n + j];
        }
        out << std::endl;
    }
    return out;
}

template <int N, typename S>
bool sqrtm(S out[N*N], const S in[N*N])
{
    S m[N*N], y[N*N];
    for (int i=0; i<N*N; ++i) {
        m[i] = y[i] = in[i];
    }

    return false;
}

template <int N, typename S>
bool invert(S invm[N*N], const S m[N*N], S &detm)
{
    for (int i=0; i<N*N; ++i)
        invm[i] = (S)0;
    
    for (int i=0; i<N; ++i)
        invm[i*(N+1)] = (S)1;
    
    S mcopy[N*N];
    for (int i=0; i<N*N; ++i)
        mcopy[i] = m[i];
    
    if (!liegroups::solve<N,N>(mcopy, invm))
        return false;

    detm = mcopy[0];
    for (int i=1; i<N; ++i)
        detm *= mcopy[i*(N+1)];

    return true;
}

int main(int argc, char *argv[])
{
    if (argc > 1)
        std::srand(atoi(argv[1]));
    typedef double S;
    std::cout.precision(19);
    
    S w[3];
    random_vec<3>(w, (S)0.1);
    S wx[3*3] = {0, w[0], w[1], -w[0], 0, w[2], -w[1], -w[2], 0};
    S exp_wx[3*3];
    liegroups::expm<3>(exp_wx, wx);

    print_mat(std::cout, wx, 3, 3) << std::endl;
    print_mat(std::cout, exp_wx, 3, 3) << std::endl;

    S sm[3*3];
    liegroups::sqrtm<3>(sm, exp_wx);
    print_mat(std::cout, sm, 3, 3) << std::endl;

    S lnm[3*3];
    if (!liegroups::logm<3>(lnm, exp_wx)) {
        std::cerr << "logm failed" << std::endl;
        return 1;
    }
    print_mat(std::cout, lnm, 3, 3) << std::endl;
    
    S invm[3*3];
    S detm;
    invert<3>(invm, exp_wx, detm);

    std::cout << "det = " << detm << std::endl;
    
    return 0;    
}
         
