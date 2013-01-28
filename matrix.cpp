#include <liegroups/matrix_impl.hpp>

template bool liegroups::LU_decompose<3,float>(float[], int[]);
template bool liegroups::LU_decompose<3,double>(double[], int[]);

template bool liegroups::sqrtm<3,float>(float[], const float[], const float);
template bool liegroups::sqrtm<3,double>(double[], const double[], const double);

template bool liegroups::expm<3,float>(float[], const float[]);
template bool liegroups::expm<3,double>(double[], const double[]);

template bool liegroups::logm<3,float>(float[], const float[]);
template bool liegroups::logm<3,double>(double[], const double[]);
