#include <liegroups/matrix_impl.hpp>


namespace liegroups {

#ifndef _MSC_VER
    template <>
    bool invert<2,float>(float invm[2*2], const float m[2*2], float *detm)
    {
        return invert2(invm, m, detm);
    }

    template <>
    bool invert<2,double>(double invm[2*2], const double m[2*2], double *detm)
    {
        return invert2(invm, m, detm);
    }

    template <>
    bool invert<3,float>(float invm[3*3], const float m[3*3], float *detm)
    {
        return invert3(invm, m, detm);
    }

    template <>
    bool invert<3,double>(double invm[3*3], const double m[3*3], double *detm)
    {
        return invert3(invm, m, detm);
    }
#else
    template <>
    bool invert<2,float>(float *invm, const float *m, float *detm)
    {
        return invert2(invm, m, detm);
    }

    template <>
    bool invert<2,double>(double *invm, const double *m, double *detm)
    {
        return invert2(invm, m, detm);
    }

    template <>
    bool invert<3,float>(float *invm, const float *m, float *detm)
    {
        return invert3(invm, m, detm);
    }

    template <>
    bool invert<3,double>(double *invm, const double *m, double *detm)
    {
        return invert3(invm, m, detm);
    }
#endif

}

#ifndef _MSC_VER
template bool liegroups::LU_decompose<3,float>(float[], int[]);
template bool liegroups::LU_decompose<3,double>(double[], int[]);

template bool liegroups::sqrtm<3,float>(float[], const float[], const float);
template bool liegroups::sqrtm<3,double>(double[], const double[], const double);

template bool liegroups::expm<3,float>(float[], const float[]);
template bool liegroups::expm<3,double>(double[], const double[]);

template bool liegroups::logm<3,float>(float[], const float[]);
template bool liegroups::logm<3,double>(double[], const double[]);

template void liegroups::mat_mult_square<3,float>(float *, float const*, float const*);
template void liegroups::mat_mult_square<3,double>(double *, double const*, double const*);

template void liegroups::mat_mult_square<6,float>(float *, float const*, float const*);
template void liegroups::mat_mult_square<6,double>(double *, double const*, double const*);
#endif
