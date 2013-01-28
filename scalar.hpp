#pragma once

#include <cmath>

namespace liegroups {

    template <class Scalar>
    struct Constants;

    template <class Scalar>
    struct ScalarFunctions;

    template <class S> S abs(S x) { return x < (S)0 ? -x : x; }
    template <class S> S min(S x, S y) { return x < y ? x : y; }
    template <class S> S max(S x, S y) { return x < y ? y : x; }
    
    template <class S> S sqrt(S x) { return ScalarFunctions<S>::sqrt(x); }
    template <class S> S pow(S x, S y) { return ScalarFunctions<S>::pow(x, y); }
    template <class S> S exp(S x) { return ScalarFunctions<S>::exp(x); }
    template <class S> S ln(S x) { return ScalarFunctions<S>::ln(x); }
    template <class S> S sin(S x) { return ScalarFunctions<S>::sin(x); }
    template <class S> S cos(S x) { return ScalarFunctions<S>::cos(x); }
    template <class S> S tan(S x) { return ScalarFunctions<S>::tan(x); }
    template <class S> S asin(S x) { return ScalarFunctions<S>::asin(x); }
    template <class S> S acos(S x) { return ScalarFunctions<S>::acos(x); }
    template <class S> S atan(S x) { return ScalarFunctions<S>::atan(x); }
    template <class S> S atan2(S y, S x) { return ScalarFunctions<S>::atan2(y,x); }
    
    template <>
    struct Constants<float>
    {
        static float epsilon() { return 5.96e-8f; }
        static float sqrt_epsilon() { return 2.44e-4f; }
    };

    template <>
    struct Constants<double>
    {
        static double epsilon() { return 1.11e-16; }
        static double sqrt_epsilon() { return 1.054e-8; }
    };

    template <>
    struct ScalarFunctions<float>
    {
        static float sqrt(float x) { return ::sqrtf(x); }
        static float pow(float x, float y) { return ::powf(x,y); }
        static float exp(float x) { return ::expf(x); }
        static float ln(float x) { return ::logf(x); }
        static float sin(float x) { return ::sinf(x); }
        static float cos(float x) { return ::cosf(x); }
        static float tan(float x) { return ::tanf(x); }
        static float asin(float x) { return ::asinf(x); }
        static float acos(float x) { return ::acosf(x); }
        static float atan(float x) { return ::atanf(x); }
        static float atan2(float y, float x) { return ::atan2f(y,x); }
    };

    template <>
    struct ScalarFunctions<double>
    {
        static double sqrt(double x) { return ::sqrt(x); }
        static double pow(double x, double y) { return ::pow(x,y); }
        static double exp(double x) { return ::exp(x); }
        static double ln(double x) { return ::log(x); }
        static double sin(double x) { return ::sin(x); }
        static double cos(double x) { return ::cos(x); }
        static double tan(double x) { return ::tan(x); }
        static double asin(double x) { return ::asin(x); }
        static double acos(double x) { return ::acos(x); }
        static double atan(double x) { return ::atan(x); }
        static double atan2(double y, double x) { return ::atan2(y,x); }
    };
   
}
