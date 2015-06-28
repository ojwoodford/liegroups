#include "spline_io.hpp"
#include <iostream>
#include "spline.hpp"
#include "se2.hpp"
#include "se2_io.hpp"
#include "so3.hpp"
#include "so3_io.hpp"
#include "se3.hpp"
#include "se3_io.hpp"

template <typename S>
static std::ostream &print_vec(std::ostream &out, const S x[], int n)
{
    const int w = out.precision() + 8;
    for (int i=0; i<n; ++i) {
        out.width(w);
        out << x[i];
    }
    return out;
}

template <class G>
std::ostream &liegroups::operator<<(std::ostream &out, const QuinticSplineControlPoint<G> &cp)
{
    out << cp.t << std::endl
        << cp.y << std::endl;
    print_vec(out, cp.dy, G::DoF) << std::endl;
    print_vec(out, cp.d2y, G::DoF) << std::endl;
    
    return out;
}

template std::ostream &liegroups::operator<< <liegroups::SE2<float> >(std::ostream&, const QuinticSplineControlPoint<SE2<float> >&);
template std::ostream &liegroups::operator<< <liegroups::SE2<double> >(std::ostream&, const QuinticSplineControlPoint<SE2<double> >&);
template std::ostream &liegroups::operator<< <liegroups::SO3<float> >(std::ostream&, const QuinticSplineControlPoint<SO3<float> >&);
template std::ostream &liegroups::operator<< <liegroups::SO3<double> >(std::ostream&, const QuinticSplineControlPoint<SO3<double> >&);
template std::ostream &liegroups::operator<< <liegroups::SE3<float> >(std::ostream&, const QuinticSplineControlPoint<SE3<float> >&);
template std::ostream &liegroups::operator<< <liegroups::SE3<double> >(std::ostream&, const QuinticSplineControlPoint<SE3<double> >&);

template <class G>
std::istream &liegroups::operator>>(std::istream &in, QuinticSplineControlPoint<G> &cp)
{
    in >> cp.t >> cp.y;
    for (int i=0; i<G::DoF; ++i) {
        in >> cp.dy[i];
    }
    for (int i=0; i<G::DoF; ++i) {
        in >> cp.d2y[i];
    }
    
    return in;
}

template std::istream &liegroups::operator>> <liegroups::SE2<float> >(std::istream&, QuinticSplineControlPoint<SE2<float> >&);
template std::istream &liegroups::operator>> <liegroups::SE2<double> >(std::istream&, QuinticSplineControlPoint<SE2<double> >&);
template std::istream &liegroups::operator>> <liegroups::SO3<float> >(std::istream&, QuinticSplineControlPoint<SO3<float> >&);
template std::istream &liegroups::operator>> <liegroups::SO3<double> >(std::istream&, QuinticSplineControlPoint<SO3<double> >&);
template std::istream &liegroups::operator>> <liegroups::SE3<float> >(std::istream&, QuinticSplineControlPoint<SE3<float> >&);
template std::istream &liegroups::operator>> <liegroups::SE3<double> >(std::istream&, QuinticSplineControlPoint<SE3<double> >&);
