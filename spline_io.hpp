#pragma once

#include <iosfwd>

namespace liegroups {

    template <class G> struct QuinticSplineControlPoint;

    template <class S>
    std::ostream &operator<<(std::ostream &out, const QuinticSplineControlPoint<S> &cp);

    template <class S>
    std::istream &operator>>(std::istream &in, QuinticSplineControlPoint<S> &cp);
}
