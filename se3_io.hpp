#pragma once

#include <iosfwd>

namespace liegroups {

    template <class S> struct SE3;
    
    template <class S>
    std::ostream &operator<<(std::ostream &out, const SE3<S> &g);

    template <class S>
    std::istream &operator>>(std::istream &in, SE3<S> &g);    
}
