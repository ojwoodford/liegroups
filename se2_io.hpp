#pragma once

#include <iosfwd>

namespace liegroups {

    template <class S> struct SE2;
    
    template <class S>
    std::ostream &operator<<(std::ostream &out, const SE2<S> &g);

    template <class S>
    std::istream &operator>>(std::istream &in, SE2<S> &g);    
}
