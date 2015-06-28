#pragma once

#include <iosfwd>

namespace liegroups {

    template <class S> struct Sim3;
    
    template <class S>
    std::ostream &operator<<(std::ostream &out, const Sim3<S> &g);

    template <class S>
    std::istream &operator>>(std::istream &in, Sim3<S> &g);    
}
