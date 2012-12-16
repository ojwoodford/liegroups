#pragma once

#include <iosfwd>

namespace liegroups {

    template <class S> struct Sim2;
    
    template <class S>
    std::ostream &operator<<(std::ostream &out, const Sim2<S> &g);

    template <class S>
    std::istream &operator>>(std::istream &in, Sim2<S> &g);    
}
