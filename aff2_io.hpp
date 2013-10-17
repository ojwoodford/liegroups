#pragma once

#include <iosfwd>

namespace liegroups {

    template <class S> struct Aff2;
    
    template <class S>
    std::ostream &operator<<(std::ostream &out, const Aff2<S> &g);

    template <class S>
    std::istream &operator>>(std::istream &in, Aff2<S> &g);    
}
