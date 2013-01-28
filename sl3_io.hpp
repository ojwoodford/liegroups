#pragma once

#include <iosfwd>

namespace liegroups {

    template <class S> struct SL3;
    
    template <class S>
    std::ostream &operator<<(std::ostream &out, const SL3<S> &g);

    template <class S>
    std::istream &operator>>(std::istream &in, SL3<S> &g);    
}
