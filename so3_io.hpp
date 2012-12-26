#pragma once

#include <iosfwd>

namespace liegroups {

    template <class S> struct SO3;
    
    template <class S>
    std::ostream &operator<<(std::ostream &out, const SO3<S> &g);

    template <class S>
    std::istream &operator>>(std::istream &in, SO3<S> &g);    
}
