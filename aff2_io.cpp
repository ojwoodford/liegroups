#include <liegroups/aff2_io.hpp>
#include <liegroups/aff2.hpp>
#include <iostream>

template <class S>
std::ostream &liegroups::operator<<(std::ostream &out, const Aff2<S> &g)
{
    const int w = out.precision() + 6;
    out.width(w); out << g.A[0];
    out.width(w); out << g.A[1];
    out.width(w); out << g.t[0] << std::endl;
    out.width(w); out << g.A[2];
    out.width(w); out << g.A[3];
    out.width(w); out << g.t[1] << std::endl;
    return out;
}

template std::ostream &liegroups::operator<<<float>(std::ostream &, const Aff2<float> &);
template std::ostream &liegroups::operator<<<double>(std::ostream &, const Aff2<double> &);

template <class S>
std::istream &liegroups::operator>>(std::istream &in, Aff2<S> &g)
{
    in >> g.A[0] >> g.A[1] >> g.t[0];
    in >> g.A[2] >> g.A[3] >> g.t[1];
    return in;
}

template std::istream &liegroups::operator>><float>(std::istream &, Aff2<float> &);
template std::istream &liegroups::operator>><double>(std::istream &, Aff2<double> &);
 
