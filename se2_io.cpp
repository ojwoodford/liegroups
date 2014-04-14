#include <liegroups/se2_io.hpp>
#include <liegroups/se2.hpp>
#include <iostream>

template <class S>
std::ostream &liegroups::operator<<(std::ostream &out, const SE2<S> &g)
{
    const int w = out.precision() + 6;
    out.width(w); out << g.r[0];
    out.width(w); out << g.r[1];
    out.width(w); out << g.t[0] << std::endl;
    out.width(w); out << -g.r[1];
    out.width(w); out << g.r[0];
    out.width(w); out << g.t[1] << std::endl;
    return out;
}

template std::ostream &liegroups::operator<<<float>(std::ostream &, const SE2<float> &);
template std::ostream &liegroups::operator<<<double>(std::ostream &, const SE2<double> &);

template <class S>
std::istream &liegroups::operator>>(std::istream &in, SE2<S> &g)
{
    S nr1, r0;
    in >> g.r[0] >> g.r[1] >> g.t[0];
    in >> nr1 >> r0 >> g.t[1];
    return in;
}

template std::istream &liegroups::operator>><float>(std::istream &, SE2<float> &);
template std::istream &liegroups::operator>><double>(std::istream &, SE2<double> &);
 
