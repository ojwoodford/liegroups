#include <liegroups/sim2_io.hpp>
#include <liegroups/se2_io.hpp>
#include <liegroups/sim2.hpp>
#include <iostream>

template <class S>
std::ostream &liegroups::operator<<(std::ostream &out, const Sim2<S> &g)
{
    const int w = out.precision() + 6;
    out.width(w); out << g.scale << std::endl << g.rigid;
    return out;
}

template std::ostream &liegroups::operator<<<float>(std::ostream &, const Sim2<float> &);
template std::ostream &liegroups::operator<<<double>(std::ostream &, const Sim2<double> &);

template <class S>
std::istream &liegroups::operator>>(std::istream &in, Sim2<S> &g)
{
    if (!(in >> g.scale))
        return in;    
       
    g.inv_scale = (S)1/g.scale;
    return (in >> g.rigid);
}

template std::istream &liegroups::operator>><float>(std::istream &, Sim2<float> &);
template std::istream &liegroups::operator>><double>(std::istream &, Sim2<double> &);
 
