#include <liegroups/sl3_io.hpp>
#include <liegroups/sl3.hpp>
#include <iostream>

template <class S>
std::ostream &liegroups::operator<<(std::ostream &out, const SL3<S> &g)
{
    const int w = out.precision() + 6;
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            out.width(w);
            out << g.H[i*3+j];
        }
        out << std::endl;
    }
    return out;
}

template std::ostream &liegroups::operator<<<float>(std::ostream &, const SL3<float> &);
template std::ostream &liegroups::operator<<<double>(std::ostream &, const SL3<double> &);

template <class S>
std::istream &liegroups::operator>>(std::istream &in, SL3<S> &g)
{
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            in >> g.H[i*3+j];
        }
    }
    rectify(g);
    return in;
}

template std::istream &liegroups::operator>><float>(std::istream &, SL3<float> &);
template std::istream &liegroups::operator>><double>(std::istream &, SL3<double> &);
 
