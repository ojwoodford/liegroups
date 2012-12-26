#include <liegroups/so3_io.hpp>
#include <liegroups/so3.hpp>
#include <iostream>

template <class S>
std::ostream &liegroups::operator<<(std::ostream &out, const SO3<S> &g)
{
    const int w = out.precision() + 6;
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            out.width(w);
            out << g.R[i*3+j];
        }
        out << std::endl;
    }
    return out;
}

template std::ostream &liegroups::operator<<<float>(std::ostream &, const SO3<float> &);
template std::ostream &liegroups::operator<<<double>(std::ostream &, const SO3<double> &);

template <class S>
std::istream &liegroups::operator>>(std::istream &in, SO3<S> &g)
{
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            in >> g.R[i*3+j];
        }
    }
    rectify(g);
    return in;
}

template std::istream &liegroups::operator>><float>(std::istream &, SO3<float> &);
template std::istream &liegroups::operator>><double>(std::istream &, SO3<double> &);
 
