#include <liegroups/se3_io.hpp>
#include <liegroups/se3.hpp>
#include <iostream>

template <class S>
std::ostream &liegroups::operator<<(std::ostream &out, const SE3<S> &g)
{
    const int w = out.precision() + 8;
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            out.width(w);
            out << g.R.R[i*3+j];
        }
        out.width(w);
        out << g.t[i];        
        out << std::endl;
    }
    return out;
}

template std::ostream &liegroups::operator<<<float>(std::ostream &, const SE3<float> &);
template std::ostream &liegroups::operator<<<double>(std::ostream &, const SE3<double> &);

template <class S>
std::istream &liegroups::operator>>(std::istream &in, SE3<S> &g)
{
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            in >> g.R.R[i*3+j];
        }
        in >> g.t[i];
    }
    rectify(g);
    return in;
}

template std::istream &liegroups::operator>><float>(std::istream &, SE3<float> &);
template std::istream &liegroups::operator>><double>(std::istream &, SE3<double> &);
 
