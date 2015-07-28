#include <liegroups/exp_coefs.hpp>
#include <liegroups/scalar.hpp>

using namespace liegroups;

template <typename S>
static void compute_exp_coefs_small(S theta_sq, ExpCoefs<S> &coefs)
{
    const S tt = theta_sq;
    coefs.B = (S)0.5 - tt*((S)(1.0/24) - tt*(S)(1.0/720));
    coefs.cos_theta = 1 - tt*coefs.B;
    coefs.C = (S)(1.0/6) - tt*((S)(1.0/120) - tt*(S)(1.0/5040));
    coefs.A = 1 - tt*coefs.C;
}

template <typename S>
static void compute_diff_exp_coefs_small(S theta_sq, S &A1, S &B1, S &C1)
{
    const S tt = theta_sq;
    A1 = (S)(-1.0/3) + tt*((S)(1.0/30) - tt*(S)(1.0/840));
    B1 = (S)(-1.0/12) + tt*((S)(1.0/180) - tt*(S)(1.0/6720));
    C1 = (S)(-1.0/60) + tt*((S)(1.0/1260) - tt*(S)(1.0/60480));
}

template <typename S>
static void compute_diff2_exp_coefs_small(S theta_sq, S &A2, S &B2, S &C2)
{
    const S tt = theta_sq;
    A2 = (S)(1.0/15) - tt*((S)(1.0/210) - tt*(S)(1.0/7560));
    B2 = (S)(1.0/90) - tt*((S)(1.0/1680) - tt*(S)(1.0/75600));
    C2 = (S)(1.0/630) - tt*((S)(1.0/15120) - tt*(S)(1.0/831600));
}

// Returns 1/theta_sq
template <typename S>
S compute_exp_coefs_large(S theta_sq, ExpCoefs<S> &coefs)
{
    const S theta = liegroups::sqrt(theta_sq);
    const S inv_tt = (S)1 / theta_sq;
    coefs.A = liegroups::sin(theta)/theta;
    coefs.cos_theta = liegroups::cos(theta);
    coefs.B = (1 - coefs.cos_theta) * inv_tt;
    coefs.C = (1 - coefs.A) * inv_tt;
    return inv_tt;
 }

template <typename S>
void liegroups::ExpCoefs<S>::compute(S theta_sq)
{
    if (theta_sq < 25 * Constants<S>::sqrt_epsilon()) {
        compute_exp_coefs_small(theta_sq, *this);
    } else {
        compute_exp_coefs_large(theta_sq, *this);
    }
}

template <typename S>
void liegroups::DiffExpCoefs<S>::compute(S theta_sq)
{
    if (theta_sq < 25 * Constants<S>::sqrt_epsilon()) {
        compute_exp_coefs_small(theta_sq, *this);
        compute_diff_exp_coefs_small(theta_sq, A1, B1, C1);
    } else {
        const S inv_tt = compute_exp_coefs_large(theta_sq, *this);
        A1 = this->C - this->B;
        B1 = (this->A - 2*this->B) * inv_tt;
        C1 = (this->B - 3*this->C) * inv_tt;
    }
}


template <typename S>
void liegroups::Diff2ExpCoefs<S>::compute(S theta_sq)
{
    if (theta_sq < 25 * Constants<S>::sqrt_epsilon()) {
        compute_exp_coefs_small(theta_sq, *this);
        compute_diff_exp_coefs_small(theta_sq, this->A1, this->B1, this->C1);
        compute_diff2_exp_coefs_small(theta_sq, A2, B2, C2);
    } else {
        const S inv_tt = compute_exp_coefs_large(theta_sq, *this);
        this->A1 = this->C - this->B;
        this->B1 = (this->A - 2*this->B) * inv_tt;
        this->C1 = (this->B - 3*this->C) * inv_tt;
        A2 = (-this->A - 3*this->A1) * inv_tt;
        B2 = (this->A1 - 4*this->B1) * inv_tt;
        C2 = (this->B1 - 5*this->C1) * inv_tt;
    }
}

template struct ExpCoefs<float>;
template struct ExpCoefs<double>;

template struct DiffExpCoefs<float>;
template struct DiffExpCoefs<double>;

template struct Diff2ExpCoefs<float>;
template struct Diff2ExpCoefs<double>;
