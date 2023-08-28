#include "vrfp.h"

#include <boost/math/constants/constants.hpp>

#include "libalglib/specialfunctions.h"

#include "integrate.h"

double curlyB(double F, double T, double sigma, double K) {
    const double logFK = std::log(F / K);
    const double ssT = 0.5 * sqr(sigma) * T;
    const double ssqrtT = sigma * std::sqrt(T);

    return (F * alglib::normaldistribution((logFK + ssT) / ssqrtT) -
            K * alglib::normaldistribution((logFK - ssT) / ssqrtT));
}

double vrfp(double B, double F, double T, double sigma, double K,
            ScalarComplexFunction phi) {
    const double logKF = std::log(K / F);
    const double ssT = 0.5 * sqr(sigma) * T;
    const double pi = boost::math::constants::pi<double>();

    auto integrant = [&](double z) {
        if (-eps < z && z <= 0) {
            z = -eps;
        } else if (0 <= z && z <= eps) {
            z = eps;
        }
        return (std::exp(-ssT * z * (z - complex_i)) - phi(z - complex_i)) /
               (z * (z - complex_i)) * std::exp(-complex_i * z * logKF);
    };

    const std::size_t nsamples = 5000;
    const double cutoff = 10;

    auto integral = integrate(integrant, -cutoff, cutoff, nsamples);

    // assert(std::abs(integral.imag()) < eps);

    return B * (curlyB(F, T, sigma, K) + F / (2 * pi) * integral.real());
}
