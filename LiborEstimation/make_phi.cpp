#include "make_phi.h"

#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "libalglib/specialfunctions.h"

#include "Model.h"
#include "integrate.h"
#include "util.h"

namespace ublas = boost::numeric::ublas;

namespace {

ScalarComplexFunction make_a(double kappa, double normBeta, double epsilon,
                             double rho) {
    return [=](std::complex<double> z) {
        return checkFinite(kappa - complex_i * z * normBeta * epsilon * rho);
    };
}

ScalarComplexFunction make_d(ScalarComplexFunction a, double normBeta,
                             double epsilon) {
    return [=](std::complex<double> z) {
        return checkFinite(std::sqrt(
            a(z) * a(z) +
            normBeta * normBeta * (complex_i * z + z * z) * epsilon * epsilon));
    };
}

ScalarComplexFunction make_g(ScalarComplexFunction a, ScalarComplexFunction d) {
    return [=](std::complex<double> z) {
        return checkFinite((a(z) + d(z)) / (a(z) - d(z)));
    };
}

std::function<std::complex<double>(std::complex<double>, double)>
make_a_tilde(ScalarComplexFunction a, ScalarComplexFunction d,
             ScalarComplexFunction g, double kappa, double theta,
             double epsilon) {
    std::complex<double> kte = kappa * theta / (epsilon * epsilon);
    return [=](std::complex<double> z, double T) {
        return checkFinite(
            kte * (T * (a(z) - d(z)) -
                   2. * std::log((std::exp(-d(z) * T) - g(z)) / (1. - g(z)))));
    };
}

std::function<std::complex<double>(std::complex<double>, double)>
make_B(ScalarComplexFunction a, ScalarComplexFunction d,
       ScalarComplexFunction g, double epsilon) {
    return [=](std::complex<double> z, double T) {
        const std::complex<double> edt = std::exp(-d(z) * T);
        return checkFinite((a(z) + d(z)) * (edt - 1.) /
                           (epsilon * epsilon * (edt - g(z))));
    };
}

ScalarComplexFunction make_phi(
    std::function<std::complex<double>(std::complex<double>, double)> A_tilde,
    std::function<std::complex<double>(std::complex<double>, double)> B,
    double v, double gamma, double T) {
    return [=](std::complex<double> z) {
        return checkFinite(std::exp(A_tilde(z, T) + v * B(z, T)) *
                           std::exp(-.5 * (complex_i * z + z * z) * gamma));
    };
}

double computeKappa(const Model &model, std::size_t p) {
    double result = model.kappa();

    for (std::size_t j = p; j < model.N(); ++j) {
        result -= model.delta() * (model.F_init()[0] + model.alpha()[j - 1]) /
                  (1. + model.delta() * model.F_init()[0]) *
                  inner_prod(model.sigma(), column(model.beta(), j - 1));
    }

    return result;
}

double computeTheta(const Model &model, double kappa) {
    return model.kappa() * model.theta() / kappa;
}

double computeGamma(const Model &model, std::size_t p, double T) {
    return norm_2(column(model.gamma(), p - 2)) *
           norm_2(column(model.gamma(), p - 2)) * T;
}

double computeNormBeta(const Model &model, std::size_t p) {
    return norm_2(column(model.beta(), p - 2));
}

double computeEpsilon(const Model &model) {
    return std::sqrt(sqr(norm_2(model.sigma())) + sqr(model.sigma_b()));
}

double computeRho(const Model &model, std::size_t p, double epsilon) {
    return inner_prod(model.sigma(), column(model.beta(), p - 2)) /
           norm_2(column(model.beta(), p - 2)) / epsilon;
}

} // anonymous namespace

ScalarComplexFunction make_phi(const Model &model, std::size_t p) {
    auto T = (p - 1) * model.delta();
    auto kappa = computeKappa(model, p);
    auto gamma = computeGamma(model, p, T);
    auto theta = computeTheta(model, kappa);
    auto v = theta;
    auto normBeta = computeNormBeta(model, p);
    auto epsilon = computeEpsilon(model);
    auto rho = computeRho(model, p, epsilon);
    auto a = make_a(kappa, normBeta, epsilon, rho);
    auto d = make_d(a, normBeta, epsilon);
    auto g = make_g(a, d);
    auto B = make_B(a, d, g, epsilon);
    auto A_tilde = make_a_tilde(a, d, g, kappa, theta, epsilon);
    return make_phi(A_tilde, B, v, gamma, T);
}
