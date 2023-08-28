#pragma once

#include "config.h"

#include <utility>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

class Model {
    std::size_t N_;
    boost::numeric::ublas::vector<double> alpha_;
    double delta_;
    boost::numeric::ublas::vector<double> F_init_;
    double CIR_init_;
    boost::numeric::ublas::matrix<double> beta_;
    boost::numeric::ublas::matrix<double> gamma_;
    double kappa_;
    double theta_;
    boost::numeric::ublas::vector<double> sigma_;
    double sigma_b_;

public:
    Model(std::size_t N, boost::numeric::ublas::vector<double> alpha,
          double delta, boost::numeric::ublas::vector<double> F_init,
          double CIR_init, boost::numeric::ublas::matrix<double> beta,
          boost::numeric::ublas::matrix<double> gamma, double kappa,
          double theta, boost::numeric::ublas::vector<double> sigma,
          double sigma_b)
        : N_(N), alpha_(std::move(alpha)), delta_(delta),
          F_init_(std::move(F_init)), CIR_init_(CIR_init),
          beta_(std::move(beta)), gamma_(std::move(gamma)), kappa_(kappa),
          theta_(theta), sigma_(std::move(sigma)), sigma_b_(sigma_b) {}

    std::size_t N() const { return N_; }
    const boost::numeric::ublas::vector<double> &alpha() const {
        return alpha_;
    }
    double delta() const { return delta_; }
    const boost::numeric::ublas::vector<double> &F_init() const {
        return F_init_;
    }
    double CIR_init() const { return CIR_init_; }
    const boost::numeric::ublas::matrix<double> &beta() const { return beta_; }
    const boost::numeric::ublas::matrix<double> &gamma() const {
        return gamma_;
    }
    double kappa() const { return kappa_; }
    double theta() const { return theta_; }
    const boost::numeric::ublas::vector<double> &sigma() const {
        return sigma_;
    }
    double sigma_b() const { return sigma_b_; }
};
