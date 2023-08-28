#include "HestonGenerator.h"

#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace ublas = boost::numeric::ublas;

ublas::vector<double> HestonGenerator::operator()() {
    std::normal_distribution<double> distribution(0, std::sqrt(delta_t_));

    ublas::vector<double> F_new(model_.N());
    F_new[0] = F_old_[0];

    ublas::vector<double> W(model_.beta().size1());
    for (std::size_t i = 0; i < W.size(); ++i) {
        W[i] = distribution(generator_);
    }

    ublas::vector<double> W_hat(model_.gamma().size1());
    for (std::size_t i = 0; i < W_hat.size(); ++i) {
        W_hat[i] = distribution(generator_);
    }

    double W_bar = distribution(generator_);
    for (std::size_t k = 1; k < model_.N(); ++k) {
        double drift = 0;
        for (std::size_t i = k + 1; i < model_.N(); ++i) {
            drift += model_.delta() * (F_old_[i] + model_.alpha()[i]) /
                     (1.0 + model_.delta() * F_old_[i]) *
                     (CIR_old_ * inner_prod(column(model_.beta(), k - 1),
                                            column(model_.beta(), i - 1)) +
                      inner_prod(column(model_.gamma(), k - 1),
                                 column(model_.gamma(), i - 1))) *
                     delta_t_;
        }

        F_new[k] =
            F_old_[k] +
            (F_old_[k] + model_.alpha()[k]) *
                (std::sqrt(CIR_old_) *
                     inner_prod(column(model_.beta(), k - 1), W) +
                 inner_prod(column(model_.gamma(), k - 1), W_hat) - drift);
    }

    auto CIR_new = std::max(
        CIR_old_ + model_.kappa() * (model_.theta() - CIR_old_) * delta_t_ +
            sqrt(CIR_old_) *
                (inner_prod(model_.sigma(), W) + model_.sigma_b() * W_bar),
        0.0);

    F_old_ = F_new;
    CIR_old_ = CIR_new;

    return F_new;
}
