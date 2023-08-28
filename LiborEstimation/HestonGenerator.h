#pragma once

#include "config.h"
#include "Model.h"
#include <random>

class HestonGenerator {
private:
    const Model &model_;
    const double delta_t_;
    std::mt19937 &generator_;
    boost::numeric::ublas::vector<double> F_old_;
    double CIR_old_;

public:
    HestonGenerator(const Model &model, double delta_t, std::mt19937 &generator)
        : model_(model), delta_t_(delta_t), generator_(generator),
          F_old_(model.F_init()), CIR_old_(model.CIR_init()) {}

    boost::numeric::ublas::vector<double> operator()();
};
