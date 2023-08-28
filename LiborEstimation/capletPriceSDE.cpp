#include "capletPriceSDE.h"

#include <boost/numeric/ublas/vector.hpp>

#include "EveryNthGenerator.h"
#include "HestonGenerator.h"
#include "Model.h"
#include "VectorGenerator.h"

namespace ublas = boost::numeric::ublas;

double computeCaplet(const ublas::vector<double> &F1, const ublas::vector<double> &F2,
              double strike, double F_init, double delta, std::size_t k) {
    assert(F1.size() == F2.size());

    auto N = F1.size();

    double result =
        std::max(F1[k] - strike, 0.0) * delta /
        pow(1 + delta * F_init, N);

    for (std::size_t i = k; i < N; ++i) {
        result *= 1 + F2[i] * delta;
    }

    return result;
}

std::vector<double> computeCaplets(const Model &model, const std::vector<std::vector<ublas::vector<double>>> &Fs, double strike) {
    std::vector<double> result;
    for (std::size_t i = 0; i < Fs.size(); ++i) {
        double caplet = 0.0;
        for (std::size_t j = 0; j < Fs[i].size(); ++j) {
            const ublas::vector<double> &F1 = i == 0 ? model.F_init() : Fs[i - 1][j];
            const ublas::vector<double> &F2 = Fs[i][j];
            caplet += computeCaplet(F1, F2, strike, model.F_init()[0], model.delta(), i);
        }
        result.push_back(caplet / Fs[i].size());
    }
    return result;
}

std::vector<std::vector<double>> capletPriceSDE(const Model &model, const std::vector<double> &strikes) {
    const std::size_t stepsPerInterval = 100;
    const std::size_t intervalsPerYear = 4;
    const std::size_t years = 2;
	const std::size_t nExperiments = 20000;
	const double delta_t = 1.0 / (stepsPerInterval * intervalsPerYear);

	std::mt19937 randomGenerator;
    HestonGenerator heston(model, delta_t, randomGenerator);
    auto everyNthGenerator = makeEveryNthGenerator(stepsPerInterval, std::move(heston));
    auto vectorGenerator = makeVectorGenerator(std::vector<decltype(everyNthGenerator)>(nExperiments, everyNthGenerator));

    std::vector<std::vector<ublas::vector<double>>> Fs;

    for (std::size_t i = 0; i < intervalsPerYear * years; ++i) {
        Fs.push_back(vectorGenerator());
    }

    std::vector<std::vector<double>> result;

    for (const auto &strike : strikes) {
        result.push_back(computeCaplets(model, Fs, strike));
    }

    return result;
}
