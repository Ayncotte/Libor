#include "capletPriceVRFP.h"

#include <cmath>

#include "Model.h"
#include "make_phi.h"
#include "vrfp.h"

std::vector<std::vector<double>> capletPriceVRFP(const Model &model, const std::vector<double> &strikes) {
	const double sigma = 2.0;

    std::vector<std::vector<double>> result(strikes.size());

	for (std::size_t p = 2; p < model.N(); ++p) {
        auto T = (p - 1) * model.delta();
        auto phi = make_phi(model, p);

		const double B = 1 / std::pow(1 + model.delta() * model.F_init()[0], p - 1);

        std::size_t nstrike = 0;
		for (const auto &striKe : strikes) {
			result[nstrike++].push_back(vrfp(B * model.delta(), model.F_init()[0], T, sigma, striKe, phi));
		}
	}

    return result;
}
