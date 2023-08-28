#include "config.h"

#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "Model.h"
#include "createTestModel.h"
#include "capletPriceSDE.h"
#include "capletPriceVRFP.h"
#include "vrfp.h"

#include <iostream>

namespace ublas = boost::numeric::ublas;

BOOST_AUTO_TEST_CASE(vrfp_test) {
	const double sigma = 1.0;
	const double sigma2 = 2.0;
	const double T = 1.0;
	const auto phi = [&](std::complex<double> z) {
		return std::exp(-0.5 * sqr(sigma2) * T * z * (z + complex_i));
	};
	const double F = 2;
	const double K = 1;
	const double B = 1;

	auto result1 = vrfp(B, F, T, sigma, K, phi);
	auto result2 = curlyB(F, T, sigma2, K);

	std::cout << "vrfp   = " << result1 << std::endl;
	std::cout << "curlyB = " << result2 << std::endl;

	BOOST_CHECK(std::abs(result1 - result2) < eps);
}

template<class T>
void printMatrix(std::ostream &out, const std::vector<std::vector<T>> &matrix) {
    for (const auto &row : matrix) {
        for (const auto &value: row) {
            out << '\t' << value;
        }
        out << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(caplet_price_test) {
    auto model = createTestModel();
	const std::vector<double> strikes = { 0.01, 0.025, 0.05, 0.075 };

    auto result1 = capletPriceSDE(model, strikes);
    auto result2 = capletPriceVRFP(model, strikes);

    BOOST_REQUIRE_EQUAL(result1.size(), strikes.size());
    BOOST_REQUIRE_EQUAL(result2.size(), strikes.size());

    std::cout << "-- SDE result --" << std::endl;
    printMatrix(std::cout, result1);

    std::cout << "-- VRFP result --" << std::endl;
    printMatrix(std::cout, result2);

    for (std::size_t i = 0; i < strikes.size(); ++i) {
        for (std::size_t j = 1; j + 1 < model.N(); ++j) {
            BOOST_CHECK(std::abs(result1[i][j] - result2[i][j - 1]) < eps);
        }
    }
}
