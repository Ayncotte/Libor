#pragma once

#include "config.h"

#include <cassert>

template<class Fun>
auto integrate(Fun fun, double a, double b, std::size_t nranges) -> decltype(fun(a)) {
	assert(a <= b);
	assert(nranges >= 1);

	auto result = 0.5 * (fun(a) + fun(b));

	for (std::size_t i = 1; i < nranges; ++i) {
		result += fun(a + (b - a) * i / nranges);
		assert(std::abs(fun(a + (b - a) * i / nranges).real()) < 100000);
	};

	return result * (b - a) / double(nranges);
}
