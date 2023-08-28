#pragma once

#include "config.h"

#include <cassert>
#include <complex>
#include <functional>

template<class T>
T make_value();

template <class T>
inline T sqr(T x) {
    return x * x;
}

const std::complex<double> complex_i(0, -1);
const double eps = 1e-3;

template <class T>
inline std::complex<double> checkFinite(const std::complex<T> &v) {
    assert(std::isfinite(v.real()));
    assert(std::isfinite(v.imag()));
    return v;
}

using ScalarComplexFunction = std::function<std::complex<double>(std::complex<double>)>;
