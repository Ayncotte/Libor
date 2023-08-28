#pragma once

#include "config.h"

#include <cassert>

#include "util.h"

template<class BaseGenerator>
class EveryNthGenerator {
    using Value = decltype(make_value<typename std::decay<BaseGenerator>::type>()());

    const std::size_t n_;
    BaseGenerator baseGenerator_;

public:
    EveryNthGenerator(std::size_t n, BaseGenerator baseGenerator)
        : n_(n), baseGenerator_(std::move(baseGenerator)) {
        assert(n > 0);
    }

    Value operator()() {
        for (std::size_t i = 1; i < n_; ++i) {
            baseGenerator_();
        }
        return baseGenerator_();
    }
};

template<class BaseGenerator>
EveryNthGenerator<BaseGenerator>
makeEveryNthGenerator(std::size_t n, BaseGenerator baseGenerator) {
    return EveryNthGenerator<BaseGenerator>(n, std::move(baseGenerator));
}
