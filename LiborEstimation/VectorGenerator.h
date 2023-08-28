#pragma once

#include "config.h"

#include <vector>

#include "util.h"

template <class BaseGenerator>
class VectorGenerator {
    using Value = std::vector<decltype(make_value<typename std::decay<BaseGenerator>::type>()())>;

    std::vector<BaseGenerator> baseGenerators_;

public:
    VectorGenerator(std::vector<BaseGenerator> baseGenerators)
        : baseGenerators_(std::move(baseGenerators)) {}

    Value operator()() {
        Value result;
        result.reserve(baseGenerators_.size());

        for (auto &generator : baseGenerators_) {
            result.push_back(generator());
        }

        return result;
    }
};

template<class BaseGenerator>
VectorGenerator<BaseGenerator>
makeVectorGenerator(std::vector<BaseGenerator> baseGenerators) {
    return VectorGenerator<BaseGenerator>(std::move(baseGenerators));
}
