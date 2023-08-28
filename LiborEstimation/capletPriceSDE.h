#pragma once

#include "config.h"

#include <vector>

class Model;

std::vector<std::vector<double>> capletPriceSDE(const Model &model, const std::vector<double> &strikes);
