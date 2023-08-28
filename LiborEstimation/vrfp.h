#pragma once

#include "config.h"

#include "util.h"

double curlyB(double F, double T, double sigma, double K);

double vrfp(double B, double F, double T, double sigma, double K,
            ScalarComplexFunction phi);
