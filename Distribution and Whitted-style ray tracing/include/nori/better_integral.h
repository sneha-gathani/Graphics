#pragma once

#include <nori/common.h>

NORI_NAMESPACE_BEGIN

class BetterIntegral {
public:
	BetterIntegral() = delete;

	static double Riemann1D(const std::function<double(double)> &f, double x0, double x1, int grid = 200);

	static double Riemann2D(const std::function<double(double, double)> &f, double x0, double y0,
							double x1, double y1, int grid = 200);
};

NORI_NAMESPACE_END