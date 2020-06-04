#include <nori/better_integral.h>

NORI_NAMESPACE_BEGIN

double BetterIntegral::Riemann1D(const std::function<double(double)>& f, double x0, double x1, int grid) {
	double space = (x1 - x0) / grid;
	double sum = 0;
	double x = x0;
	for (auto i = 0; i < grid; ++i) {
		sum += f(x) * space;
		x += space;
	}
	return sum;
}

double BetterIntegral::Riemann2D(const std::function<double(double, double)>& f, double x0, double y0, double x1, double y1, int grid)
{
	auto integrate = [&] (double y) {
		return Riemann1D(std::bind(f, std::placeholders::_1, y), x0, x1, grid);
	};
	double value = Riemann1D(integrate, y0, y1, grid);
	return value;
}

NORI_NAMESPACE_END