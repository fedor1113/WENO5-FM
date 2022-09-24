#ifndef TDRK3_5_H
#define TDRK3_5_H

#include <algorithm>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimestepTDRK3_5(
		std::ranges::common_range auto& u,
		std::ranges::common_range auto& dflux,
		std::ranges::common_range auto& ddflux,
		std::ranges::common_range auto& y2,
		std::ranges::common_range auto& dy2,
		std::ranges::common_range auto& y3,
		std::ranges::common_range auto& dy3,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_size,
		auto&& calcdSpace,
		auto&& calcddSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* Two-derivative 3-stage 5th-order Runge–Kutta
	 * method to discretize a modified method-of-lines (MMOL) ODE
	 * du/dt = L[u], where L is some spatial operator
	 * and L'[u] is its derivative. See Seal et al., 2014.
	 */

	// std::slice Nint(3, nSize, 1);
	dflux.resize(std::ranges::size(u));

	auto iv1 = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(std::ranges::size(u))
	);
	std::vector<std::size_t> iv(std::ranges::size(u));
	for (auto idx : iv1) iv[idx] = idx;

	// ------------------------First Stage----------------------------
	dflux = calcdSpace(u, t, dx, max_eigenvalues,
		n_size, opts_args...);
	ddflux = calcddSpace(u, dflux, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv), std::ranges::end(iv),
				[dt, &y2, &u, &dflux, &ddflux](std::size_t k) {
		y2[k] = u[k] + (40. * dt * dflux[k]
						+ 8. * dt * dt * ddflux[k]) * 0.01;
	});

	updateGhostPoints(y2);

	// ------------------------Second Stage---------------------------
	y3 = calcdSpace(y2, t, dx, max_eigenvalues,
		n_size, opts_args...);
	dy2 = calcddSpace(y2, y3, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv), std::ranges::end(iv),
				[dt, &y3, &u, &dflux, &ddflux, &y2, &dy2](std::size_t k) {
		y3[k] = (u[k] + dt * dflux[k] + (
					- 25. * dt * dt * (ddflux[k] * dflux[k])
					+ 75. * dt * dt * (dy2[k] * y2[k])) * 0.01);
	});

	updateGhostPoints(y3);


	// ------------------------Third Stage----------------------------
	y2 = calcdSpace(y3, t, dx, max_eigenvalues,
					n_size, opts_args...);
	dy3 = calcddSpace(y3, y2, t, dx, max_eigenvalues,
					  n_size, opts_args...);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv), std::ranges::end(iv),
				[dt, &u, &dflux, &ddflux, &dy2, &dy3](
					std::size_t k) {
		u[k] += dt * dflux[k] + dt * dt * (ddflux[k] * 0.125
			+ dy2[k] * (25. / 72.)
			+ dy3[k] * (1. / 36.));
	});

	updateGhostPoints(u);

	// return U;
	// U = res;
}


#endif // TDRK3_5_H
