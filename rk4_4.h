#ifndef RK4_4_H
#define RK4_4_H

#include <algorithm>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimestepRK4_4(
		std::ranges::common_range auto& U,
		std::ranges::common_range auto& dflux,
		std::ranges::common_range auto& Y2,
		std::ranges::common_range auto& Y3,
		std::ranges::common_range auto& Y4,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_ghost_points,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* 4th Order 4 Stage Explicit Runge-Kutta Scheme
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 *
	 * From the original Shu-Osher paper on the efficient
	 * implementation of ENO schemes, 1998.
	 */

	const std::size_t n_size = std::ranges::size(U)
			- 2 * n_ghost_points;

	dflux.resize(std::ranges::size(U));
	Y2.resize(std::ranges::size(U));
	Y3.resize(std::ranges::size(U));
	Y4.resize(std::ranges::size(U));

//	auto interior_view = std::views::drop(n_ghost_points)
//			| std::views::take(n_size)
//			| std::ranges::views::common;
//	auto interior_view = std::ranges::views::common;

	// ------------------------First Stage----------------------------
	dflux = calcdSpace(U, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(U/* | interior_view*/),
				std::ranges::end(U/* | interior_view*/),
				std::ranges::begin(dflux/* | interior_view*/),
				std::ranges::begin(Y2/* | interior_view*/),
				[dt](const auto u, const auto df) {
		return u + dt * df * 0.5;
	});

	updateGhostPoints(Y2);


	// ------------------------Second Stage---------------------------
	dflux = calcdSpace(Y2, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(Y2/* | interior_view*/),
				std::ranges::end(Y2/* | interior_view*/),
				std::ranges::begin(dflux/* | interior_view*/),
				std::ranges::begin(Y3/* | interior_view*/),
				[dt](const auto u, const auto df) {
		return u + dt * df * 0.5;
	});

	updateGhostPoints(Y3);


	// ------------------------Third Stage----------------------------
	dflux = calcdSpace(Y3, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(Y3/* | interior_view*/),
				std::ranges::end(Y3/* | interior_view*/),
				std::ranges::begin(dflux/* | interior_view*/),
				std::ranges::begin(Y4/* | interior_view*/),
				[dt](const auto u, const auto df) {
		return u + dt * df;
	});

	updateGhostPoints(Y4);


	// ------------------------Fourth Stage---------------------------
	dflux = calcdSpace(Y4, t, dx, max_eigenvalues,
		n_size, opts_args...);

	auto iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(std::ranges::size(U))
	);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, &U, &Y2, &Y3, &Y4, &dflux](std::size_t k) {
		U[k] = (-U[k] + Y2[k] + 2. * Y3[k] + Y4[k]) * (1. / 3.)
				+ dt * dflux[k] * (1. / 6.);
	});

	updateGhostPoints(U);

	// return U;
	// U = res;
}

#endif // RK4_4_H
