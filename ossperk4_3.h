#ifndef OSSPERK4_3_H
#define OSSPERK4_3_H

#include <algorithm>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimestepoSSPERK4_3(
		std::ranges::common_range auto& u_curr,
		std::ranges::common_range auto& dflux,
		std::ranges::common_range auto& dflux_temp,
		std::ranges::common_range auto& u1,
		std::ranges::common_range auto& u2,
		std::ranges::common_range auto& u3,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_ghost_points,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* Optimal 3rd Order 4 Stage Explicit Total Variation Dimming
	 * / Diminishing (Strong Stability Preserving)
	 * Runge-Kutta Scheme (TVD ERK3 / oSSPERK(4,3))
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 * This method is due to the work of Ruuth and Spitery (2002).
	 *
	 * This Ruuth-Spitery 3rd Order SSP purely multistage method
	 * is provably optimal among its kind
	 * in terms of its effective SSP ('CFL') coefficient:
	 * SSP coefficient = 2., so C_eff = 2./4. = 0.5.
	 * Note that it outperforms oSSPERK(3, 3) in this metric.
	 */

	const std::size_t n_size = std::ranges::size(u_curr)
			- 2 * n_ghost_points;

	dflux.resize(std::ranges::size(u_curr));
	dflux_temp.resize(std::ranges::size(u_curr));

//	auto interior_view = std::views::drop(n_ghost_points)
//			| std::views::take(n_size)
//			| std::ranges::views::common;
//	auto interior_view = std::ranges::views::common;

	// ------------------------First Stage----------------------------
	dflux = calcdSpace(u_curr, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(u_curr),
				std::ranges::end(u_curr),
				std::ranges::begin(dflux),
				std::ranges::begin(u1),
				[dt](const auto u0, const auto df) {
		return u0 + 0.5 * dt * df;
	});

	updateGhostPoints(u1);


	// ------------------------Second Stage---------------------------
	dflux_temp = calcdSpace(u1, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(u1),
				std::ranges::end(u1),
				std::ranges::begin(dflux_temp),
				std::ranges::begin(u2),
				[dt](const auto u1, const auto df) {
		return u1 + 0.5 * dt * df;
	});

	updateGhostPoints(u2);


	// ------------------------Third Stage----------------------------
	dflux_temp = calcdSpace(u2, t, dx, max_eigenvalues,
		n_size, opts_args...);


	auto iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(std::ranges::size(u_curr))
	);
	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv),
				std::ranges::end(iv),
				[dt, &u_curr, &u2, &u3, &dflux_temp](std::size_t k) {
		u3[k] = (4. * u_curr[k]
				 + 2. * u2[k]
				 + dt * dflux_temp[k]) * (1./6.);
	});

	updateGhostPoints(u3);

	// ------------------------Fourth Stage---------------------------
	dflux_temp = calcdSpace(u3, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(u3),
				std::ranges::end(u3),
				std::ranges::begin(dflux_temp),
				std::ranges::begin(u_curr),
				[dt](const auto u3, const auto df3) {
		return u3 + 0.5 * dt * df3;
	});

	updateGhostPoints(u_curr);
}

#endif // OSSPERK4_3_H
