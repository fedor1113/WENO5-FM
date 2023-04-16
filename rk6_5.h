#ifndef RK6_5_H
#define RK6_5_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimestepRK6_5(
		std::ranges::common_range auto& U,
		std::ranges::common_range auto& L0,
		std::ranges::common_range auto& Y1,
		std::ranges::common_range auto& L1,
		std::ranges::common_range auto& Y2,
		std::ranges::common_range auto& L2,
		std::ranges::common_range auto& Y3,
		std::ranges::common_range auto& L3,
		std::ranges::common_range auto& Y4,
		std::ranges::common_range auto& L4,
		std::ranges::common_range auto& Y5,
		std::ranges::common_range auto& L5,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_ghost_points,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* 5th Order 6 Stage Explicit Runge-Kutta Scheme (ERK(6,5))
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 *
	 * This is an embedded RK(6, 5) / SSP(3, 3) method introduced
	 * in Colin Barr Macdonald's thesis:
	 * see 'Constructing High-Order Runge-Kutta Methods with
	 * Embedded Strong-Stability-Preserving Pairs' (2001)
	 * and successfully employed by Henrick et al. with WENO5-M
	 * to resolve detonation waves.
	 */

	const std::size_t n_size = std::ranges::size(U)
			- 2 * n_ghost_points;

	auto iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(n_ghost_points/*0*/))
				| std::views::take(n_size/*std::ranges::size(U)*/)
	);

	L0.resize(std::ranges::size(U));
	L1.resize(std::ranges::size(U));
	L2.resize(std::ranges::size(U));
	L3.resize(std::ranges::size(U));
	L4.resize(std::ranges::size(U));
	L5.resize(std::ranges::size(U));

//	auto interior_view = std::views::drop(n_ghost_points)
//			| std::views::take(n_size)
//			| std::ranges::views::common;
//	auto interior_view = std::ranges::views::common;

	// ------------------------First Stage----------------------------
	L0 = calcdSpace(U, t, dx, max_eigenvalues,
		n_size, opts_args...);
	updateGhostPoints(L0);

	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(U/* | interior_view*/),
				std::ranges::end(U/* | interior_view*/),
				std::ranges::begin(L0/* | interior_view*/),
				std::ranges::begin(Y1/* | interior_view*/),
				[dt](const auto u, const auto df) {
		return u + dt * df;
	});

	updateGhostPoints(Y1);


	// ------------------------Second Stage-----------------------------
	L1 = calcdSpace(Y1, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::for_each(
//				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, &U, &L0, &L1, &Y2](std::size_t k) {
		Y2[k] =  U[k] + (L0[k] + L1[k]) * 0.25 * dt;
	});

	updateGhostPoints(Y2);


	// ------------------------Third Stage------------------------------
	L2 = calcdSpace(Y2, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::for_each(
//				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, &U, &L0, &L1, &L2, &Y3](std::size_t k) {
		Y3[k] =  U[k] + (2046. * L0[k]
					- 454. * L1[k]
					+ 1533. * L2[k]) * dt / 15625.;
	});

	updateGhostPoints(Y3);


	// ------------------------Fourth Stage-----------------------------
	L3 = calcdSpace(Y3, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::for_each(
//				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, &U, &L0, &L1, &L2, &L3, &Y4](std::size_t k) {
		Y4[k] =  U[k] + (-2217. * L0[k]
					+ 1533. * L1[k]
					- 566. * L2[k]
					+ 12500. * L3[k]) * dt / 16875.;
	});

	updateGhostPoints(Y4);


	// ------------------------Fifth Stage------------------------------
	L4 = calcdSpace(Y4, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::for_each(
//				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, &U, &L0, &L1, &L2, &L3, &L4, &Y5](std::size_t k) {
		Y5[k] =  U[k]
				+ (11822. * L0[k]
					- 6928. * L1[k]
					- 4269. * L2[k]
					- 12500. * L3[k]
					+ 33750. * L4[k]) * dt / 21875.;
	});

	updateGhostPoints(Y5);


	// ------------------------Final Stage------------------------------
	L5 = calcdSpace(Y5, t, dx, max_eigenvalues,
		n_size, opts_args...);

	std::for_each(
//				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, &U, &L0, &L1, &L2, &L3, &L4, &L5](std::size_t k) {
//		auto U_old = U[k];
		U[k] += (14. * L0[k]
					+ 125. * L3[k]
					+ 162. * L4[k]
					+ 35. * L5[k]) * dt / 336.;
//		if (std::abs(std::max(U[k] - U_old) / std::max(U[k])) > 0.5)
//			U[k] += (L0[k] + L1[k] + 4. * L2[k]) * dt / 6.;
	});

	updateGhostPoints(U);
}

#endif // RK6_5_H
