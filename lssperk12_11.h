#ifndef LSSPERK12_11_H
#define LSSPERK12_11_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimesteplSSPERK12_11(
		std::ranges::common_range auto& u,
		std::ranges::common_range auto& dflux,
		std::array<
			std::reference_wrapper<std::valarray<T>
		>, 28> interim_us,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_ghost_points,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* 9th Order 10 Stage Explicit linear Strong Stability Preserving
	 * Runge-Kutta Scheme (lSSPERK(10,9))
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 *
	 * Only SSP for linear problems!!!
	 */

	const std::size_t n_size = std::ranges::size(u)
			- 2 * n_ghost_points;

	auto iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(n_ghost_points/*0*/))
				| std::views::take(n_size/*std::ranges::size(U)*/)
	);

	interim_us[0].get() = u;

	for (std::size_t k : std::ranges::iota_view{
			static_cast<std::size_t>(1), static_cast<std::size_t>(12)}) {
		dflux = calcdSpace(
					interim_us[k-1].get(),
					t, dx, max_eigenvalues,
					n_size, opts_args...);

		std::transform(
					std::execution::par_unseq,
					std::ranges::begin(interim_us[k-1].get()),
					std::ranges::end(interim_us[k-1].get()),
					std::ranges::begin(dflux),
					std::ranges::begin(interim_us[k].get()),
					[dt](const auto u, const auto df) {
			return u + 0.5 * dt * df;
		});
		updateGhostPoints(interim_us[k].get());
	}

	dflux = calcdSpace(
				interim_us[11].get(),
				t, dx, max_eigenvalues,
				n_size, opts_args...);


	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[&](std::size_t k) {
		u[k] = 1151./8505. * (interim_us[0].get())[k]
					+ 134./495. * (interim_us[1].get())[k]
					+ 142./525. * (interim_us[2].get())[k]
					+ 44./243. * (interim_us[3].get())[k]
					+ 4./45. * (interim_us[4].get())[k]
					+ 4./105. * (interim_us[5].get())[k]
					+ 4./405. * (interim_us[6].get())[k]
					+ 8./1575. * (interim_us[7].get())[k]
//					+ 0. * (interim_us[8].get())[k]
					+ 4./8505. * (interim_us[9].get())[k]
//					+ 0. * (interim_us[10].get())[k]
					+ 2./467775. * (
						(interim_us[11].get())[k] + 0.5 * dt * dflux[k]);
	});
	updateGhostPoints(u);
}

#endif // LSSPERK12_11_H
