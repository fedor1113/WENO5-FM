#ifndef EBDF5_H
#define EBDF5_H

#include <algorithm>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"
// #include "_vector4.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimestep_eBDF5(
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
		std::ranges::common_range auto& Ytemp,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_ghost_points,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* 5th Order Explicit 5 Step extrapolated
	 * Backward Differentiation Formula (eBDF5).
	 * At least it only has one stage.
	 *
	 * Stable with WENO with CFL <= 0.16.
	 */

	const std::size_t n_size = std::ranges::size(U)
			- 2 * n_ghost_points;

	auto interior_view_indices = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(n_ghost_points/*0*/))
				| std::views::take(n_size/*std::ranges::size(U)*/)
	);

	L0 = calcdSpace(U, t, dx, max_eigenvalues,
		n_size, opts_args...);

	Ytemp = U;
//	Ltemp = L0;

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(interior_view_indices),
				std::ranges::end(interior_view_indices),
				[&](std::size_t k) {
		U[k] = ((300. * (U[k] - Y1[k])
					+ 200. * Y2[k]
					- 75. * Y3[k]
					+ 12. * Y4[k])
				+ dt * 10. * (30. * L0[k]
					- 60. * L1[k]
					+ 60. * L2[k]
					- 30. * L3[k]
					+ 6. * L4[k])) * (1./137.);

//		Y4[k] = Y3[k]; L4[k] = L3[k];
//		Y3[k] = Y2[k]; L3[k] = L2[k];
//		Y2[k] = Y1[k]; L2[k] = L1[k];
//		Y1[k] = Ytemp[k]; L1[k] = L0[k];
	});

	updateGhostPoints(U);

	Y4 = Y3; L4 = L3;
	Y3 = Y2; L3 = L2;
	Y2 = Y1; L2 = L1;
	Y1 = Ytemp; L1 = L0;
}

#endif // EBDF5_H
