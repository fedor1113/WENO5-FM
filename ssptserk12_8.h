#ifndef SSPTSERK12_8_H
#define SSPTSERK12_8_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimestepSSPTSERK12_8(
		std::ranges::common_range auto& u_curr,
		std::ranges::common_range auto& u_prev,
		std::ranges::common_range auto& dflux,
		std::array<
			std::reference_wrapper<std::valarray<T>
		>, 13> interim_us,
		std::array<
			std::reference_wrapper<std::valarray<T>
		>, 13> interim_fluxes,
		numeric_val t, numeric_val dt, numeric_val dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_ghost_points,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* 8th Order 2 Step 12 Stage Explicit
	 * Strong Stability Preserving
	 * Runge-Kutta Scheme (SSPTSERK(12,8))
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 *
	 * Inefficient implementation! Esp. mem.
	 */

	const std::size_t n_size = std::ranges::size(u_curr)
			- 2 * n_ghost_points;

	numeric_val theta = 4.796147528566197e-05;
	numeric_val one_over_C = 1./0.9416;
	constexpr const static std::array<const numeric_val, 13> d {
		1.000000000000000, 0.000000000000000, 0.036513886685777,
		0.000000000000000, 0.004205435886220, 0.000457751617285,
		0.000000000000000, 0.007407526543898, 0.000486094553850,
		0.000000000000000, 0.000000000000000, 0.000000000000000,
		0.000000000000000
	};
	constexpr const static std::array<const numeric_val, 13> eta {
		0.000000000000000, 0.033190060418244, 0.001567085177702,
		0.014033053074861, 0.017979737866822, 0.094582502432986,
		0.082918042281378, 0.020622633348484, 0.033521998905243,
		0.092066893962539, 0.076089630105122, 0.070505470986376,
		0.072975312278165
	};
	constexpr const static std::array<std::array<const numeric_val, 12>, 13> q {{
		{{
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.017683145596548, 0.154785324942633, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.001154189099465, 0.000000000000000, 0.200161251441789,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.113729301017461, 0.000000000000000,
			0.057780552515458, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.061188134340758, 0.000000000000000,
			0.000000000000000, 0.165254103192244, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000065395819685, 0.068824803789446, 0.008642531617482,
			0.000000000000000, 0.000000000000000, 0.229847794524568,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.133098034326412, 0.000000000000000,
			0.000000000000000, 0.005039627904425, 0.000000000000000,
			0.252990567222936, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.080582670156691, 0.000000000000000,
			0.000000000000000, 0.069726774932478, 0.000000000000000,
			0.000000000000000, 0.324486261336648, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000042696255773, 0.038242841051944, 0.000000000000000,
			0.029907847389714, 0.022904196667572, 0.095367316002296,
			0.176462398918299, 0.000000000000000, 0.120659479468128,
			0.000000000000000, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.071728403470890, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.281349762794588, 0.000000000000000, 0.000000000000000,
			0.166819833904944, 0.000000000000000, 0.000000000000000
		}},
		{{
			0.000116117869841, 0.053869626312442, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.327578464731509, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.157699899495506, 0.000000000000000
		}},
		{{
			0.000019430720566, 0.009079504342639, 0.000000000000000,
			0.000000000000000, 0.130730221736770, 0.000000000000000,
			0.149446805276484, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.314802533082027
		}}
	}};  // 15 digits after the point

	interim_fluxes[0].get() = calcdSpace(u_prev, t, dx, max_eigenvalues,
				n_size, opts_args...);
	dflux = calcdSpace(u_curr, t, dx, max_eigenvalues,
				n_size, opts_args...);
	interim_fluxes[1].get() = dflux;
	interim_us[0].get() = u_prev;
	interim_us[1].get() = u_curr;

	auto iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(std::ranges::size(u_curr))
	);
	for (std::size_t k = 1; k <= 12; ++ k) {
		std::for_each(
//					std::execution::par_unseq,
					std::ranges::begin(iv),
					std::ranges::end(iv),
					[&](std::size_t pt_idx) {
			(interim_us[k].get())[pt_idx] = d[k] * u_prev[pt_idx]
						+ (1. - d[k]) * u_curr[pt_idx];

			for (std::size_t j = 0; j < k; ++ j)
				(interim_us[k].get())[pt_idx] += q[k][j]
							* (-u_curr[pt_idx]
							   + (interim_us[j].get())[pt_idx]
							   + dt * one_over_C
									* (interim_fluxes[j].get())[pt_idx]);
		});

		updateGhostPoints(interim_us[k].get());

		(interim_fluxes[k].get()) = calcdSpace(
					(interim_us[k].get()),
					t, dx, max_eigenvalues,
					n_size, opts_args...);

//		updateGhostPoints(interim_us[k].get());
	}

	std::for_each(
//				std::execution::par_unseq,
				std::ranges::begin(iv),
				std::ranges::end(iv),
				[&](std::size_t pt_idx) {
		T temp;
		temp = theta * u_prev[pt_idx]
					+ (1. - theta) * u_curr[pt_idx];

		for (std::size_t j = 0; j <= 12; ++ j)
			temp += eta[j]
						* (-u_curr[pt_idx]
						   + (interim_us[j].get())[pt_idx]
						   + dt * one_over_C
								* (interim_fluxes[j].get())[pt_idx]);

		u_prev[pt_idx] = u_curr[pt_idx];
		u_curr[pt_idx] = temp;
	});

	updateGhostPoints(u_curr);
}

#endif // SSPTSERK12_8_H
