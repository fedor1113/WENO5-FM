#ifndef SSPTSERK8_5_H
#define SSPTSERK8_5_H

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
		>, 9> interim_us,
		std::array<
			std::reference_wrapper<std::valarray<T>
		>, 9> interim_fluxes,
		T t, T dt, T dx,
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
	 * Only SSP for linear problems!!!
	 */

	const std::size_t n_size = std::ranges::size(u_curr)
			- 2 * n_ghost_points;

	T theta = 0.;
	T one_over_C = 1./ 3.5794;
	constexpr const static std::array<const T, 9> d {
		1.000000000000000, 0.000000000000000, 0.000000000000000,
		0.000000000000000, 0.000000000000000, 0.000000000000000,
		0.000000000000000, 0.003674184820260, 0.000000000000000
	};
	constexpr const static std::array<const T, 9> eta {
		0.000000000000000, 0.000000000000000, 0.179502832154858,
		0.073789956884809, 0.000000000000000, 0.000000000000000,
		0.017607159013167, 0.000000000000000, 0.729100051947166
	};
	constexpr const static std::array<std::array<const T, 8>, 9> q {{
		{{
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000
		}},
		{{
			0.085330772947643, 0.914669227052357, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000
		}},
		{{
			0.058121281984411, 0.000000000000000, 0.941878718015589,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.036365639242841, 0.000000000000000,
			0.802870131352638, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.491214340660555, 0.000000000000000,
			0.000000000000000, 0.508785659339445, 0.000000000000000,
			0.000000000000000, 0.000000000000000
		}},
		{{
			0.000000000000000, 0.566135231631241, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.433864768368758,
			0.000000000000000, 0.000000000000000
		}},
		{{
			0.020705281786630, 0.091646079651566, 0.000000000000000,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.883974453741544, 0.000000000000000
		}},
		{{
			0.008506650138784, 0.110261531523242, 0.030113037742445,
			0.000000000000000, 0.000000000000000, 0.000000000000000,
			0.000000000000000, 0.851118780595529
		}},
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
	for (std::size_t k = 1; k <= 8; ++ k) {
		std::for_each(
					std::execution::par_unseq,
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

		updateGhostPoints(interim_us[k].get());
	}

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv),
				std::ranges::end(iv),
				[&](std::size_t pt_idx) {
		T temp;
		temp = theta * u_prev[pt_idx]
					+ (1. - theta) * u_curr[pt_idx];

		for (std::size_t j = 0; j <= 8; ++ j)
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

#endif // SSPTSERK8_5_H
