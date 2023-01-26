#ifndef ENO3_H
#define ENO3_H

#include <algorithm>
#include <cstddef>
// #include <iostream>
#include <ranges>
#include <span>
#include <valarray>
#include <vector>

#include "arithmeticwith.h"
#include "weno5.h"


// template <ArithmeticWith<numeric_val> T>
// T undivided_difference(T u0, T u2) {
// 	return (u2 - u0);  // ?!
// }


// T divided_difference(T x0, T x2, T u0, T u2) {
// 	return undivided_difference(u0, u2) / undivided_difference(x0, x2);
// }


// template <ArithmeticWith<numeric_val> T>
// T eno3_divided_difference(T x0, T x1, T x2, T u0, T u1, T u2) {
// 	// return ((u1 - u2) / (x2 - x1) - (u1 - u0) / (x1 - x0)) / (x2 - x0);
// 	return divided_difference(
// 		x0, x2,
// 		divided_difference(x1, x2, u1, u2),
// 		divided_difference(x0, x1, u0, u1));
// }


template <ArithmeticWith<numeric_val> T>
bool chooseENO2Stencil(const std::ranges::sized_range auto&& u_stencil) {
	/* TO DO */

	if (std::abs(u_stencil[1] - u_stencil[0])
			<= std::abs(u_stencil[2] - u_stencil[1]))
		return true;

	return false;
}


template <ArithmeticWith<numeric_val> T>
T computeENO2ReconstructionKernel(
		const std::ranges::sized_range auto&& u_stencil,
		bool use_left_substencil) {
	/* 2nd order ENO reconstructions of f(j) from all the 2 2-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 4 values [j-1, j+0, j+1, j+2, ...] for '+'
	 *               (or [j+2, j+1, j+0, j-1, ...] for '-').
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 * Choose a stencil with `which_stencil`:
	 * 	 true: [j-1, j+0          ]
	 * 	false: [          j+1, j+2]
	 */

	if (use_left_substencil)
		return (-1. * u_stencil[0] + 3. * u_stencil[1]) * 0.5;

	return (1. * u_stencil[1] + 1. * u_stencil[2]) * 0.5;
}


template <ArithmeticWith<numeric_val> T>
short chooseENO3Stencil(const std::ranges::sized_range auto&& u_stencil) {
	/* TO DO */

	if (std::abs(u_stencil[2] - u_stencil[1])
			<= std::abs(u_stencil[3] - u_stencil[2])) {
		if (std::abs(u_stencil[2] - 2 * u_stencil[1] + u_stencil[0])
				<= std::abs(u_stencil[3] - 2 * u_stencil[2] + u_stencil[1])
				) {
			return 2;
		} else {
			return 1;
		}
	}

	if (std::abs(u_stencil[3] - 2 * u_stencil[2] + u_stencil[1])
			<= std::abs(u_stencil[4] - 2 * u_stencil[3] + u_stencil[2]))
		return 1;

	return 0;
}


template <ArithmeticWith<numeric_val> T>
T computeENO3ReconstructionKernel(
		const std::ranges::sized_range auto&& u_stencil,
		short which_stencil) {
	/* 3rd order ENO reconstructions of f(j) from all the 3 3-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *               (or [j+3, j+2, j+1, j+0, j-1, ...] for '-').
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 * Choose a stencil with `which_stencil`:
	 * 	2: [j-2, j-1, j+0          ]
	 * 	1: [     j-1, j+0, j+1     ]
	 * 	0: [          j+0, j+1, j+2]
	 */
	//		{{2./6, -7./6, 11./6, 0., 0., 0.}},
	//		{{0., -1./6, 5./6, 2./6, 0., 0.}},
	//		{{0., 0., 2./6, 5./6, -1./6, 0.}}

	switch (which_stencil) {
		case 2: return (2. * u_stencil[0]
						- 7. * u_stencil[1]
						+ 11. * u_stencil[2]) / 6.; break;
		case 1: return (-1. * u_stencil[1]
						+ 5. * u_stencil[2]
						+ 2. * u_stencil[3]) / 6.; break;
		case 0: return (2. * u_stencil[2]
						+ 5. * u_stencil[3]
						- 1. * u_stencil[4]) / 6.; break;
		default: return 0.;
	}
}


template <ArithmeticWith<numeric_val> T>
T computeENO3ReconstructionKernelRev(
		const std::ranges::sized_range auto&& u_stencil,
		short which_stencil) {
	/* 3rd order ENO reconstructions of f(j) from all the 3 3-element
	 * substencils of `f_stencil` (f_minus or reversed f_plus:
	 * receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '-'
	 *               (or [j+3, j+2, j+1, j+0, j-1, ...] for '+').
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 * Choose a stencil with `which_stencil`:
	 * 	2: [j-2, j-1, j+0          ]
	 * 	1: [     j-1, j+0, j+1     ]
	 * 	0: [          j+0, j+1, j+2]
	 */
	//		{{0., 0., 0., 11./6, -7./6, 2./6}},
	//		{{0., 0., 2./6, 5./6, -1./6, 0.}},
	//		{{0., -1./6, 5./6, 2./6, 0., 0.}}

	switch (which_stencil) {
		case 2: return (-1. * u_stencil[0]
						+ 5. * u_stencil[1]
						+ 2. * u_stencil[2]) / 6.; break;
		case 1: return (2. * u_stencil[1]
						+ 5. * u_stencil[2]
						- 1. * u_stencil[3]) / 6.; break;
		case 0: return (11. * u_stencil[2]
						- 7. * u_stencil[3]
						+ 2. * u_stencil[4]) / 6.; break;
		default: return 0.;
	}
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDENO3(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_size) {
	/* Simple finite-difference essentially non-oscillatory 3-rd order
	 * (FD ENO-3) reconstruction (flux-split version).
	 */

	const std::size_t stencil_size = 5;
	const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = _actual_stencil_size / 2;  // 3

	// r = (order + 1) / 2 = 3
	const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;  // 3
	const std::size_t maxi = n_ghost_cells + n_size - 1;
	auto shifted_index_range = std::ranges::iota_view{
		mini - 1,
		maxi + 1
	};

	T uhatminus = 0.;
	T uhatplus = 0.;

	auto j_it_p = std::ranges::begin(f_plus);  // f_plus
	auto j_it_m = std::ranges::begin(f_minus);  // f_minus

	// std::vector<T> u_plus(_actual_stencil_size);
	// std::vector<T> u_minus(_actual_stencil_size);

	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size - 1);
	std::advance(j_it_m, mini - 1 + half_size + 1 - stencil_size - 1);
	auto u_plus = std::ranges::views::counted(j_it_p, 6);
	auto u_minus = std::ranges::views::counted(j_it_m, 6)
					| std::ranges::views::reverse;
	short which_stencil_pl;
	short which_stencil_mn;

	for (std::size_t j : shifted_index_range) {
		j_it_p = std::ranges::begin(f_plus);  // u_plus
		j_it_m = std::ranges::begin(f_minus);  // u_minus
		std::advance(j_it_p, j + half_size + 1 - stencil_size - 1);
		std::advance(j_it_m, j + half_size + 1 - stencil_size - 1);
		u_plus = std::ranges::views::counted(j_it_p, 6);
		u_minus = std::ranges::views::counted(j_it_m, 6)
					| std::ranges::views::reverse;

		which_stencil_pl = chooseENO3Stencil<T>(
						std::ranges::views::counted(
							std::ranges::begin(u_plus), 5)
						);

		uhatplus = computeENO3ReconstructionKernel<T>(
					std::ranges::views::counted(
						std::ranges::begin(u_plus), 5),
					which_stencil_pl
					);

		which_stencil_mn = chooseENO3Stencil<T>(
						std::ranges::subrange(
							std::ranges::begin(u_minus),
							std::ranges::end(u_minus) - 1)
						);

//		uhatminus = computeENO3ReconstructionKernelRev<T>(
//					std::ranges::views::counted(
//						std::ranges::begin(u_plus)+1, 5),
//					which_stencil
//				);
		uhatminus = computeENO3ReconstructionKernel<T>(
					std::ranges::subrange(
						std::begin(u_minus),
						std::end(u_minus) - 1
					),
					which_stencil_mn/*2 - which_stencil_mn*/
					);

		numerical_flux[j] = uhatplus + uhatminus;

//		std::cout << (u_minus_rec[j]-u_plus_rec[j])*(u[j]-u[j+1])
//					<< "\n";
	}
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFVENO3(
		const std::ranges::common_range auto&& u,
		T t,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_size) {
	/* Simple finite-volume essentially non-oscillatory 3-rd order
	 * (FV ENO-3) reconstruction.
	 */

	const std::size_t stencil_size = 5;
	const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = _actual_stencil_size / 2;  // 3

	// r = (order + 1) / 2 = 3
	const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;  // 3
	const std::size_t maxi = n_ghost_cells + n_size - 1;
	auto shifted_index_range = std::ranges::iota_view{
		mini - 1,
		maxi + 1
	};

	T uhatminus = 0.;
	T uhatplus = 0.;

	auto j_it_p = std::ranges::begin(u);  // f_plus
	// auto j_it_m = std::ranges::begin(u);  // f_minus

	// std::vector<T> u_plus(_actual_stencil_size);
	// std::vector<T> u_minus(_actual_stencil_size);

	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size - 1);
	auto u_plus = std::ranges::views::counted(j_it_p, 6);
	auto u_minus = u_plus | std::ranges::views::reverse;
	short which_stencil = chooseENO3Stencil<T>(
				std::ranges::views::counted(std::ranges::begin(u_plus), 5)
				);

	for (std::size_t j : shifted_index_range) {
		j_it_p = std::ranges::begin(u);  // u_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size - 1);
		u_plus = std::ranges::views::counted(j_it_p, 6);
		u_minus = u_plus | std::ranges::views::reverse;

		uhatplus = computeENO3ReconstructionKernel<T>(
					std::ranges::views::counted(
						std::ranges::begin(u_plus), 5),
					which_stencil
					);

		which_stencil = chooseENO3Stencil<T>(
					std::ranges::views::counted(
						std::ranges::begin(u_plus) + 1, 5)
					);

//		uhatminus = computeENO3ReconstructionKernelRev<T>(
//					std::ranges::views::counted(
//						std::ranges::begin(u_plus)+1, 5),
//					which_stencil
//				);
		uhatminus = computeENO3ReconstructionKernel<T>(
					std::ranges::subrange(
						std::begin(u_minus),
						std::end(u_minus) - 1
					),
					2 - which_stencil
					);

		u_plus_rec[j] = uhatplus;
		u_minus_rec[j] = uhatminus;

//		std::cout << (u_minus_rec[j]-u_plus_rec[j])*(u[j]-u[j+1])
//					<< "\n";
	}
}


template <ArithmeticWith<numeric_val> T, ArithmeticWith<numeric_val> VT,
			std::size_t N=5>
void calcHydroStageCharWiseFDENO3(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		const std::ranges::common_range auto&& flux,
		std::ranges::common_range auto&& numerical_flux,
		T t,
		auto&& project,
//		auto&& deproject,
		T alpha,
		std::size_t n_ghost_cells = (N + 1) / 2) {
	const unsigned order = N;
	const std::size_t stencil_size = order;
	const std::size_t _actual_stencil_size = stencil_size + 1;
	const std::size_t half_size = order / 2;

	const std::size_t r = (order + 1) / 2;
	assert(n_ghost_cells >= r);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(
				numerical_flux) - n_ghost_cells - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	auto shifted_index_range = std::ranges::common_view(
				std::views::iota(mini - 1)
					| std::views::take(maxi + 1 - (mini - 1)/* + 1*/));
	VT fhatminus = static_cast<VT>(0.);
	VT fhatplus = static_cast<VT>(0.);

	auto j_it_q = std::ranges::begin(u);
	auto j_it_f = std::ranges::begin(flux);

	std::advance(j_it_q, mini - 1 + half_size + 1 - stencil_size);
	std::advance(j_it_f, mini - 1 + half_size + 1 - stencil_size);
	auto f_stencil = std::ranges::views::counted(j_it_f,
												 _actual_stencil_size);
	auto q_stencil = std::ranges::views::counted(j_it_q,
												 _actual_stencil_size);

	std::valarray<VT> u_plus(_actual_stencil_size);
	std::valarray<VT> u_minus(_actual_stencil_size);

	auto components = {
		std::make_pair(0, &Vector4<T>::x),
		std::make_pair(1, &Vector4<T>::y),
		std::make_pair(2, &Vector4<T>::z)
//		&Vector4<T>::w
	};

	short which_stencil_pl;
	short which_stencil_mn;

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(shifted_index_range),
				std::ranges::end(shifted_index_range),
				[&](std::size_t j) {
		j_it_q = std::ranges::begin(u);
		j_it_f = std::ranges::begin(flux);
		std::advance(j_it_q, j + half_size + 1 - stencil_size);
		q_stencil = std::ranges::views::counted(j_it_q,
												_actual_stencil_size);

		std::advance(j_it_f, j + half_size + 1 - stencil_size);
		f_stencil = std::ranges::views::counted(j_it_f,
												_actual_stencil_size);

		auto proj_u_j = [j, &q_avg, &project](auto u) -> decltype(u) {
			const auto q = q_avg[j];
			return project(q, u);
		};

		for (std::size_t k = 0; k < _actual_stencil_size; ++ k)
			u_plus[k] = (proj_u_j(f_stencil[k])
					+ /*1.1 * */alpha * proj_u_j(q_stencil[k])) * 0.5;

		for (std::size_t k = 0; k < _actual_stencil_size; ++ k)
			u_minus[k] = (proj_u_j(f_stencil[_actual_stencil_size - k - 1])
					- /*1.1 * */alpha * proj_u_j(
						q_stencil[_actual_stencil_size - k - 1])) * 0.5;

		for (auto& comp : components) {
			which_stencil_pl = chooseENO3Stencil<T>(
							std::ranges::views::counted(
								std::ranges::begin(u_plus), order)
						| std::ranges::views::transform(
							comp.second)
							);
			fhatplus[comp.first] = computeENO3ReconstructionKernel<T>(
				std::ranges::views::counted(
							std::ranges::begin(u_plus), order)
						| std::ranges::views::transform(
							comp.second), which_stencil_pl
			);

			which_stencil_mn = chooseENO3Stencil<T>(
							std::ranges::views::counted(
								std::ranges::begin(u_minus), order)
						| std::ranges::views::transform(
							comp.second)
							);
			fhatminus[comp.first] = computeENO3ReconstructionKernel<T>(
				std::ranges::views::counted(
							std::ranges::begin(u_minus), order)
						| std::ranges::views::transform(
							comp.second), which_stencil_mn
			);
		}

		numerical_flux[j] = fhatplus + fhatminus;
	});
}


template <ArithmeticWith<numeric_val> T, ArithmeticWith<numeric_val> VT>
void calcHydroStageCharWiseFDENO2(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		const std::ranges::common_range auto&& flux,
		std::ranges::common_range auto&& numerical_flux,
		T t,
		auto&& project,
//		auto&& deproject,
		T alpha,
		std::size_t n_ghost_cells = 2) {
	const unsigned order = 2;
	const std::size_t stencil_size = 3;
	const std::size_t _actual_stencil_size = stencil_size + 1;
	const std::size_t half_size = order / 2;  // 1

	const std::size_t r = (order + 1) / 2;  // 1
	assert(n_ghost_cells > r);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(
				numerical_flux) - n_ghost_cells - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	auto shifted_index_range = std::ranges::common_view(
				std::views::iota(mini - 1)
					| std::views::take(maxi + 1 - (mini - 1)/* + 1*/));
	VT fhatminus = static_cast<VT>(0.);
	VT fhatplus = static_cast<VT>(0.);

	auto j_it_q = std::ranges::begin(u);
	auto j_it_f = std::ranges::begin(flux);

	std::advance(j_it_q, mini - 1 + half_size + 1 - stencil_size);
	std::advance(j_it_f, mini - 1 + half_size + 1 - stencil_size);
	auto f_stencil = std::ranges::views::counted(j_it_f,
												 _actual_stencil_size);
	auto q_stencil = std::ranges::views::counted(j_it_q,
												 _actual_stencil_size);

	std::valarray<VT> u_plus(_actual_stencil_size);
	std::valarray<VT> u_minus(_actual_stencil_size);

	auto components = {
		std::make_pair(0, &Vector4<T>::x),
		std::make_pair(1, &Vector4<T>::y),
		std::make_pair(2, &Vector4<T>::z)
//		&Vector4<T>::w
	};

	bool which_stencil_pl;
	bool which_stencil_mn;

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(shifted_index_range),
				std::ranges::end(shifted_index_range),
				[&](std::size_t j) {
		j_it_q = std::ranges::begin(u);
		j_it_f = std::ranges::begin(flux);
		std::advance(j_it_q, j + half_size + 1 - stencil_size);
		q_stencil = std::ranges::views::counted(j_it_q,
												_actual_stencil_size);

		std::advance(j_it_f, j + half_size + 1 - stencil_size);
		f_stencil = std::ranges::views::counted(j_it_f,
												_actual_stencil_size);

		auto proj_u_j = [j, &q_avg, &project](auto u) -> decltype(u) {
			const auto q = q_avg[j];
			return project(q, u);
		};

		for (std::size_t k = 0; k < _actual_stencil_size; ++ k)
			u_plus[k] = (proj_u_j(f_stencil[k])
					+ /*1.1 * */alpha * proj_u_j(q_stencil[k])) * 0.5;

		for (std::size_t k = 0; k < _actual_stencil_size; ++ k)
			u_minus[k] = (proj_u_j(f_stencil[_actual_stencil_size - k - 1])
					- /*1.1 * */alpha * proj_u_j(
						q_stencil[_actual_stencil_size - k - 1])) * 0.5;

		for (auto& comp : components) {
			which_stencil_pl = chooseENO2Stencil<T>(
							std::ranges::views::counted(
								std::ranges::begin(u_plus), stencil_size)
						| std::ranges::views::transform(
							comp.second)
							);
			fhatplus[comp.first] = computeENO2ReconstructionKernel<T>(
				std::ranges::views::counted(
							std::ranges::begin(u_plus), stencil_size)
						| std::ranges::views::transform(
							comp.second), which_stencil_pl
			);

			which_stencil_mn = chooseENO2Stencil<T>(
							std::ranges::views::counted(
								std::ranges::begin(u_minus), stencil_size)
						| std::ranges::views::transform(
							comp.second)
							);
			fhatminus[comp.first] = computeENO2ReconstructionKernel<T>(
				std::ranges::views::counted(
							std::ranges::begin(u_minus), stencil_size)
						| std::ranges::views::transform(
							comp.second), /*1 - */which_stencil_mn
			);
		}

		numerical_flux[j] = fhatplus + fhatminus;
	});
}


template <ArithmeticWith<numeric_val> T, ArithmeticWith<numeric_val> VT,
			std::size_t N=5>
void calcHydroStageCharWiseFDMPENO3(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		const std::ranges::common_range auto&& flux,
		std::ranges::common_range auto&& numerical_flux,
		T t,
		auto&& project,
//		auto&& deproject,
		T alpha,
		std::size_t n_ghost_cells = (N + 1) / 2) {
	const unsigned order = N;
	const std::size_t stencil_size = order;
	const std::size_t _actual_stencil_size = stencil_size + 1;
	const std::size_t half_size = order / 2;

	const std::size_t r = (order + 1) / 2;
	assert(n_ghost_cells >= r);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(
				numerical_flux) - n_ghost_cells - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	auto shifted_index_range = std::ranges::common_view(
				std::views::iota(mini - 1)
					| std::views::take(maxi + 1 - (mini - 1)/* + 1*/));
	VT fhatminus = static_cast<VT>(0.);
	VT fhatplus = static_cast<VT>(0.);

	auto j_it_q = std::ranges::begin(u);
	auto j_it_f = std::ranges::begin(flux);

	std::advance(j_it_q, mini - 1 + half_size + 1 - stencil_size);
	std::advance(j_it_f, mini - 1 + half_size + 1 - stencil_size);
	auto f_stencil = std::ranges::views::counted(j_it_f,
												 _actual_stencil_size);
	auto q_stencil = std::ranges::views::counted(j_it_q,
												 _actual_stencil_size);

	std::valarray<VT> u_plus(_actual_stencil_size);
	std::valarray<VT> u_minus(_actual_stencil_size);

	auto components = {
		std::make_pair(0, &Vector4<T>::x),
		std::make_pair(1, &Vector4<T>::y),
		std::make_pair(2, &Vector4<T>::z)
//		&Vector4<T>::w
	};

	short which_stencil_pl;
	short which_stencil_mn;

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(shifted_index_range),
				std::ranges::end(shifted_index_range),
				[&](std::size_t j) {
		j_it_q = std::ranges::begin(u);
		j_it_f = std::ranges::begin(flux);
		std::advance(j_it_q, j + half_size + 1 - stencil_size);
		q_stencil = std::ranges::views::counted(j_it_q,
												_actual_stencil_size);

		std::advance(j_it_f, j + half_size + 1 - stencil_size);
		f_stencil = std::ranges::views::counted(j_it_f,
												_actual_stencil_size);

		auto proj_u_j = [j, &q_avg, &project](auto u) -> decltype(u) {
			const auto q = q_avg[j];
			return project(q, u);
		};

		for (std::size_t k = 0; k < _actual_stencil_size; ++ k)
			u_plus[k] = (proj_u_j(f_stencil[k])
					+ /*1.1 * */alpha * proj_u_j(q_stencil[k])) * 0.5;

		for (std::size_t k = 0; k < _actual_stencil_size; ++ k)
			u_minus[k] = (proj_u_j(f_stencil[_actual_stencil_size - k - 1])
					- /*1.1 * */alpha * proj_u_j(
						q_stencil[_actual_stencil_size - k - 1])) * 0.5;

		for (auto& comp : components) {
			which_stencil_pl = chooseENO3Stencil<T>(
							std::ranges::views::counted(
								std::ranges::begin(u_plus), order)
						| std::ranges::views::transform(
							comp.second)
							);
			fhatplus[comp.first] = computeENO3ReconstructionKernel<T>(
				std::ranges::views::counted(
							std::ranges::begin(u_plus), order)
						| std::ranges::views::transform(
							comp.second), which_stencil_pl
			);
			fhatplus[comp.first] = MPLimiterMM(
						u_plus[stencil_size / 2 - 1][comp.first],
						u_plus[stencil_size / 2 + 0][comp.first],
						u_plus[stencil_size / 2 + 1][comp.first],
						u_plus[stencil_size / 2 + 2][comp.first],
						fhatplus[comp.first]);

			which_stencil_mn = chooseENO3Stencil<T>(
							std::ranges::views::counted(
								std::ranges::begin(u_minus), order)
						| std::ranges::views::transform(
							comp.second)
							);
			fhatminus[comp.first] = computeENO3ReconstructionKernel<T>(
				std::ranges::views::counted(
							std::ranges::begin(u_minus), order)
						| std::ranges::views::transform(
							comp.second), which_stencil_mn
			);
			fhatminus[comp.first] = MPLimiterMM(
						u_minus[stencil_size / 2 - 1][comp.first],
						u_minus[stencil_size / 2 + 0][comp.first],
						u_minus[stencil_size / 2 + 1][comp.first],
						u_minus[stencil_size / 2 + 2][comp.first],
						fhatminus[comp.first]);
		}

		numerical_flux[j] = fhatplus + fhatminus;
	});
}

#endif // ENO3_H
