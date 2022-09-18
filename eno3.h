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
// #include "weno5.h"


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
		j_it_p = std::ranges::begin(u);  // u_plus
		j_it_m = std::ranges::begin(u);  // u_minus
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
					2 - which_stencil_mn
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

#endif // ENO3_H
