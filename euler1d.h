#ifndef EULER1D_H
#define EULER1D_H

#include <algorithm>
#include <execution>
#include <ranges>

#include "Eigen/Dense"

#include "_vector4.h"
#include "arithmeticwith.h"
#include "eno3.h"
// #include "eos.h"
// #include "exactsolver.h"
// #include "hllc_solver.h"
#include "lf_flux.h"
// #include "miegruneisen.h"
#include "roe.h"
// #include "rk6_5.h"
#include "ssprk33.h"
// #include "ssprk10_4.h"
#include "ssptserk12_8.h"
#include "ebdf5.h"
// #include "eulerforward.h"
#include "weno5.h"

//template <ArithmeticWith<numeric_val> T>
//T average(T left, T right) {
//	/* A simple arithmetic mean average of 2 values. */
//	return (left + right) * 0.5;
//}


template <ArithmeticWith<numeric_val> T,
			ArithmeticWith<numeric_val> T2=numeric_val>
T average(T left, T right) {
	/* A simple arithmetic mean average of 2 values. */

	T2 sqrt_left_rho = std::sqrt(std::abs(left[0]));
	T2 sqrt_right_rho = std::sqrt(std::abs(right[0]));
	T2 avg_rho = sqrt_left_rho * sqrt_right_rho;
	T2 sqrt_rho_sum = sqrt_left_rho + sqrt_right_rho;

	T2 left_v = left[1] / left[0];
	T2 right_v = right[1] / right[0];
	T2 avg_rho_v = avg_rho
			* (sqrt_left_rho * left_v + sqrt_right_rho * right_v)
				/ sqrt_rho_sum;

	T2 left_E = left[2] / left[0];
	T2 right_E = right[2] / right[0];
	T2 avg_E = (sqrt_left_rho * left_E + sqrt_right_rho * right_E)
				/ sqrt_rho_sum;
	T2 avg_rho_E = avg_rho * avg_E;

	return {avg_rho, avg_rho_v, avg_rho_E, 0.};
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcPhysFlux(
		const std::valarray<Vector4<T>>& u_arr/*, T gamma = 1.4*/) {
	/* Calculate physical fluxes of conservative variables
	 * in all points in the computation domain u_arr.
	 */

	std::valarray<Vector4<T>> res(std::ranges::size(u_arr));

	std::transform(std::ranges::begin(u_arr),
				   std::ranges::end(u_arr),
				   std::ranges::begin(res),
				   [/*&gamma*/](Vector4<T> v) {
		return calcPhysicalFluxFromConservativeVec<T>(v/*, gamma*/);
	});

	return res;  // f_arr
}


template <ArithmeticWith<numeric_val> T>
T calcMaxWaveSpeedD(
		const std::ranges::common_range auto& u_arr/*,
		T gamma = 1.4*/) {
	/* Calculate |df/du| for 1D Euler eq'ns. */

	return std::ranges::max(std::ranges::transform_view(
		std::as_const(u_arr),
		[/*gamma*/](const auto& u_arr_vec_pt) -> T {
			return (std::sqrt(calcSquareSoundSpeed(
								u_arr_vec_pt[0],
								u_arr_vec_pt[1],
								u_arr_vec_pt[2]/*, gamma*/))
					+ std::abs(u_arr_vec_pt[1] / u_arr_vec_pt[0]));
		}
	));
}


template <ArithmeticWith<numeric_val> T>
auto calcExactEulerFlux = [](auto& u_plus, auto& u_minus, auto gamma) {
	std::valarray<Vector4<T>> res(std::ranges::size(u_plus));
	std::transform(
				std::ranges::begin(u_plus), std::ranges::end(u_plus),
				std::ranges::begin(u_minus),
				std::ranges::begin(res),
				[gamma](const auto u_pl, const auto u_mn) {
					return calcExactFlux(
							u_pl[0], u_pl[1], u_pl[2],
							u_mn[0], u_mn[1], u_mn[2], gamma);
				});

	return res;
};


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxComponentWiseFDWENO5(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		std::size_t number_of_ghost_points = 3,
		T gamma = 1.4,
		T eps = 1e-40,
		T p = 2.) {
	std::valarray<Vector4<T>> res = calcPhysFlux(u/*, gamma*/);
//	std::array<std::valarray<T>, 2> monotone_flux_components = {
//		std::valarray<T>(std::ranges::size(u)),
//		std::valarray<T>(std::ranges::size(u)),
//	};
//	SplitFluxes monotone_flux_components {
//		std::valarray<T>(std::ranges::size(u)),
//		std::valarray<T>(std::ranges::size(u))
//	};

	// auto iv = std::ranges::iota_view{0, 4};
	auto components = {
		&Vector4<T>::x,
		&Vector4<T>::y,
		&Vector4<T>::z,
		&Vector4<T>::w
	};
	std::for_each(
//			std::execution::par_unseq,
			std::ranges::begin(components),
			std::ranges::end(components),
			[&](auto kth_vector_component) {
		auto [f_plus, f_minus] = splitFluxAsLaxFriedrichs(
					u | std::ranges::views::transform(
						kth_vector_component),
					res | std::ranges::views::transform(
						kth_vector_component), lam[0]
					);

		calcHydroStageFDWENO5FM<T>(
			std::ranges::views::all(f_plus),
			std::ranges::views::all(f_minus), t,
			res | std::ranges::views::transform(kth_vector_component),
			number_of_ghost_points, eps, p
			);
//		calcHydroStageFDWENO7FM<T>(
//			std::ranges::views::all(f_plus),
//			std::ranges::views::all(f_minus), t,
//			res | std::ranges::views::transform(kth_vector_component),
//			number_of_ghost_points, eps, p
//			);
	});

	return res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxCharacteristicWiseFDWENO5(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		std::size_t ghost_point_number = 3,
		T gamma = 1.4,
		T eps = 1e-40,
		T p = 2.) {
//	const std::size_t n_size = std::ranges::size(u)
//			- 2 * ghost_point_number;

	std::valarray<Vector4<T>> avg(std::ranges::size(u));
	std::valarray<Vector4<T>> proj_u(std::ranges::size(u));
	std::valarray<Vector4<T>> proj_f(std::ranges::size(u));
	std::valarray<Vector4<T>> res(std::ranges::size(u));
	std::valarray<Vector4<T>> flux = calcPhysFlux(u/*, gamma*/);

//	auto interior_view = std::ranges::views::drop(ghost_point_number)
//			| std::ranges::views::take(n_size)
//			| std::ranges::views::common;
//	auto interior_view_shifted_by_neg_1
//			= std::ranges::views::drop(ghost_point_number - 1)
//				| std::ranges::views::take(n_size)
//				| std::ranges::views::common;

//	auto u_view = std::ranges::views::all(u);
//	auto avg_view = std::ranges::views::all(avg);
//	auto u_proj_view = std::ranges::views::all(proj_u);
//	auto f_proj_view = std::ranges::views::all(proj_f);
//	auto res_view = std::ranges::views::all(res);

//	std::transform(
//				std::execution::par_unseq,
//				std::ranges::begin(u_view) | interior_view,
//				std::ranges::end(u_view) | interior_view,
//				std::ranges::begin(u_view) | interior_view_shifted_by_neg_1,
//				std::ranges::begin(avg_view) | interior_view,
//				[](auto q_l, auto q_r) {
//		return average<Vector4<T>>(q_l, q_r);
//	});
	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(u) + 1,
				std::ranges::end(u),
				std::ranges::begin(u),
				std::ranges::begin(avg) + 1,
				[](auto q_l, auto q_r) {
		return average<Vector4<T>>(q_l, q_r);
	});
	// updateGhostPointsTransmissive(avg);
	// updateGhostPointsPeriodic(avg);

	auto project = [gamma](
			Vector4<T> q_ast, Vector4<T> vec) -> Vector4<T> {
		return projectOntoCharacteristics<T>(q_ast, vec/*, gamma*/);
	};

//	auto deproject = [gamma](
//			Vector4<T> q_ast, Vector4<T> vec) -> Vector4<T> {
//		return projectCharacteristicVariablesBackOntoConserved<T>(q_ast, vec, gamma);
//	};

//	calcHydroStageCharWiseFDWENO5FM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, /*deproject,*/ lam[0],
//				ghost_point_number, eps, p);
	calcHydroStageCharWiseFDMPWENO5FM<T, Vector4<T>>(
				std::ranges::views::all(u),
				std::ranges::views::all(avg),
				std::ranges::views::all(flux),
				res, t,
				project, /*deproject,*/ lam[0],
				ghost_point_number, eps, p);

//	calcHydroStageCharWiseFDWENO7FM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, lam[0],
//				ghost_point_number, eps, p);

//	std::transform(
//				std::execution::par_unseq,
//				std::ranges::begin(avg_view) | interior_view,
//				std::ranges::end(avg_view) | interior_view,
//				std::ranges::begin(res_view) | interior_view,
//				std::ranges::begin(res_view) | interior_view,
//				[gamma](auto q_ast, auto f) {
//		return projectCharacteristicVariablesBackOntoConserved<Vector4<T>>(
//					q_ast, f, gamma);
//	});

	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(avg),
				std::ranges::end(avg),
				std::ranges::begin(res),
				std::ranges::begin(res),
				[gamma](auto q_ast, auto f) {
		return projectCharacteristicVariablesBackOntoConserved<T>(
					q_ast, f/*, gamma*/);
	});
	// updateGhostPointsPeriodic(res);

	return res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxCharacteristicWiseFDWENO7(
		const std::ranges::common_range auto& u,
		/*T dx, T dt, */T t, const std::ranges::common_range auto& lam,
		std::size_t ghost_point_number = 4,
		T gamma = 1.4,
		T eps = 1e-40,
		T p = 2.) {
//	const std::size_t n_size = std::ranges::size(u)
//			- 2 * ghost_point_number;

	std::valarray<Vector4<T>> avg(std::ranges::size(u));
	std::valarray<Vector4<T>> proj_u(std::ranges::size(u));
	std::valarray<Vector4<T>> proj_f(std::ranges::size(u));
	std::valarray<Vector4<T>> res(std::ranges::size(u));
	std::valarray<Vector4<T>> flux = calcPhysFlux(u/*, gamma*/);

	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(u) + 1,
				std::ranges::end(u),
				std::ranges::begin(u),
				std::ranges::begin(avg) + 1,
				[](auto q_l, auto q_r) {
		return average<Vector4<T>>(q_l, q_r);
	});
	// updateGhostPointsTransmissive(avg);
	// updateGhostPointsPeriodic(avg);

	auto project = [gamma](
			Vector4<T> q_ast, Vector4<T> vec) -> Vector4<T> {
		return projectOntoCharacteristics<T>(q_ast, vec/*, gamma*/);
	};

//	calcHydroStageCharWiseFDWENO7FM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, lam[0],
//				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDWENO7SM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, lam[0],
//				ghost_point_number, eps, p);
//calcHydroStageCharWiseFDMPWENO7SM<T, Vector4<T>>(
//			std::ranges::views::all(u),
//			std::ranges::views::all(avg),
//			std::ranges::views::all(flux),
//			res, t,
//			project, lam[0],
//			ghost_point_number, eps, p);
	calcHydroStageCharWiseFDMPWENO7SM<T, Vector4<T>>(
				std::ranges::views::all(u),
				std::ranges::views::all(avg),
				std::ranges::views::all(flux),
				res, t,
				project, lam[0],
				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDPORWENO7SM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, dx, dt,
//				project, lam[0],
//				ghost_point_number, eps, p);

	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(avg),
				std::ranges::end(avg),
				std::ranges::begin(res),
				std::ranges::begin(res),
				[gamma](auto q_ast, auto f) {
		return projectCharacteristicVariablesBackOntoConserved<T>(
					q_ast, f/*, gamma*/);
	});
	// updateGhostPointsPeriodic(res);

	return res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxCharacteristicWiseFDWENO9(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		std::size_t ghost_point_number = 5,
		T gamma = 1.4,
		T eps = 1e-100,
		T p = 2.) {
//	const std::size_t n_size = std::ranges::size(u)
//			- 2 * ghost_point_number;

	std::valarray<Vector4<T>> avg(std::ranges::size(u));
	std::valarray<Vector4<T>> proj_u(std::ranges::size(u));
	std::valarray<Vector4<T>> proj_f(std::ranges::size(u));
	std::valarray<Vector4<T>> res(std::ranges::size(u));
	std::valarray<Vector4<T>> flux = calcPhysFlux(u/*, gamma*/);

	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(u) + 1,
				std::ranges::end(u),
				std::ranges::begin(u),
				std::ranges::begin(avg) + 1,
				[](auto q_l, auto q_r) {
		return average<Vector4<T>>(q_l, q_r);
	});
	// updateGhostPointsTransmissive(avg);
	// updateGhostPointsPeriodic(avg);

	auto project = [gamma](
			Vector4<T> q_ast, Vector4<T> vec) -> Vector4<T> {
		return projectOntoCharacteristics<T>(q_ast, vec/*, gamma*/);
	};
//	auto project = [gamma](
//			Vector4<T> left, Vector4<T> right, Vector4<T> vec) -> Vector4<T> {
//		return projectOntoCharacteristicsRoe<T>(left, right, vec);
//	};

//	calcHydroStageCharWiseFDMPWENO5FM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, /*deproject,*/ lam[0],
//				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDWENO5FIM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, /*deproject,*/ lam[0],
//				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDMPWENO5FIM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, /*deproject,*/ lam[0],
//				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDWENO9FM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, lam[0],
//				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDMPWENO7S<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, lam[0],
//				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDMPWENO9S<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, lam[0],
//				ghost_point_number, eps, p);
	calcHydroStageCharWiseFDMPWENO9SM<T, Vector4<T>>(
				std::ranges::views::all(u),
				std::ranges::views::all(avg),
				std::ranges::views::all(flux),
				res, t,
				project, lam[0],
				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDMPWENO11SM<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, lam[0],
//				ghost_point_number, eps, p);
//	calcHydroStageCharWiseFDMPWENO11S<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, lam[0],
//				ghost_point_number, eps, p);

	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(avg),
				std::ranges::end(avg),
				std::ranges::begin(res),
				std::ranges::begin(res),
				[gamma](auto q_ast, auto f) {
		return projectCharacteristicVariablesBackOntoConserved<T>(
					q_ast, f/*, gamma*/);
	});
//	for (std::size_t k = 1; k < std::ranges::size(u) - 1; ++ k) {
//		res[k] = projectCharacteristicVariablesBackOntoConservedRoe(
//					u[k], u[k + 1], res[k]);
//	}
	// updateGhostPointsPeriodic(res);

	return res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxComponentWiseFV(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		std::size_t n_size,
		auto&& calcStage,
		auto&& calcFlux) {
	std::valarray<Vector4<T>> u_plus(Vector4<T>::ZERO, u.size());
	std::valarray<Vector4<T>> u_minus(Vector4<T>::ZERO, u.size());

	// auto iv = std::ranges::iota_view{0, 4};
	auto components = {
		&Vector4<T>::x,
		&Vector4<T>::y,
		&Vector4<T>::z,
		&Vector4<T>::w
	};
	std::for_each(
			std::execution::par_unseq,
			std::ranges::begin(components),
			std::ranges::end(components),
			[&](auto kth_vector_component) {
		calcStage(
			u | std::ranges::views::transform(kth_vector_component),
			t,
			u_plus
				| std::ranges::views::transform(kth_vector_component),
			u_minus
				| std::ranges::views::transform(kth_vector_component),
			n_size
			);
	});

	std::valarray<Vector4<T>> res = calcFlux(u_plus, u_minus);

	return res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxComponentWiseFVWENO5(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		std::size_t n_ghost_points,
		T gamma = 1.4,
		T eps = 1e-40,
		T p = 2.) {

	const std::size_t n_size = std::ranges::size(u)
			- 2 * n_ghost_points;

//	auto calcExactFluxLambda = [gamma](auto& u_plus, auto& u_minus) {
//		return calcExactEulerFlux<T>(u_plus, u_minus, gamma);
//	};

	auto calcLFLambda = [gamma, &lam](auto& u_plus, auto& u_minus) {
		std::valarray<Vector4<T>> res(std::ranges::size(u_plus));
		calcLaxFriedrichsNumericalFlux(u_plus, u_minus, res,
			[gamma, &lam](const Vector4<T>& u) {
				return calcPhysicalFluxFromConservativeVec<T>(u/*, gamma*/);
			},
			lam[0]);

		return res;
	};

//	auto calcRoeLambda = [gamma](auto& u_plus, auto& u_minus) {
//		std::valarray<Vector4<T>> res(std::ranges::size(u_plus));
//		calcRoeFlux(u_plus, u_minus, res, gamma);

//		return res;
//	};

//	auto calcHLLCLambda = [gamma, &lam](auto& u_plus, auto& u_minus) {
//		std::valarray<Vector4<T>> res(std::ranges::size(u_plus));
//		calcHLLCFlux(u_plus, u_minus, res, gamma);

//		return res;
//	};

	return calcFluxComponentWiseFV(
				u, t, lam, n_size,
				[eps, p, n_ghost_points](const auto&& u,
						auto t,
						auto&& u_plus_rec,
						auto&& u_minus_rec,
						auto n_size) {
					calcHydroStageFVWENO5FM<T>(
								std::forward<decltype(u)>(u),
								t,
								std::forward<decltype(u_plus_rec)>(
									u_plus_rec),
								std::forward<decltype(u_minus_rec)>(
									u_minus_rec),
								n_ghost_points,
								eps,
								p);
				},
				/*calcRoeLambda*/
				calcLFLambda
				/*calcHLLCLambda*/
				/*calcExactFluxLambda*/);
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxComponentWiseFVENO3(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		std::size_t n_size, T gamma = 1.4) {

//	auto calcExactFluxLambda = [gamma](auto& u_plus, auto& u_minus) {
//		return calcExactEulerFlux<T>(u_plus, u_minus, gamma);
//	};

	auto calcLFLambda = [gamma, &lam](auto& u_plus, auto& u_minus) {
		std::valarray<Vector4<T>> res(std::ranges::size(u_plus));
		calcLaxFriedrichsNumericalFlux(u_plus, u_minus, res,
			[gamma, &lam](const Vector4<T> u) {
				return calcPhysicalFluxFromConservativeVec<T>(u/*, gamma*/);
			},
			lam[0]);

		return res;
	};

	return calcFluxComponentWiseFV(
				u, t, lam, n_size,
				[](const auto&& u,
						auto t,
						auto&& u_plus_rec,
						auto&& u_minus_rec,
						auto n_size) {
					calcHydroStageFVENO3<T>(
								std::forward<decltype(u)>(u),
								t,
								std::forward<decltype(u_plus_rec)>(
									u_plus_rec),
								std::forward<decltype(u_minus_rec)>(
									u_minus_rec),
								n_size);
				},
				calcLFLambda/*calcExactFluxLambda*/);
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxCharacteristicWiseFDENO3(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		std::size_t ghost_point_number = 3,
		T gamma = 1.4,
		T eps = 1e-40,
		T p = 2.) {
//	const std::size_t n_size = std::ranges::size(u)
//			- 2 * ghost_point_number;

	std::valarray<Vector4<T>> avg(std::ranges::size(u));
	std::valarray<Vector4<T>> proj_u(std::ranges::size(u));
	std::valarray<Vector4<T>> proj_f(std::ranges::size(u));
	std::valarray<Vector4<T>> res(std::ranges::size(u));
	std::valarray<Vector4<T>> flux = calcPhysFlux(u/*, gamma*/);

//	auto interior_view = std::ranges::views::drop(ghost_point_number)
//			| std::ranges::views::take(n_size)
//			| std::ranges::views::common;
//	auto interior_view_shifted_by_neg_1
//			= std::ranges::views::drop(ghost_point_number - 1)
//				| std::ranges::views::take(n_size)
//				| std::ranges::views::common;

//	auto u_view = std::ranges::views::all(u);
//	auto avg_view = std::ranges::views::all(avg);
//	auto u_proj_view = std::ranges::views::all(proj_u);
//	auto f_proj_view = std::ranges::views::all(proj_f);
//	auto res_view = std::ranges::views::all(res);

//	std::transform(
//				std::execution::par_unseq,
//				std::ranges::begin(u_view) | interior_view,
//				std::ranges::end(u_view) | interior_view,
//				std::ranges::begin(u_view) | interior_view_shifted_by_neg_1,
//				std::ranges::begin(avg_view) | interior_view,
//				[](auto q_l, auto q_r) {
//		return average<Vector4<T>>(q_l, q_r);
//	});
	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(u) + 1,
				std::ranges::end(u),
				std::ranges::begin(u),
				std::ranges::begin(avg) + 1,
				[](auto q_l, auto q_r) {
		return average<Vector4<T>>(q_l, q_r);
	});
	// updateGhostPointsTransmissive(avg);
	// updateGhostPointsPeriodic(avg);

	auto project = [gamma](
			Vector4<T> q_ast, Vector4<T> vec) -> Vector4<T> {
		return projectOntoCharacteristics<T>(q_ast, vec/*, gamma*/);
	};

//	auto deproject = [gamma](
//			Vector4<T> q_ast, Vector4<T> vec) -> Vector4<T> {
//		return projectCharacteristicVariablesBackOntoConserved<T>(q_ast, vec, gamma);
//	};

	calcHydroStageCharWiseFDENO3<T, Vector4<T>>(
				std::ranges::views::all(u),
				std::ranges::views::all(avg),
				std::ranges::views::all(flux),
				res, t,
				project, /*deproject,*/ lam[0],
				ghost_point_number);
//	calcHydroStageCharWiseFDMPENO3<T, Vector4<T>>(
//				std::ranges::views::all(u),
//				std::ranges::views::all(avg),
//				std::ranges::views::all(flux),
//				res, t,
//				project, /*deproject,*/ lam[0],
//				ghost_point_number);
//	calcHydroStageCharWiseFDENO2<T, Vector4<T>>(
//			std::ranges::views::all(u),
//			std::ranges::views::all(avg),
//			std::ranges::views::all(flux),
//			res, t,
//			project, /*deproject,*/ lam[0],
//			ghost_point_number);

//	std::transform(
//				std::execution::par_unseq,
//				std::ranges::begin(avg_view) | interior_view,
//				std::ranges::end(avg_view) | interior_view,
//				std::ranges::begin(res_view) | interior_view,
//				std::ranges::begin(res_view) | interior_view,
//				[gamma](auto q_ast, auto f) {
//		return projectCharacteristicVariablesBackOntoConserved<Vector4<T>>(
//					q_ast, f, gamma);
//	});

	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(avg),
				std::ranges::end(avg),
				std::ranges::begin(res),
				std::ranges::begin(res),
				[gamma](auto q_ast, auto f) {
		return projectCharacteristicVariablesBackOntoConserved<T>(
					q_ast, f/*, gamma*/);
	});
	// updateGhostPointsPeriodic(res);

	return res;
}


template <ArithmeticWith<numeric_val> T, typename... Args>
std::valarray<Vector4<T>> calcdSpaceEu1D(
	const std::valarray<Vector4<T>>& U,
	T t,
	T dx,
	const std::ranges::common_range auto& lam,
	std::size_t ghost_point_number,
	auto&& calcFlux,
	auto&& addSource,
	Args... opts
) {
	const std::size_t n_size = U.size() - 2 * ghost_point_number;
	// const std::size_t ghost_point_number = 4;

	std::valarray<Vector4<T>> dflux(Vector4<T>(0., 0., 0., 0.), U.size());
	std::valarray<Vector4<T>> lf = calcFlux(U, t, lam, n_size, opts...);

	std::slice Nweno(ghost_point_number, n_size, 1);
	std::slice Nweno_shifted_by_neg_1(ghost_point_number - 1, n_size, 1);
//	std::slice Nweno(1, U.size()-1, 1);
//	std::slice Nweno_shifted_by_neg_1(0, U.size()-1, 1);

	std::valarray<Vector4<T>> f_mn = lf[Nweno_shifted_by_neg_1];
	std::valarray<Vector4<T>> f_pl = lf[Nweno];


	dflux[Nweno] = -(f_pl - f_mn) / Vector4<T>(dx);

	// dflux += source terms...
	addSource(U, dflux);

	return dflux;
}


template <ArithmeticWith<numeric_val> T>
void integrateRiemannProblem(
	std::ranges::common_range auto& u,
	std::ranges::common_range auto& flux,
	T t0, T dx, std::size_t ghost_point_number,
	T t_fin,
	auto&& timeStepFunction,
	auto&& inittimeStepFunction,
	T cfl = 0.4
) {
	/* Time Operator of the Riemann problem: perform the time loop
	 * and solve it, storing the result in `u` and the numerical flux
	 * needed for the calculation in `flux`.
	 */

//	std::valarray<Vector4<T>> y2(Vector4<T>::ZERO,
//		std::ranges::size(u));
//	std::valarray<Vector4<T>> y3(Vector4<T>::ZERO,
	const std::size_t n_size = std::ranges::size(u)
			- 2 * ghost_point_number;
//		std::ranges::size(u));
//	std::array<std::reference_wrapper<
//		std::valarray<Vector4<T>>
//	>, 2> fluxes = {
//		std::ref(y2), std::ref(y3)
//	};
	std::valarray<Vector4<T>> y1(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy1(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y2(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy2(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y3(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy3(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y4(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy4(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y5(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy5(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y6(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy6(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y7(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy7(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y8(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy8(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y9(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy9(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y10(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy10(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y11(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy11(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y12(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy12(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y13(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy13(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> y14(Vector4<T>::ZERO, std::ranges::size(u));
	std::valarray<Vector4<T>> dy14(Vector4<T>::ZERO, std::ranges::size(u));
//	std::valarray<Vector4<T>> f_tmp(Vector4<T>::ZERO, std::ranges::size(u));
	std::array<
				std::reference_wrapper<std::valarray<Vector4<T>>
			>, 28> fluxes = {
		std::ref(y1), std::ref(dy1),
		std::ref(y2), std::ref(dy2),
		std::ref(y3), std::ref(dy3),
		std::ref(y4), std::ref(dy4),
		std::ref(y5), std::ref(dy5),
		std::ref(y6), std::ref(dy6),
		std::ref(y7), std::ref(dy7),
		std::ref(y8), std::ref(dy8),
		std::ref(y9), std::ref(dy9),
		std::ref(y10), std::ref(dy10),
		std::ref(y11), std::ref(dy11),
		std::ref(y12), std::ref(dy12),
		std::ref(y13), std::ref(dy13),
		std::ref(y14), std::ref(dy14)
	};

	auto find_max_lam = [](const decltype(u)& u, T dt) -> T {
		return calcMaxWaveSpeedD<T>(u);
	};
	auto dt_upd = [](T dx, T lam) -> T {
		// return std::pow(.1 * dx / lam, /*5./3.*/10.);
		return std::pow(dx, 8./3.);
	};

	auto find_max_lam_init = find_max_lam;

	T dt = .1 * std::pow(dx, /*5./3.*/10.) / 10.;
	std::valarray<T> max_lam = {find_max_lam(u, dt)};
	dt = dt_upd(dx, max_lam[0]);

	std::size_t init_steps = 2 - 1;
	fluxes[0].get() = u;

	for (std::size_t k = 2; k <= 2 * init_steps; k += 2) {
		inittimeStepFunction(
			u, flux, fluxes,
			t0, dt, dx, max_lam, ghost_point_number);
		fluxes[k].get() = u; fluxes[k - 1].get() = flux;
		max_lam[0] = find_max_lam_init(u, dt);
		dt = dt_upd(dx, max_lam[0]); t0 += dt;
	}

//	inittimeStepFunction(
//		u, flux, fluxes,
//		t0, dt, dx, max_lam, ghost_point_number);
//	y2 = u; dy3 = flux;
//	max_lam[0] = find_max_lam(u, dt);
//	dt = dt_upd(dx, max_lam[0]); t0 += dt;

//	inittimeStepFunction(
//		u, flux, fluxes,
//		t0, dt, dx, max_lam, ghost_point_number);
//	y1 = u; dy2 = flux;
//	max_lam[0] = find_max_lam(u, dt);
//	dt = dt_upd(dx, max_lam[0]); t0 += dt;

//	inittimeStepFunction(
//		u, flux, fluxes,
//		t0, dt, dx, max_lam, ghost_point_number);
//	dy1 = flux;
//	max_lam[0] = find_max_lam(u, dt);
//	dt = dt_upd(dx, max_lam[0]); t0 += dt;

	timeOperator<T>(
		u, flux, fluxes, t0, dx, ghost_point_number, t_fin,
		timeStepFunction,
		find_max_lam,
		cfl
	);
}


template <ArithmeticWith<numeric_val> T>
std::function<Vector4<T>(Vector4<T>/*, T*/)> primitiveToConservativeU
	= primitiveToConservative<T>;


template <ArithmeticWith<numeric_val> T>
void prepareRiemannProblem(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	T gamma,
	T rho_left, T v_left, T p_left, T e_left, T rhoE_left,
	T rho_right, T v_right, T p_right, T e_rightR, T rhoE_right,
	T q0, T l_min, T l_max,
	std::function<Vector4<T>(Vector4<T>/*, T*/)> primitiveToConservativeU,
	std::size_t n_ghost_points = 5
) {
	/* Fill u_init with the Riemann problem data
	 * at mesh_size number of nodes using Vector4<T>
	 * to store the data vector field. Fill x
	 * with the corresponding node coordinates.
	 */

	const std::size_t mesh_size = std::ranges::size(u_init)
			- 2 * n_ghost_points;
	std::size_t computational_domain_size = mesh_size;
	std::size_t full_mesh_size = computational_domain_size
		+ 2*n_ghost_points;
	T dx = (l_max - l_min) / (mesh_size - 1.);  // [L]

	x = std::valarray<T>(0., full_mesh_size);
	u_init = std::valarray<Vector4<T>>(
				Vector4<T>::ZERO,
				full_mesh_size);

	std::size_t k = 0;
	for (k = 0; k < n_ghost_points; ++ k)
		x[k] = l_min/* + dx * 0.5*/ - dx * (n_ghost_points - k);

	x[n_ghost_points] = l_min/* + dx * 0.5*/;
	for (k = n_ghost_points + 1; k < full_mesh_size; ++ k)
		x[k] = x[k-1] + dx;

	std::size_t x0_index = 0;
	while (x[x0_index] < q0)
		++ x0_index;

	Vector4<T> vec(rho_left, v_left, p_left, 0.);
	vec = primitiveToConservativeU(vec/*, gamma*/);
	for (k = 0; k < x0_index; ++ k) {
		u_init[k] = vec;
	}

	vec = primitiveToConservativeU(
				Vector4<T>(rho_right, v_right, p_right, 0.)/*,
				gamma*/);
	for (k = x0_index; k < full_mesh_size; ++ k) {
		u_init[k] = vec;
	}
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> solve1DRiemannProblemForEulerEq(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	T gamma,
	T rho_left, T v_left, T p_left, T e_left, T rhoE_left,
	T rho_right, T v_right, T p_right, T e_rightR, T rhoE_right,
	T q0, // Initial coordinate of the discontinuity
	T t0, T t_max, T l_min, T l_max,
	std::function<Vector4<T>(Vector4<T>/*, T*/)> primitiveToConservativeU,
	auto&& calcdSpace,
	auto&& updateGhostPoints,
	std::size_t n_ghost_points = 5, T cfl = 0.4
) {
	/* Solve a given Riemann problem for 1D Euler equations. */

	const std::size_t mesh_size = std::ranges::size(u_init)
			- 2 * n_ghost_points;

	prepareRiemannProblem<T>(
		u_init, x, gamma,
		rho_left, v_left, p_left, e_left, rhoE_left,
		rho_right, v_right, p_right, e_rightR, rhoE_right,
		q0, l_min, l_max, primitiveToConservativeU, n_ghost_points
	);
	updateGhostPoints(u_init);

	std::valarray<Vector4<T>> flux(Vector4<T>::ZERO,
		std::ranges::size(u_init));

	integrateRiemannProblem<T>(u_init, flux,
		t0, (l_max-l_min) / (mesh_size - 1.), n_ghost_points, t_max,
		[&calcdSpace, &updateGhostPoints](
			std::valarray<Vector4<T>>& u,
			std::valarray<Vector4<T>>& dflux,
			std::array<
				std::reference_wrapper<std::valarray<Vector4<T>>
			>, 28>& fluxes,
			T t, T dt, T dx,
			const std::valarray<T>& lam,
			std::size_t n_ghost_points = 5
		) {
			/*advanceTimestepSSPRK10_4<T>(
				u, dflux, fluxes[0].get(),
				t, dt, dx, lam, n_size,
				calcdSpace, updateGhostPoints);*/
//			advanceTimestepTVDRK3<T>(
//				u, dflux, fluxes[3].get(),
//				fluxes[0].get(), fluxes[1].get(),
//				t, dt, dx, lam, n_ghost_points,
//				calcdSpace, updateGhostPoints);
//			advanceTimestepRK6_5<T>(
//				u,
//				dflux,
//				fluxes[0].get(), fluxes[1].get(),
//				fluxes[2].get(), fluxes[3].get(),
//				fluxes[4].get(), fluxes[5].get(),
//				fluxes[6].get(), fluxes[7].get(),
//				fluxes[8].get(), fluxes[9].get(),
//				t, dt, dx, lam,
//				n_ghost_points, calcdSpace, updateGhostPoints);
//			advanceTimestep_eBDF5<T>(
//				u, dflux,
//				fluxes[0].get(), fluxes[1].get(),
//				fluxes[2].get(), fluxes[3].get(),
//				fluxes[4].get(), fluxes[5].get(),
//				fluxes[6].get(), fluxes[7].get(),
//				fluxes[8].get(),
//				t, dt, dx, lam,
//				n_ghost_points, calcdSpace, updateGhostPoints);
			std::array<std::reference_wrapper<
					std::valarray<Vector4<T>>>, 13> us = {
				fluxes[2],
				fluxes[4],
				fluxes[6],
				fluxes[8],
				fluxes[10],
				fluxes[12],
				fluxes[14],
				fluxes[16],
				fluxes[18],
				fluxes[20],
				fluxes[22],
				fluxes[24],
				fluxes[26]
			};
			std::array<std::reference_wrapper<
					std::valarray<Vector4<T>>>, 13> fs = {
				fluxes[3],
				fluxes[5],
				fluxes[7],
				fluxes[9],
				fluxes[11],
				fluxes[13],
				fluxes[15],
				fluxes[17],
				fluxes[19],
				fluxes[21],
				fluxes[23],
				fluxes[25],
				fluxes[27]
			};
			advanceTimestepSSPTSERK12_8(
				u, fluxes[0].get(), dflux,
				us, fs,
				t, dt, dx, lam,
				n_ghost_points, calcdSpace, updateGhostPoints);
			/*EulerForward<T>(
				u, dflux,
				t, dt, dx, lam, n_ghost_points,
				calcdSpace, updateGhostPoints);*/
		},
		[&calcdSpace, &updateGhostPoints](
			std::valarray<Vector4<T>>& u,
			std::valarray<Vector4<T>>& dflux,
			std::array<
				std::reference_wrapper<std::valarray<Vector4<T>>
			>, 28>& fluxes,
			T t, T dt, T dx,
			const std::valarray<T>& lam,
			std::size_t n_ghost_points = 5
		) {
			advanceTimestepTVDRK3<T>(
				u, dflux, fluxes[10].get(),
				fluxes[8].get(), fluxes[9].get(),
				t, dt, dx, lam, n_ghost_points,
				calcdSpace, updateGhostPoints);
		}, cfl);

	return u_init;
}


template <typename T>
void addEmptySource(auto& u, auto& flux, auto&& x) {
	return;
}


template <ArithmeticWith<numeric_val> T>
void prepareHighGradientLaserProblem(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	std::function<Vector4<T>(Vector4<T>)> primitiveToConservativeU,
	T q0 = 1000.e-9, T q1 = 1050.e-9, T l_min = 0., T l_max = 1250.e-9,
	std::size_t n_ghost_points = 6/*, T gamma = 1.4*/
) {
	/* ... */

	const std::size_t mesh_size = std::ranges::size(u_init)
			- 2 * n_ghost_points;
	std::size_t computational_domain_size = mesh_size;
	std::size_t full_mesh_size = computational_domain_size
		+ 2*n_ghost_points;
	T dx = (l_max - l_min) / (mesh_size-1.);  // [L]

	x = std::valarray<T>(0., full_mesh_size);
	u_init = std::valarray<Vector4<T>>(
				Vector4<T>::ZERO,
				full_mesh_size);

	std::size_t k = 0;
	for (k = 0; k < n_ghost_points; ++ k)
		x[k] = l_min/* + dx * 0.5*/ - dx * (n_ghost_points - k);

	x[n_ghost_points] = l_min/* + dx * 0.5*/;
	for (k = n_ghost_points + 1; k < full_mesh_size; ++ k)
		x[k] = x[k-1] + dx;

	std::size_t x0_index = 0;
	while (x[x0_index] < q0)
		++ x0_index;

	Vector4<T> vec(2413., 0., 0., 0.);
	vec = primitiveToConservativeU(vec);
	for (k = 0; k < x0_index; ++ k) {
		u_init[k] = vec;
	}

	std::size_t x1_index = x0_index;
	while (x[x1_index] < q1 && x1_index < full_mesh_size)
		++ x1_index;

	vec = Vector4<T>(2413., 0., 35.6e9, 0.);
	vec = primitiveToConservativeU(vec);
	for (k = x0_index; k < x1_index; ++ k) {
		u_init[k] = vec;
	}

	vec = Vector4<T>(2., 0., 0., 0.);
	vec = primitiveToConservativeU(vec);
	for (k = x1_index; k < full_mesh_size; ++ k) {
		u_init[k] = vec;
	}
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> solve1DHighGradientLaserProblem(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	std::function<Vector4<T>(Vector4<T>/*, T*/)> primitiveToConservativeU,
	std::size_t ghost_point_number = 6, T cfl = 0.4, T t_max = 1.e-12,
	T gamma = 1.4,
	T q0 = 1000.e-9, T q1 = 1050.e-9, T t0 = 0.,
	T l_min = 0., T l_max = 1250.e-9
) {
	/* ... */

	const std::size_t mesh_size = std::ranges::size(u_init)
			- 2 * ghost_point_number;

	prepareHighGradientLaserProblem<T>(
		u_init, x, primitiveToConservativeU,
		q0, q1, l_min, l_max, ghost_point_number
	);

	std::valarray<Vector4<T>> flux(Vector4<T>::ZERO,
		std::ranges::size(u_init));

	auto calcdSpace = [gamma, ghost_point_number](
			std::valarray<Vector4<T>>& u,
			T t, T dx, const std::valarray<T>& max_eigenvalues,
			std::size_t n_size) {
		return calcdSpaceEu1D<T>(
			u, t, dx, max_eigenvalues, ghost_point_number,
			[gamma, ghost_point_number, dx](
					const std::valarray<Vector4<T>>& u,
					T t, const std::valarray<T>& lam,
					std::size_t n_size,
					T eps = 1e-100, T p = 2.) {
//				return calcFluxComponentWiseFDWENO5<T>(
//					u, t, lam, ghost_point_number, gamma, eps, p);
//				return calcFluxComponentWiseFVWENO5<T>(
//					u, t, lam, ghost_point_number, gamma, eps, p);
//				return calcFluxCharacteristicWiseFDWENO5<T>(
//					u, t, lam, ghost_point_number, gamma, eps, p);
//				T dt = 0.00054;
//				return calcFluxCharacteristicWiseFDWENO7<T>(
//					u, t, lam, ghost_point_number, gamma, eps, p);
				return calcFluxCharacteristicWiseFDWENO9<T>(
					u, t, lam, ghost_point_number, gamma, eps, p);
//				return calcFluxCharacteristicWiseFDENO3<T>(
//					u, t, lam, ghost_point_number, gamma, eps, p);
//				return calcFluxComponentWiseFVENO3<T>(
//					u, t, lam, n_size, gamma);
			},
			[](
					const std::valarray<Vector4<T>>& u,
					std::valarray<Vector4<T>>& f,
					std::valarray<Vector4<T>>&& x = {}) {
				addEmptySource<T>(u, f, x);
			},
			1e-100, 2.
		);
	};

	auto updateGhostPoints = [ghost_point_number](
			std::valarray<Vector4<T>>& u) {
		updateGhostPointsTransmissive(u, ghost_point_number);
		// updateGhostPointsPeriodic(u);
	};

	integrateRiemannProblem<T>(u_init, flux,
		t0, (l_max-l_min) / (mesh_size-1.), mesh_size, t_max,
	   [&calcdSpace, &updateGhostPoints](
		   std::valarray<Vector4<T>>& u,
		   std::valarray<Vector4<T>>& dflux,
		   std::array<
			   std::reference_wrapper<std::valarray<Vector4<T>>
		   >, 28>& fluxes,
		   T t, T dt, T dx,
		   const std::valarray<T>& lam,
		   std::size_t n_ghost_points = 6
	   ) {
		   /*advanceTimestepSSPRK10_4<T>(
			   u, dflux, fluxes[0].get(),
			   t, dt, dx, lam, n_size,
			   calcdSpace, updateGhostPoints);*/
//			advanceTimestepTVDRK3<T>(
//				u, dflux, fluxes[3].get(),
//				fluxes[0].get(), fluxes[1].get(),
//				t, dt, dx, lam, n_ghost_points,
//				calcdSpace, updateGhostPoints);
//			advanceTimestepRK6_5<T>(
//				u,
//				dflux,
//				fluxes[0].get(), fluxes[1].get(),
//				fluxes[2].get(), fluxes[3].get(),
//				fluxes[4].get(), fluxes[5].get(),
//				fluxes[6].get(), fluxes[7].get(),
//				fluxes[8].get(), fluxes[9].get(),
//				t, dt, dx, lam,
//				n_ghost_points, calcdSpace, updateGhostPoints);
//			advanceTimestep_eBDF5<T>(
//				u, dflux,
//				fluxes[0].get(), fluxes[1].get(),
//				fluxes[2].get(), fluxes[3].get(),
//				fluxes[4].get(), fluxes[5].get(),
//				fluxes[6].get(), fluxes[7].get(),
//				fluxes[8].get(),
//				t, dt, dx, lam,
//				n_ghost_points, calcdSpace, updateGhostPoints);
		   std::array<std::reference_wrapper<
				   std::valarray<Vector4<T>>>, 13> us = {
			   fluxes[2],
			   fluxes[4],
			   fluxes[6],
			   fluxes[8],
			   fluxes[10],
			   fluxes[12],
			   fluxes[14],
			   fluxes[16],
			   fluxes[18],
			   fluxes[20],
			   fluxes[22],
			   fluxes[24],
			   fluxes[26]
		   };
		   std::array<std::reference_wrapper<
				   std::valarray<Vector4<T>>>, 13> fs = {
			   fluxes[3],
			   fluxes[5],
			   fluxes[7],
			   fluxes[9],
			   fluxes[11],
			   fluxes[13],
			   fluxes[15],
			   fluxes[17],
			   fluxes[19],
			   fluxes[21],
			   fluxes[23],
			   fluxes[25],
			   fluxes[27]
		   };
		   advanceTimestepSSPTSERK12_8(
			   u, fluxes[0].get(), dflux,
			   us, fs,
			   t, dt, dx, lam,
			   n_ghost_points, calcdSpace, updateGhostPoints);
		   /*EulerForward<T>(
			   u, dflux,
			   t, dt, dx, lam, n_ghost_points,
			   calcdSpace, updateGhostPoints);*/
	   },
	   [&calcdSpace, &updateGhostPoints](
		   std::valarray<Vector4<T>>& u,
		   std::valarray<Vector4<T>>& dflux,
		   std::array<
			   std::reference_wrapper<std::valarray<Vector4<T>>
		   >, 28>& fluxes,
		   T t, T dt, T dx,
		   const std::valarray<T>& lam,
		   std::size_t n_ghost_points = 6
	   ) {
		   advanceTimestepTVDRK3<T>(
			   u, dflux, fluxes[10].get(),
			   fluxes[8].get(), fluxes[9].get(),
			   t, dt, dx, lam, n_ghost_points,
			   calcdSpace, updateGhostPoints);
	   }, cfl);

	return u_init;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> solve1DRiemannProblemForEulerEq(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	T gamma,
	T rho_left, T v_left, T p_left,
	T rho_right, T v_right, T p_right,
	T q0, // Initial coordinate of the discontinuity
	T t0, T t_max, T l_min, T l_max,
	std::function<Vector4<T>(Vector4<T>/*, T*/)> primitiveToConservativeU,
	std::size_t ghost_point_number = 5, T cfl = 0.4
) {
	/* Solve a given Riemann problem for 1D Euler equations. */

//	const std::size_t mesh_size = std::ranges::size(u_init)
//			- 2 * ghost_point_number;

	T e_left = p_left / (gamma - 1.) / rho_left;
	T e_right = p_right / (gamma - 1.) / rho_right;
	T E_left = (e_left + (v_left*v_left)/2.);
	T E_right = (e_right + (v_right*v_right)/2.);

	return solve1DRiemannProblemForEulerEq(
		u_init, x, gamma,
		rho_left, v_left, p_left, e_left, E_left,
		rho_right, v_right, p_right, e_right, E_right,
		q0, t0, t_max, l_min, l_max,
		primitiveToConservativeU,
		[gamma, ghost_point_number](
				std::valarray<Vector4<T>>& u,
				T t, T dx, const std::valarray<T>& max_eigenvalues,
				std::size_t n_size) {
			return calcdSpaceEu1D<T>(
				u, t, dx, max_eigenvalues, ghost_point_number,
				[gamma, ghost_point_number, dx](
						const std::valarray<Vector4<T>>& u,
						T t, const std::valarray<T>& lam,
						std::size_t n_size,
						T eps = 1e-40, T p = 2.) {
//					return calcFluxComponentWiseFDWENO5<T>(
//						u, t, lam, ghost_point_number, gamma, eps, p);
//					return calcFluxComponentWiseFVWENO5<T>(
//						u, t, lam, ghost_point_number, gamma, eps, p);
//					return calcFluxCharacteristicWiseFDWENO5<T>(
//						u, t, lam, ghost_point_number, gamma, eps, p);
//					T dt = 0.00054;
//					return calcFluxCharacteristicWiseFDWENO7<T>(
//						u, t, lam, ghost_point_number, gamma, eps, p);
					return calcFluxCharacteristicWiseFDWENO9<T>(
						u, t, lam, ghost_point_number, gamma, eps, p);
//					return calcFluxCharacteristicWiseFDENO3<T>(
//						u, t, lam, ghost_point_number, gamma, eps, p);
//					return calcFluxComponentWiseFVENO3<T>(
//						u, t, lam, n_size, gamma);
				},
				[](
						const std::valarray<Vector4<T>>& u,
						std::valarray<Vector4<T>>& f,
						std::valarray<Vector4<T>>&& x = {}) {
					addEmptySource<T>(u, f, x);
				},
				/*1e-4*//*1e-6*//*1e-8*//*1e-40*/1e-100,
				/*1.*/2./*3.*//*4.*//*5.*//*6.*/
			);
		},
		[ghost_point_number](std::valarray<Vector4<T>>& u) {
			updateGhostPointsTransmissive(u, ghost_point_number);
			// updateGhostPointsPeriodic(u);
		},
		ghost_point_number, cfl
	);
}

#endif // EULER1D_H
