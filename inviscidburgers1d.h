#include <algorithm>
#include <array>
#include <execution>
#include <ranges>
#include <span>
#include <string>
#include <valarray>
// #include <vector>

#include "arithmeticwith.h"
#include "weno5.h"
#include "eno3.h"
#include "lf_flux.h"
#include "rk6_5.h"
#include "ssprk33.h"
#include "ssprk10_4.h"
#include "tdrk3_5.h"


template <ArithmeticWith<numeric_val> T>
T calcInviscidBurgersFlux(T u) {
	/* Calculate inviscid Bateman-Burgers' flux.
	 */

//	return u * u * u - std::sin(u);
	return u * u * 0.5;
}


template <ArithmeticWith<numeric_val> T>
T calcLinearAdvectionFlux(T u, T alpha = 1.) {
	/* Calculate linear advection flux.
	 */

	return alpha * u;
}


template <ArithmeticWith<numeric_val> T>
T calcInviscidBurgersFluxDerivative(T u) {
	/* Calculate inviscid Bateman-Burgers' flux derivative.
	 */

	return u;
}


template <ArithmeticWith<numeric_val> T>
T calcLinearAdvectionDerivative(T u, T alpha = 1.) {
	/* Calculate linear advection flux derivative.
	 */

	return alpha;
}


template <ArithmeticWith<numeric_val> T>
T calcInviscidBurgersMaxWaveSpeed(
		const /* std::ranges::common_range auto */std::valarray<T>& u_arr
		) {
	/* Calculate max |df/du| for the inviscid Bateman-Burgers' eq'n. */

	// return *std::ranges::max_element(std::abs(u_arr));
//	return std::ranges::max(3. * u_arr * u_arr - std::sin(u_arr));
	return std::abs(calcInviscidBurgersFluxDerivative(u_arr)).max();
}


template <ArithmeticWith<numeric_val> T>
T calcLinearAdvectionMaxWaveSpeed(
		const /* std::ranges::common_range auto */std::valarray<T>& u_arr
		) {
	/* Calculate max |df/du| for the linear advection eq'n. */

	return 1.;
}


template <ArithmeticWith<numeric_val> T>
T BurgersSource(T x) {
	/* Define the source term S for the Burgers'
	 * equation u_t + F_x = S.
	 */

	return 0.;
}



template <ArithmeticWith<numeric_val> T>
void calcInviscidBurgersPreciseGodunovNumericalFlux(
		const std::ranges::common_range auto& u_p,
		const std::ranges::common_range auto& u_m,
		std::ranges::common_range auto& res_f,
		auto&& calcPhysFlux) {
	/* Compute Godunov flux for the inviscid Bateman-Burgers'
	 * (Hopf) equation.
	 */

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(u_p),
				std::ranges::end(u_p),
				std::ranges::begin(u_m),
				std::ranges::begin(res_f),
				[&calcPhysFlux](T u_p_pt, T u_m_pt) {
		T res = 0.;
		if (u_m_pt <= u_p_pt && u_m_pt > 0.) {
			res = calcPhysFlux(u_m_pt);
		} else if (u_m_pt <= u_p_pt && u_m_pt < 0. && u_p_pt > 0.) {
			res = calcPhysFlux(0.);
		} else if (u_m_pt <= u_p_pt && u_p_pt < 0.) {
			res = calcPhysFlux(u_p_pt);
		} else if (u_m_pt > u_p_pt && u_p_pt > 0.) {
			res = calcPhysFlux(u_m_pt);
		} else if (u_m_pt > u_p_pt && u_p_pt < 0. && u_m_pt > 0.) {
			if (std::abs(u_p_pt) > std::abs(u_m_pt)) {
				res = calcPhysFlux(u_p_pt);
			} else {
				res = calcPhysFlux(u_m_pt);
			}
		} else if (u_m_pt > u_p_pt && u_m_pt < 0.)
			res = calcPhysFlux(u_p_pt);

		return res;
	});
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> calcFluxInviscidBurgersFVENO3(
		const std::ranges::common_range auto U,
		T t, const std::ranges::common_range auto& lam,
		std::size_t n_size) {
	std::valarray<T> res(0., std::ranges::size(U));

	// ----FD variant----
	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(U),
				std::ranges::end(U),
				std::ranges::begin(res),
				[](T u) {
		return calcInviscidBurgersFlux<T>(u);
		// return calcLinearAdvectionFlux<T>(u);
	});
//	for (std::size_t k = 0; k < std::ranges::size(U); ++ k)
//		res[k] = calcInviscidBurgersFlux<T>(U[k]);

	// std::array<std::valarray<T>, 2> fs {
	// 			std::valarray<T>(0., U.size()),
	// 			std::valarray<T>(0., U.size())
	// };

	auto [monotone_flux_component_pl,
		monotone_flux_component_mn]  = splitFluxAsLaxFriedrichs<T>(
			U, res, lam[0]);

	// updateGhostPointsPeriodic(monotone_flux_component_pl);
	// updateGhostPointsPeriodic(monotone_flux_component_mn);
	updateGhostPointsTransmissive(monotone_flux_component_pl);
	updateGhostPointsTransmissive(monotone_flux_component_mn);

	calcHydroStageFDWENO5FM<T>(
				std::ranges::views::all(monotone_flux_component_pl),
				std::ranges::views::all(monotone_flux_component_mn),
				t, res, 3, 1e-40, 2.);
//	calcHydroStageFDENO3<T>(
//				std::ranges::views::all(monotone_flux_component_pl),
//				std::ranges::views::all(monotone_flux_component_mn),
//				t, res, n_size);

	// ----FV variant----
//	std::valarray<T> u_plus = std::valarray<T>(0., U.size());
//	std::valarray<T> u_minus = std::valarray<T>(0., U.size());
//	calcHydroStageFVWENO5FM<T>(
//				std::ranges::views::all(U),
//				t, u_plus, u_minus,
//				3, 1e-40, 2.);
//	calcHydroStageFVENO3<T>(std::ranges::views::all(U),
//							t, u_plus, u_minus, n_size);

//	calcLaxFriedrichsNumericalFlux(u_plus, u_minus, res,
//		[](const T u) {
//			return calcInviscidBurgersFlux<T>(u);
//			// return calcLinearAdvectionFlux<T>(u);
//		},
//		lam[0]);

//	calcInviscidBurgersPreciseGodunovNumericalFlux<T>(
//				u_plus, u_minus, res,
//				[](const T u) {
//		return calcInviscidBurgersFlux<T>(u);
//		// return calcLinearAdvectionFlux<T>(u);
//	});
//	updateGhostPointsTransmissive(res);

	return res;
}


template <ArithmeticWith<numeric_val> T, typename... Args>
std::valarray<T> calcdSpaceInviscidBurgers(
	std::span<T> const u, std::span<T> const x, T t, T dx,
	const std::ranges::common_range auto& lam,
	std::size_t n_size,
	auto&& calcFlux,
	auto&& addSource,
	Args... opts
) {
	std::valarray<T> dflux(0., u.size());
	std::valarray<T> lf = calcFlux(u, t, lam, n_size, opts...);

	const std::size_t ghost_point_number = 3;

//	std::slice Nweno(ghost_point_number, n_size, 1);
//	std::slice Nweno_shifted_by_neg_1(ghost_point_number-1, n_size, 1);

//	std::valarray<T> f_mn = lf[Nweno_shifted_by_neg_1];
//	std::valarray<T> f_pl = lf[Nweno];


//	dflux[Nweno] = -(f_pl - f_mn) / dx;
	auto interior_view = std::views::drop(ghost_point_number)
			| std::views::take(n_size)
			| std::ranges::views::common;
	auto interior_view_shifted_by_neg_1
			= std::views::drop(ghost_point_number - 1)
				| std::views::take(n_size)
				| std::ranges::views::common;

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(lf | interior_view),
				std::ranges::end(lf | interior_view),
				std::ranges::begin(lf | interior_view_shifted_by_neg_1),
				std::ranges::begin(dflux | interior_view),
				[dx](const auto f_pl, const auto f_mn) {
		return -(f_pl - f_mn) / dx;
	});

	// std::cout << dx << "\n";

	// dflux += source term:
	addSource(u, dflux, x);

	return dflux;
}


template <ArithmeticWith<numeric_val> T>
void integrate1DInviscidBurgersProblem(
	std::ranges::common_range auto& u_sol,
	std::ranges::common_range auto& flux,
	T t0, T dx, std::size_t n_size, T t_fin,
	auto&& timeStepFunction,
	T cfl = 0.4
) {
	/* Time Operator of the inviscid Bateman-Burgers' problem:
	 * perform the time loop and solve it, storing the result
	 * in `u_sol` and the numerical flux needed for the
	 * calculation in `flux`.
	 */

	std::valarray<T> y1(0., std::ranges::size(u_sol));
	std::valarray<T> dy1(0., std::ranges::size(u_sol));
	std::valarray<T> y2(0., std::ranges::size(u_sol));
	std::valarray<T> dy2(0., std::ranges::size(u_sol));
	std::valarray<T> y3(0., std::ranges::size(u_sol));
	std::valarray<T> dy3(0., std::ranges::size(u_sol));
	std::valarray<T> y4(0., std::ranges::size(u_sol));
	std::valarray<T> dy4(0., std::ranges::size(u_sol));
	std::valarray<T> y5(0., std::ranges::size(u_sol));
	std::valarray<T> dy5(0., std::ranges::size(u_sol));
	std::array<std::reference_wrapper<std::valarray<T>>, 10> fluxes = {
		std::ref(y1), std::ref(dy1),
		std::ref(y2), std::ref(dy2),
		std::ref(y3), std::ref(dy3),
		std::ref(y4), std::ref(dy4),
		std::ref(y5), std::ref(dy5)
	};

	timeOperator<T>(
		u_sol, flux, fluxes, t0, dx, n_size, t_fin,
		timeStepFunction,
		[](const decltype(u_sol)& u, T dt) {
			T D = calcInviscidBurgersMaxWaveSpeed<T>(u);
			// T D = calcLinearAdvectionMaxWaveSpeed<T>(u);

			return /* 1.1 * */ D;
		},
		cfl
	);
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> solve1DInviscidBurgersProblem(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	T t0 = 0., T t_max = 0.5, T l_min = -1., T l_max = +1.,
	std::size_t mesh_size = 201, T cfl = 0.4,
	char type = 0
) {
	T t = t0;
	const std::size_t number_of_ghost_points = 3;

	// std::size_t computational_domain_size = mesh_size;
	std::size_t full_mesh_size = mesh_size + 2 * number_of_ghost_points;

	T dx = (l_max - l_min) / (mesh_size - 1.);  // [L]
	// T dx = (l_max - l_min) / mesh_size;  // [L]

	x = std::valarray<T>(0., full_mesh_size);
	u_init = std::valarray<T>(0., full_mesh_size);

	std::size_t k = 0;
	for (k = 0; k < number_of_ghost_points; ++ k)
		x[k] = l_min - dx * (number_of_ghost_points - k);

	x[number_of_ghost_points] = l_min;
	for (k = number_of_ghost_points + 1; k < full_mesh_size; ++ k)
		x[k] = x[k-1] + dx;
	// x[5+3] = 0.;

	for (k = number_of_ghost_points;
			k < mesh_size + number_of_ghost_points;
			++ k) {
		switch (type) {
			case 0:
				u_init[k] = 0.;
			break;
			case 10:
				u_init[k] = (std::abs(x[k]) < 1 / 3.0) * 1.
					+ (std::abs(x[k]) >= 1 / 3.0) * (0.);
			break;
			case 20:
				u_init[k] = 0.5 + std::sin(
							std::numbers::pi_v<T> * x[k]);
			break;
			case 1:  // Toro-1 for linear advection
				u_init[k] = 1.0 * std::exp(-8.0 * x[k] * x[k]);
			break;
			case 2:  // Toro-2 for linear advection
				u_init[k] = (x[k] >= 0.3) * (x[k] <= 0.7) * 1.;
			break;
			case 53:  // Henrick et al. 5.3. Linear advection example
				u_init[k] = std::sin(
							std::numbers::pi_v<T> * x[k] - std::sin(
								std::numbers::pi_v<T> * x[k]
								) / std::numbers::pi_v<T>);
//				u_init[k] = sinq(
//							M_PIq * x[k] - sinq(
//								M_PIq * x[k]
//								) / M_PIq);
			break;
			case 21:  // Evstigneev's first test for Hopf's eq'n
				u_init[k] = std::pow(std::sin(x[k]), 9);
			break;
		}
	}
	updateGhostPointsTransmissive(u_init);
	// updateGhostPointsPeriodic(u_init);
//	u_init[0] = -2;
//	u_init[1] = -1;
//	u_init[2] = 0;
//	u_init[3] = 1;
//	u_init[4] = 2;
//	u_init[5] = 3;
//	u_init[6] = 4;
//	u_init[7] = 5;
//	u_init[8] = 6;
//	u_init[9] = 7;
	dx = (l_max - l_min) / (mesh_size - 1.);
	// dx = (l_max - l_min) / mesh_size;

	auto spaceOp = [&x](
			std::span<T> const u,
			T t, T dx, const std::valarray<T>& max_eigenvalues,
			T n_size) {
		return calcdSpaceInviscidBurgers<T>(
				std::ranges::views::all(u), std::span{x},
				t, dx, max_eigenvalues, n_size,
				[](
					std::span<T> u,
					T t, const std::valarray<T>& lam,
					std::size_t n_size/* ,
					T eps = 1e-40, T p = 2. */
				) {
					return calcFluxInviscidBurgersFVENO3<T>(
						u, t, lam, n_size);
				},
				[](
					std::span<T> const u,
					std::valarray<T>& f,
					std::span<T> x
				) {
					// std::ranges::transform(x, f, std::begin(f),
					// 		[](const auto x_el,
					// 				const auto f_el) {
					// 			return f_el + BurgersSource<T>(
					// 					x_el);
					// });
					return;
				}/* ,
				1e-40, 2. */
		);
	};

	auto updateGhostPoints = [&x, number_of_ghost_points, dx](
			std::valarray<T>& u) {
		updateGhostPointsTransmissive(u);
		// updateGhostPointsPeriodic(u);
	};

	std::valarray<T> flux(0., std::ranges::size(u_init));
	std::valarray<T> ddflux(0., std::ranges::size(u_init));
	std::valarray<T> ddflux_temp(1., std::ranges::size(u_init));

	auto secondDerivative = [&ddflux_temp](
			std::span<T> const u, std::span<T> const du,
			T t, T dx, const std::valarray<T>& max_eigenvalues,
			T n_size) {
//		auto a0 = std::ranges::views::transform(
//					[](auto u_pt) {
//						return calcInviscidBurgersFluxDerivative(u_pt);
//					});
		ddflux_temp = std::valarray<T>(std::ranges::size(ddflux_temp));
		ddflux_temp = calcFluxInviscidBurgersFVENO3<T>(
					du, t, max_eigenvalues, n_size);

		std::transform(
					std::execution::par_unseq,
					std::ranges::begin(du),
					std::ranges::end(du),
					std::ranges::begin(ddflux_temp),
					std::ranges::begin(ddflux_temp),
					[](auto ddu, auto du) {
						return ddu * du;
					});

		return ddflux_temp;
	};

	integrate1DInviscidBurgersProblem<T>(u_init, flux,
		t0, dx, mesh_size, t_max,
		[&spaceOp, &secondDerivative, &updateGhostPoints, &ddflux](
			std::valarray<T>& u,
			std::valarray<T>& dflux,
			std::array<
				std::reference_wrapper<std::valarray<T>
			>, 10> fluxes,
			T t, T dt, T dx,
			const std::valarray<T>& lam,
			std::size_t n_size
		) {
			advanceTimestepTVDRK3<T>(
				u, dflux, fluxes[0].get(), fluxes[1].get(),
				t, dt, dx, lam,
				3, spaceOp, updateGhostPoints);
//			advanceTimestepSSPRK10_4<T>(
//				u, dflux, fluxes[0].get(),
//				t, dt, dx, lam,
//				n_size, spaceOp, updateGhostPoints);
//			 advanceTimestepTDRK3_5<T>(
//				u, dflux, ddflux,
//				fluxes[0].get(), fluxes[1].get(),
//				fluxes[2].get(), fluxes[3].get()cccccc,
//				t, dt, dx, lam,
//				n_size, spaceOp, secondDerivative,
//				updateGhostPoints);
//			advanceTimestepRK6_5<T>(
//				u,
//				dflux,
//				fluxes[0].get(), fluxes[1].get(),
//				fluxes[2].get(), fluxes[3].get(),
//				fluxes[4].get(), fluxes[5].get(),
//				fluxes[6].get(), fluxes[7].get(),
//				fluxes[8].get(), fluxes[9].get(),
//				t, dt, dx, lam,
//				3, spaceOp, updateGhostPoints);
		}, cfl);

	return u_init;
}
