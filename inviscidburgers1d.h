#include <algorithm>
#include <array>
#include <execution>
#include <ranges>
#include <span>
#include <valarray>
// #include <vector>

#include "weno5.h"
#include "eno3.h"
#include "lf_flux.h"
#include "ssprk33.h"
// #include "tdrk3_5.h"


template <ArithmeticWith<numeric_val> T>
T calcInviscidBurgersFlux(T u) {
	/* Calculate inviscid Bateman-Burgers' flux.
	 */

    return u * u * 0.5;
}


template <ArithmeticWith<numeric_val> T>
T calcInviscidBurgersFluxDerivative(T u) {
	/* Calculate inviscid Bateman-Burgers' flux derivative.
	 */

	return u;
}


template <ArithmeticWith<numeric_val> T>
T calcInviscidBurgersMaxWaveSpeed(
		const /* std::ranges::common_range auto */std::valarray<T>& u_arr
		) {
	/* Calculate max |df/du| for the inviscid Bateman-Burgers' eq'n. */

	// return *std::ranges::max_element(std::abs(u_arr));
	return std::abs(calcInviscidBurgersFluxDerivative(u_arr)).max();
}


template <typename T>
T BurgersSource(T x) {
	/* Define the source term S for the Burgers'
	 * equation u_t + F_x = S.
	 */

	return 0.;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> calcFluxInviscidBurgersFVENO3(
		const std::ranges::common_range auto U,
		T t, const std::ranges::common_range auto& lam,
		std::size_t n_size) {
	// std::valarray<Vector4<T>> res = calcPhysFlux(U);

	std::valarray<T> u_plus = std::valarray<T>(0., U.size());
	std::valarray<T> u_minus = std::valarray<T>(0., U.size());

	calcHydroStageFVWENO5FM<T>(std::ranges::views::all(U),
							t, u_plus, u_minus, 3, 1e-40, 2.);
//	calcHydroStageFVENO3<T>(std::ranges::views::all(U),
//							t, u_plus, u_minus, n_size);

	std::valarray<T> res(0., std::ranges::size(U));
	calcLaxFriedrichsNumericalFlux(u_plus, u_minus, res,
		[](const T u) {
			return calcInviscidBurgersFlux<T>(u);
		},
		lam[0]);

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

	std::slice Nweno(ghost_point_number, n_size, 1);
	std::slice Nweno_shifted_by_neg_1(ghost_point_number-1, n_size, 1);

	std::valarray<T> f_mn = lf[Nweno_shifted_by_neg_1];
	std::valarray<T> f_pl = lf[Nweno];


	dflux[Nweno] = -(f_pl - f_mn) / dx;
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

	std::valarray<T> y2(0., std::ranges::size(u_sol));
	std::valarray<T> dy2(0., std::ranges::size(u_sol));
	std::valarray<T> y3(0., std::ranges::size(u_sol));
	std::valarray<T> dy3(0., std::ranges::size(u_sol));
	std::array<std::reference_wrapper<std::valarray<T>>, 4> fluxes = {
		std::ref(y2), std::ref(dy2), std::ref(y3), std::ref(dy3)
	};

	timeOperator<T>(
		u_sol, flux, fluxes, t0, dx, n_size, t_fin,
		timeStepFunction,
		[](const decltype(u_sol)& u, T dt) {
			T D = calcInviscidBurgersMaxWaveSpeed<T>(u);

			return D;
		},
		cfl
	);
}


template <typename T>
std::valarray<T> solve1DInviscidBurgersProblem(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	T t0 = 0., T t_max = 0.5, T l_min = -1., T l_max = +1.,
	std::size_t mesh_size = 201, T cfl = 0.4
) {
	T t = t0;
	const std::size_t number_of_ghost_points = 3;

	// std::size_t computational_domain_size = mesh_size;
	std::size_t full_mesh_size = mesh_size + 2 * number_of_ghost_points;

	T dx = (l_max - l_min) / (mesh_size - 1.);  // [L]

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
//		u_init[k] = (std::abs(x[k]) < 1 / 3.0) * 1.
//			+ (std::abs(x[k]) >= 1 / 3.0) * (0.);
		u_init[k] = 0.5 + std::sin(std::numbers::pi_v<T> * x[k]);
	}
	updateGhostPointsTransmissive<T>(u_init);
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
	dx = (l_max - l_min) / (mesh_size);

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
		updateGhostPointsTransmissive<T>(u);
	};

	std::valarray<T> flux(0., std::ranges::size(u_init));
	std::valarray<T> ddflux(0., std::ranges::size(u_init));
	std::valarray<T> ddflux_temp(0., std::ranges::size(u_init));

	auto secondDerivative = [&ddflux_temp](
			std::span<T> const u, std::span<T> const du,
			T t, T dx, const std::valarray<T>& max_eigenvalues,
			T n_size) {
//		auto a0 = std::ranges::views::transform(
//					[](auto u_pt) {
//						return calcInviscidBurgersFluxDerivative(u_pt);
//					});

		std::transform(
					std::execution::par_unseq,
					std::ranges::begin(du),
					std::ranges::end(du),
					std::ranges::begin(ddflux_temp),
					[](auto dup) { return dup * dup; }
					);

		return ddflux_temp;
	};

	integrate1DInviscidBurgersProblem<T>(u_init, flux,
		t0, dx, mesh_size, t_max,
		[&spaceOp, &secondDerivative, &updateGhostPoints, &ddflux](
			std::valarray<T>& u,
			std::valarray<T>& dflux,
			std::array<
				std::reference_wrapper<std::valarray<T>
			>, 4> fluxes,
			T t, T dt, T dx,
			const std::valarray<T>& lam,
			std::size_t n_size
		) {
			advanceTimestepTVDRK3<T>(
				u, dflux, fluxes[0].get(), fluxes[1].get(),
				t, dt, dx, lam,
				n_size, spaceOp, updateGhostPoints);
//			advanceTimestepTDRK3_5<T>(
//				u, dflux, ddflux,
//				fluxes[0].get(), fluxes[1].get(),
//				fluxes[2].get(), fluxes[3].get(),
//				t, dt, dx, lam,
//				n_size, spaceOp, secondDerivative,
//				updateGhostPoints);
		}, cfl);

	return u_init;
}
