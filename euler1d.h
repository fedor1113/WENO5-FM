#ifndef EULER1D_H
#define EULER1D_H

#include <algorithm>
#include <execution>
#include <ranges>

#include "_vector4.h"
#include "arithmeticwith.h"
#include "weno5.h"
#include "eno3.h"
#include "exactsolver.h"
#include "ssprk33.h"
// #include "eulerforward.h"


template <ArithmeticWith<numeric_val> T>
T gete(T rho, T p, T gamma) {
	if (rho != 0.)
		return p / (gamma - 1.) / rho;

	return 0.;
}


template <ArithmeticWith<numeric_val> T>
T getP(T rho, T e, T gamma) {
	return (gamma - 1.) * rho * e;
}


//template <ArithmeticWith<numeric_val> T>
//T eFromConservative(T rho, T j, T rhoE) {
//	return (rhoE - 0.5 * j*j / rho) / rho;
//}


template <ArithmeticWith<numeric_val> T>
Vector4<T> calcPhysicalFlux(T rho, T u, T p, T last, T gamma) {
	/* Calculate for a Vector4<T> of conserved variables for
	 * the 1D Euler equations its corresponding flux.
	 */

	if (rho == 0) return Vector4<T>::ZERO;

	T e = gete(rho, p, gamma);
	return Vector4<T>(rho * u, p + rho*u*u,
				u*(p + rho*(e + 0.5*u*u)),
				u*last);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> primitiveToConservative(Vector4<T> u, T gamma = 1.4) {
	// Conservative variables
	return Vector4<T>(u[0], u[0] * u[1],
		u[2] / (gamma - 1.) + 0.5 * u[0] * u[1] * u[1], u[3]);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> conservativeToPrimitive(Vector4<T> q, T gamma) {
	// Primitive variables
	T rho = q[0];
	T u = q[1] / rho;
	T E = q[2] / rho;
	T p = (gamma - 1.) * rho * (E - 0.5*u*u);
	T e = gete(rho, p, gamma);

	return Vector4<T>(rho, u, p, e);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> calcPhysicalFluxFromConservativeVec(
		Vector4<T> u,
		T gamma) {
	/* Calculate for a Vector4<T> of conserved variables
	 * for the 1D Euler equations (rho, j=rho*v, rhoE=rho(e+v^2/2), smth)
	 * its corresponding flux.
	 */

//	return calcPhysicalFlux(u[0],
//			u[1] / u[0],
//			getP(u[0], eFromConservative(u[0], u[1], u[2])));
	Vector4<T> prim = conservativeToPrimitive(u, gamma);

	return calcPhysicalFlux(prim[0], prim[1], prim[2], prim[3],
		gamma);
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcPhysFlux(
		const std::valarray<Vector4<T>>& u_arr, T gamma = 1.4) {
	/* Calculate physical fluxes of conservative variables
	 * in all points in the computation domain u_arr.
	 */

	std::valarray<Vector4<T>> res(std::ranges::size(u_arr));

	std::transform(std::ranges::begin(u_arr),
				   std::ranges::end(u_arr),
				   std::ranges::begin(res),
				   [&gamma](Vector4<T> v) {
		return calcPhysicalFluxFromConservativeVec<T>(v, gamma);
	});

	return res;  // f_arr
}


template <ArithmeticWith<numeric_val> T>
T calcSquareSoundSpeed(T rho, T rho_v, T rho_E, T gamma = 1.4) {
	/* Compute the square of sound speed. */

	return gamma * (gamma - 1.)
			* (rho_E - rho_v * rho_v * 0.5 / rho) / rho;
}


template <ArithmeticWith<numeric_val> T>
T calcMaxWaveSpeedD(
		const std::ranges::common_range auto& u_arr,
		T gamma = 1.4) {
	/* Calculate |df/du| for 1D Euler eq'ns. */

	return std::ranges::max(std::ranges::transform_view(
		std::as_const(u_arr),
		[gamma](const auto& u_arr_vec_pt) -> T {
			return (std::sqrt(std::abs(calcSquareSoundSpeed(
								u_arr_vec_pt[0],
								u_arr_vec_pt[1],
								u_arr_vec_pt[2], gamma)))
					+ std::abs(u_arr_vec_pt[1] / u_arr_vec_pt[0]));
		}
	));
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxComponentWiseWENO5(
		const std::ranges::common_range auto& U,
		T t, const std::ranges::common_range auto& lam,
		std::size_t n_size,
		T eps = 1e-40,
		T p = 2.) {
	std::valarray<Vector4<T>> res = calcPhysFlux(U);

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
		calcHydroStageFDWENO5FM<T>(
			U | std::ranges::views::transform(kth_vector_component),
			t, lam[0],
			res | std::ranges::views::transform(kth_vector_component),
			n_size, eps, p
			);
	});

	return res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<Vector4<T>> calcFluxComponentWiseFVENO3(
		const std::ranges::common_range auto& U,
		T t, const std::ranges::common_range auto& lam,
		std::size_t n_size, T gamma = 1.4) {
	// std::valarray<Vector4<T>> res = calcPhysFlux(U);

	std::valarray<Vector4<T>> u_plus(Vector4<T>::ZERO, U.size());
	std::valarray<Vector4<T>> u_minus(Vector4<T>::ZERO, U.size());

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
		calcHydroStageENO3<T>(
			U | std::ranges::views::transform(kth_vector_component),
			t,
			u_plus | std::ranges::views::transform(kth_vector_component),
			u_minus | std::ranges::views::transform(kth_vector_component),
			n_size
			);
	});

	std::valarray<Vector4<T>> res(Vector4<T>::ZERO, U.size());
//	calcLaxFriedrichsNumericalFlux(u_plus, u_minus, res,
//		[gamma](const Vector4<T> u) {
//			return calcPhysicalFluxFromConservativeVec<T>(u, gamma);
//		},
//		lam[0]);
	std::transform(
				std::ranges::begin(u_plus), std::ranges::end(u_plus),
				std::ranges::begin(u_minus),
				std::ranges::begin(res),
				[gamma](const auto u_pl, const auto u_mn) {
					return calcExactFlux<T>(
							u_pl[0], u_pl[1], u_pl[2],
							u_mn[0], u_mn[1], u_mn[2], gamma);
				});

	return res;
}


template <ArithmeticWith<numeric_val> T, typename... Args>
std::valarray<Vector4<T>> calcdSpaceEu1D(
	const std::valarray<Vector4<T>>& U,
	T t,
	T dx,
	const std::ranges::common_range auto& lam,
	std::size_t nSize,
	auto&& calcFlux,
	auto&& addSource,
	Args... opts
) {
	std::valarray<Vector4<T>> dflux(Vector4<T>(0., 0., 0., 0.), U.size());
	std::valarray<Vector4<T>> lf = calcFlux(U, t, lam, nSize, opts...);

	const std::size_t ghost_point_number = 3;

	std::slice Nweno(ghost_point_number, nSize, 1);
	std::slice Nweno_shifted_by_neg_1(ghost_point_number - 1, nSize, 1);
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
	T t0, T dx, std::size_t n_size,
	T t_fin,
	auto&& timeStepFunction,
	T cfl = 0.4
) {
	/* Time Operator of the Riemann problem: perform the time loop
	 * and solve it, storing the result in `u` and the numerical flux
	 * needed for the calculation in `flux`.
	 */

	std::valarray<Vector4<T>> y2(Vector4<T>::ZERO,
		std::ranges::size(u));
	std::valarray<Vector4<T>> y3(Vector4<T>::ZERO,
		std::ranges::size(u));
	std::array<std::reference_wrapper<
		std::valarray<Vector4<T>>
	>, 2> fluxes = {
		std::ref(y2), std::ref(y3)
	};

	timeOperator<T>(
		u, flux, fluxes, t0, dx, n_size, t_fin,
		timeStepFunction,
		[](const decltype(u)& u, T dt) { return calcMaxWaveSpeedD<T>(u); },
		cfl
	);
}


template <ArithmeticWith<numeric_val> T>
std::function<Vector4<T>(Vector4<T>, T)> primitiveToConservativeU
	= primitiveToConservative<T>;


template <ArithmeticWith<numeric_val> T>
void prepareRiemannProblem(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	T gamma,
	T rho_left, T v_left, T p_left, T e_left, T rhoE_left,
	T rho_right, T v_right, T p_right, T e_rightR, T rhoE_right,
	T q0, T l_min, T l_max,
	std::function<Vector4<T>(Vector4<T>, T)> primitiveToConservativeU,
	std::size_t mesh_size
) {
	/* Fill u_init with the Riemann problem data
	 * at mesh_size number of nodes using Vector4<T>
	 * to store the data vector field. Fill x
	 * with the corresponding node coordinates.
	 */

	std::size_t computational_domain_size = mesh_size;
	const std::size_t n_ghost_points = 3;
	std::size_t full_mesh_size = computational_domain_size
		+ 2*n_ghost_points;
	T dx = (l_max - l_min) / (mesh_size-1);  // [L]

	x = std::valarray<T>(0., full_mesh_size);
	u_init = std::valarray<Vector4<T>>(
				Vector4<T>::ZERO,
				full_mesh_size);

	std::size_t k = 0;
	for (k = 0; k < n_ghost_points; ++ k)
		x[k] = l_min - dx * (n_ghost_points - k);

	x[n_ghost_points] = l_min;
	for (k = n_ghost_points + 1; k < full_mesh_size; ++ k)
		x[k] = x[k-1] + dx;

	std::size_t x0_index = 0;
	while (x[x0_index] < q0)
		++ x0_index;

	Vector4<T> vec(rho_left, v_left, p_left, 0.);
	vec = primitiveToConservativeU(vec, gamma);
	for (k = 0; k < x0_index; ++ k) {
		u_init[k] = vec;
	}

	vec = primitiveToConservativeU(
				Vector4(rho_right, v_right, p_right, 0.),
				gamma);
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
	std::function<Vector4<T>(Vector4<T>, T)> primitiveToConservativeU,
	auto&& updateGhostPoints,
	auto&& calcdSpace,
	std::size_t mesh_size, T cfl = 0.4
) {
	/* Solve a given Riemann problem for 1D Euler equations. */

	prepareRiemannProblem<T>(
		u_init, x, gamma,
		rho_left, v_left, p_left, e_left, rhoE_left,
		rho_right, v_right, p_right, e_rightR, rhoE_right,
		q0, l_min, l_max, primitiveToConservativeU, mesh_size
	);

	std::valarray<Vector4<T>> flux(Vector4<T>::ZERO,
		std::ranges::size(u_init));

	integrateRiemannProblem<T>(u_init, flux,
		t0, (l_max-l_min) / (mesh_size-1), mesh_size, t_max,
		[&calcdSpace, &updateGhostPoints](
			std::valarray<Vector4<T>>& u,
			std::valarray<Vector4<T>>& dflux,
			std::array<
				std::reference_wrapper<std::valarray<Vector4<T>>
			>, 2>& fluxes,
			T t, T dt, T dx,
			const std::valarray<T>& lam,
			std::size_t n_size
		) {
			advanceTimestepTVDRK3<T>(
				u, dflux, fluxes[0].get(), fluxes[1].get(),
				t, dt, dx, lam, n_size,
				updateGhostPoints, calcdSpace);
			/*EulerForward<T>(
				u, dflux,
				t, dt, dx, lam, n_size,
				updateGhostPoints, calcdSpace);*/
		}, cfl);

	return u_init;
}


template <typename T>
void addEmptySource(auto& u, auto& flux, auto&& x) {
	return;
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
	std::function<Vector4<T>(Vector4<T>, T)> primitiveToConservativeU,
	std::size_t mesh_size = 201, T cfl = 0.4
) {
	/* Solve a given Riemann problem for 1D Euler equations. */

	double e_left = p_left / (gamma - 1.) / rho_left;
	double e_right = p_right / (gamma - 1.) / rho_right;
	double E_left = (e_left + (v_left*v_left)/2.);
	double E_right = (e_right + (v_right*v_right)/2.);

	return solve1DRiemannProblemForEulerEq(
		u_init, x, gamma,
		rho_left, v_left, p_left, e_left, E_left,
		rho_right, v_right, p_right, e_right, E_right,
		q0, t0, t_max, l_min, l_max,
		primitiveToConservativeU,
		[gamma](
				std::valarray<Vector4<T>>& u,
				T t, T dx, const std::valarray<T>& max_eigenvalues,
				std::size_t n_size) {
			return calcdSpaceEu1D<T>(
				u, t, dx, max_eigenvalues, n_size,
				[gamma](
						const std::valarray<Vector4<T>>& u,
						T t, const std::valarray<T>& lam,
						std::size_t n_size/*,
						T eps = 1e-40, T p = 2.*/) {
//					return calcFluxComponentWiseWENO5<T>(
//						u, t, lam, n_size, eps, p);
					return calcFluxComponentWiseFVENO3<T>(
						u, t, lam, n_size, gamma);
				},
				[](
						const std::valarray<Vector4<T>>& u,
						std::valarray<Vector4<T>>& f,
						std::valarray<Vector4<T>>&& x = {}) {
					addEmptySource<T>(u, f, x);
				}/*,
				1e-40, 2.*/
			);
		},
		[](std::valarray<Vector4<T>>& u) {
			updateGhostPointsTransmissive<T>(u, 3);
		},
		mesh_size, cfl
	);
}

#endif // EULER1D_H
