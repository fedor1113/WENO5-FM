#include <algorithm>
#include <execution>
#include <iostream>
#include <span>
#include <ranges>
#include <valarray>
#include <vector>

#include "arithmeticwith.h"
#include "extrapolation.h"
#include "lf_flux.h"
#include "rk6_5.h"
// #include "ssprk33.h"
// #include "ssprk10_4.h"
// #include "tdrk3_5.h"
#include "weno5.h"


template <typename T>
T calcInitSpeedAheadOfDetu_a(T x_lab, T amp, T wn, short sign = +1) {
	/* Define the initial speed beyond the shock wave, u_a. */

	return amp / (1.-amp) * (1. + sign * std::sin(wn * x_lab));
}


template <typename T>
T calcDshock(
		T u_s,
		T amp,
		// std::optional<T> x_lab = std::nullopt,
		// std::optional<T> wn = std::nullopt,
		T x_lab,
		T wn,
		short sign = +1) {
	/* Define the speed of the shock wave, D,
	 * using the Rankine-Hugoniot condition for the Burgers eq'n.
	 */

	// if (x_lab)
	return (u_s + calcInitSpeedAheadOfDetu_a(x_lab, amp, wn, sign)) * 0.5;

	// if (!x_lab)
	// 	return (u_s + amp/(1.-amp)) * 0.5;
}


template <typename T>
T xiForKFR(T u_s, T alpha, T amp) {
	/* Define the induction function ξ dependence
	 * on the shock speed, u_s.
	 */

	return -std::pow((1. - amp) * u_s, -alpha);
}


template <typename T>
T reactiveBurgersSource(T x, T u_s, T alpha, T beta, T amp) {
	/* Define the source term S for the reactive Burgers
	 * (KFR) equation u_t + F_x = S.
	 */
	T ainv = (8.
			  * std::sqrt(std::numbers::pi_v<T> * beta)
			  * (1. + std::erf(1. / std::sqrt(4. * beta))));
	T source = std::exp(
		-std::pow(x - xiForKFR(u_s, alpha, amp),  2) / (4. * beta)
		) / ainv;

	return source;
}


template <typename T>
T BurgersFlux(T u) {
	/* Define the flux function F in the Burgers
	 * equation u_t + F_x = 0.
	 */

	return 0.5 * (u * u);
}


template <typename T>
T BurgersFluxDerivative(T u) {
	/* Return the derivative of the flux function F
	 * in the Burgers equation u_t + F_x = 0.
	 */

	return u;
}


template <typename T>
T reactiveBurgersFlux(
		T u,
		T u_s,
		T amp,
		T x_lab = 0.,
		T wn = 0.,
		short sign = +1) {
	/* Define the flux function F in the reactive Burgers
	 * (KFR) equation u_t + F_x = S.
	 */

	T D = calcDshock(u_s, amp, x_lab, wn, sign);

	return BurgersFlux(u) - D * u;
}


template <typename T>
T reactiveBurgersFluxDerivative(
		T u,
		T u_s,
		T amp,
		T x_lab = 0.,
		T wn = 0.,
		short sign = +1) {
	/* Define the flux function F in the reactive Burgers
	 * (KFR) equation u_t + F_x = S.
	 */

	T D = calcDshock(u_s, amp, x_lab, wn, sign);

	return BurgersFluxDerivative(u) - D;
}


template <typename T>
T reactiveBurgersFluxDerivative(T u, T D) {
	/* Define the flux function F in the reactive Burgers
	 * (KFR) equation u_t + F_x = S.
	 */

	return BurgersFluxDerivative(u) - D;
}


template <typename T>
T kfrSteadyStateSol(T x, T amp, T beta, T epsilon = 1e-6) {
	/* Calculate steady-state solution u0
	 * for the reactive Burgers (KFR) eq'n.
	 */

	T var1 = 1. + std::erf((x + 1.) / std::sqrt(4. * beta));
	T var2 = 1. + std::erf((0. + 1.) / std::sqrt(4. * beta));

	T u0 = (amp / (1. - amp)
			+ 0.5 * (1. + std::sqrt(var1 / var2)) + epsilon * std::sin(x));

	return u0;
}


// template <ArithmeticWith<numeric_val> T, std::size_t N>
// void polyfit(std::span<T> const, std::span<T> const, std::span<T, N>);

template
void polyfit<numeric_val, 5>(std::span<numeric_val> const argument_data,
			std::span<numeric_val> const function_value_data,
			std::span<numeric_val, 5 + 1> fitted_polynomial_coefficients,
			Eigen::Matrix<numeric_val, Eigen::Dynamic, 5 + 1>& vandermonde_mat);


Eigen::Matrix<numeric_val, Eigen::Dynamic, 5 + 1> VANDERMONDE_MAT(5 + 1, 5 + 1);


template <typename T>
void updateGhostPointsForDetonationProfile(
		std::span<T> u_interior,
		std::span<T> u_left,
		std::span<T> u_right,
		std::span<T> x,
		T dx) {
	/* Set outflow/transparent/transmissive boundary conditions
	 * in ghost points (fictitious cells) for u on (grid) x with
	 * interior (i. e. computational-domain) points inside u_init,
	 * ghost points on the left at u_left and on the right at u_right
	 * (assuming a perturbation on the right and constancy on the left).
	 *
	 * Waves will hopefully go through the boundaries
	 * as if the boundaries were not there.
	 *
	 * No waves travelling in from outside the domain can
	 * (presumably) influence the solution within it,
	 * so the ghost cells may be updated by
	 * extrapolating from the interior solution,
	 * assuming continuity at the boundary:
	 *
	 *                                         shock, x=0
	 *                  x=-L                     :
	 *    |     |     |  :  |                 |  :  |     |     |
	 *    |     |     |  :  |                 |  :  |     |     |
	 * o-----o-----o-----X-----X--  ...  --X-----X-----o-----o-----o
	 * [---u_left--]<--------------u_int-------------->[--u_right--]
	 *-0     1     2     3     4 interior N+2   N+3   N+4   N+5   N+6->

	 * The boundaries are shown as ':'.
	 * The physical domain [-L..0] is from x[3]=-L to x[N+3]=0.
	 * The shock is located at x[N+3] = 0.
	 * Internal points are from i=3 to i=N+3.

	 * The state right after the shock (at i=N+4:N+6)
	 * is given by parabolic extrapolation.

	 * Right before the left end of our domain -L
	 * we continue Ym for a bit assuming its approximate constancy there
	 * (constant / zero-order extrapolation). Mass, therefore, won't be
	 * conserved.
	 */

	// Transmissive b.c.s

	// --------------------------Left continuation----------------------
	std::ranges::for_each(u_left,
						  [&u_interior](auto& u_l) {
		u_l = u_interior.front();
	});
	// U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];
	// order 0, copying the value from the nearest non-ghost cell

	// At these ghost nodes a zero gradient condition is
	// enforced. Formally, this introduces spurious waves at the boundary.
	// However, as a check, the forward characteristic emanating from this
	// boundary was calculated, and it was guaranteed that the domain was
	// sufficiently large so as to prevent corruption of the shock and
	// reaction zone structure from this downstream acoustic noise.



	// --------------------------Right continuation---------------------
	const std::size_t order = 5;  // right boundary extrapolation order

//	std::ranges::for_each(u_right,
//						  [&u_interior](auto& u_r) {
//		u_r = u_interior.back();
//	});

	std::valarray<T> coefficients(order + 1);
	std::span<T, order + 1> coef_span = std::span<T, order + 1>{
		std::begin(coefficients), order + 1};

	std::span<T, order + 1> nextrap = x.subspan(
		x.size() - u_right.size() - 1 - order
	).template first<order + 1>();
	// Nextrap contains the points used for extrapolation
	// std::cout << nextrap.back() << "\n";

	polyfit<T, order>(nextrap, u_interior.last(order + 1), coef_span,
				VANDERMONDE_MAT);
	// Least squares polynomial extrapolation.

	// next we could use Horner's scheme here:
	std::ranges::transform(x.last(u_right.size()), std::begin(u_right),
						  [&coef_span](const auto& x) {
		T res = 0.;
		for (T coef : coef_span | std::views::reverse)
			res = coef + res * x;

		return res;
	});
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> calcFDSplitFlux(
		const std::ranges::common_range auto& u,
		std::ranges::common_range auto& u_flux,
		T t, const std::ranges::common_range auto& lam,
		std::size_t number_of_ghost_points = 3,
		T eps = 1e-40,
		T p = 2.) {


	auto [monotone_flux_component_pl,
			monotone_flux_component_mn] = splitFluxAsLaxFriedrichs(
				std::views::all(u), u_flux, lam[0]);

	calcHydroStageFDWENO5FM<T>(
		std::views::all(monotone_flux_component_pl),
		std::views::all(monotone_flux_component_mn),
		t, u_flux, number_of_ghost_points, eps, p);

//	const std::size_t u_size = std::ranges::size(u);
//	std::valarray<T> res(u_size);
//	T u_s = u[u_size - 1 - 3];

//	std::array<std::valarray<T>, 2> u_reconstruction {
//		std::valarray<T>(u_size),
//		std::valarray<T>(u_size)
//	};

//	calcHydroStageFVWENO5FM<T>(
//		std::views::all(u), t,
//		u_reconstruction[0], u_reconstruction[1],
//		3, eps, p);

//	res = calcLaxFriedrichsNumericalFlux<T>(
//				std::views::all(u_reconstruction[0]),
//				std::views::all(u_reconstruction[1]),
//				res, [u_s, amp, wn, &x0_lab](const auto& u_pt) {
//		return reactiveBurgersFlux(u_pt, u_s, amp, x0_lab, wn, +1);
//	}, lam[0]);

	return u_flux;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> calcFluxForKFR(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		T& x0_lab,
		T amp,
		T wn,
		T eps = 1e-40,
		T p = 2.) {
	const std::size_t u_size = std::ranges::size(u);
	std::valarray<T> res(u_size);
	T u_s = u[u_size - 1 - 3];
	std::transform(
		std::execution::par_unseq,
		std::ranges::begin(u),
		std::ranges::end(u),
		std::begin(res),
		[u_s, amp, wn, &x0_lab](const auto& u_pt) {
			return reactiveBurgersFlux(u_pt, u_s, amp, x0_lab, wn, +1);
	});

	return calcFDSplitFlux(u, res, t, lam, 3, eps, p);
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> calcFluxJacobianForKFR(
		const std::ranges::common_range auto& u,
		T t, const std::ranges::common_range auto& lam,
		T& x0_lab,
		T amp,
		T wn,
		T eps = 1e-40,
		T p = 2.) {
	const std::size_t u_size = std::ranges::size(u);
	std::valarray<T> res(u_size);
	T u_s = u[u_size - 1 - 3];
	std::transform(
		std::execution::par_unseq,
		std::ranges::begin(u),
		std::ranges::end(u),
		std::begin(res),
		[u_s, amp, wn, &x0_lab](const auto& u_pt) {
			return reactiveBurgersFluxDerivative(
						u_pt, u_s, amp, x0_lab, wn, +1);
	});

	return calcFDSplitFlux(u, res, t, lam, 3, eps, p);
}


template <ArithmeticWith<numeric_val> T>
T sourceTimeDerForKFR(T x, T u_s, T alpha, T beta, T amp) {
	T xi = xiForKFR(u_s, alpha, amp);
	T xinv = (x * (-1. / xi) + 1.);

	return 0.5 * alpha * xi * xi * xinv
			* std::pow(reactiveBurgersSource(
						   static_cast<T>(0.), u_s, alpha, beta, amp),
					   xinv * xinv) / beta / u_s;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> calcSpatialOperatorTimeDerivativeForKFR(
		const std::ranges::common_range auto& u,
		const std::ranges::common_range auto& du,
		const std::ranges::common_range auto& x,
		T t, const std::ranges::common_range auto& lam,
		T& x0_lab,
		T amp,
		T wn,
		T alpha,
		T beta,
		T eps = 1e-40,
		T p = 2.) {
	const std::size_t u_size = std::ranges::size(u);

	T u_s = u[u_size - 1 - 3];
	T du_s = du[u_size - 1 - 3];
	T D = calcDshock<T>(u_s, amp, std::ref(x0_lab), wn);
	T du_a_over_dt = D * amp * wn * std::cos(wn * x0_lab) / (1. - amp);
	std::valarray<T> res(u_size);

	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(u),
				std::ranges::end(u),
				std::begin(res),
				[u_s, amp, x0_lab, wn](auto u_pt) {
		return reactiveBurgersFluxDerivative(
					u_pt, u_s, amp, x0_lab, wn, +1);
	});
	res = calcFDSplitFlux(u, res, t, lam, 3, eps, p);
	// L_h [dF[u_h] / du]

	auto iv1 = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(std::ranges::size(u))
	);
	std::vector<std::size_t> iv(std::ranges::size(u));
	for (auto idx : iv1) iv[idx] = idx;

	std::for_each(
				std::execution::par_unseq,
				iv.begin(), iv.end(),
				[&](std::size_t k) {
		res[k] = res[k] * du[k]
				- u[k] * 0.5 * (du_s + du_a_over_dt)
				+ du_s * sourceTimeDerForKFR<T>(x[k], u_s, alpha, beta, amp);
	});

	return res;
}


template <ArithmeticWith<numeric_val> T, typename... Args>
std::valarray<T> calcdSpaceDet(
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

	// std::slice Nweno(ghost_point_number, n_size, 1);
	// std::slice Nweno_shifted_by_neg_1(ghost_point_number-1, n_size, 1);

	// std::valarray<T> f_mn = lf[Nweno_shifted_by_neg_1];
	// std::valarray<T> f_pl = lf[Nweno];
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

	// dflux[Nweno] = -(f_pl - f_mn) / dx;

	// dflux += source term:
	addSource(u, dflux, x);

	return dflux;
}


template <typename T>
void prepare1DDetonationProfileProblem(
		std::span<T> u_full,
		std::span<T> u_interior,
		std::span<T> u_left,
		std::span<T> u_right,
		std::span<T> x,
		T dx,
		// T q0,
		T l_min,  //= -10.,
		T l_max,  //= 0.,
		T alpha,
		T beta,
		T amp,
		T epsilon) {
	// ...reac

	std::size_t k = 0;
	for (k = 0; k < std::ranges::size(u_left); ++ k)
		x[k] = l_min - dx * (std::ranges::size(u_left) - k);

	x[std::ranges::size(u_left)] = l_min;
	for (k = std::ranges::size(u_left) + 1;
		 k < std::ranges::size(u_left)
			+ std::ranges::size(u_interior)
			+ std::ranges::size(u_right);
		 ++ k)
		x[k] = x[k-1] + dx;

	std::size_t x0_index = 0;
	T q0 = l_max;
	while (x[x0_index] < q0)
		++ x0_index;

	// std::size_t ishock = std::ranges::size(u_interior) + 3;
	// the shock state, i.e. the right boundary

	// ------------------Set the initial profile and BC------------------
	// initial profiles for steady solutions (could be other Cauchy data)
	T u_s = 1. / (1. - amp);
	// source = reactive_burgers_source(x, u_s, alpha, beta, amp);
	std::ranges::transform(std::as_const(x), std::begin(u_full),
						   [amp, beta, epsilon](const auto x) {
		return kfrSteadyStateSol(x, amp, beta, epsilon);
	});
	// steady-state solution u0

	updateGhostPointsForDetonationProfile(
		u_interior, u_left, u_right, x, dx
	);
}


template <ArithmeticWith<numeric_val> T>
void integrate1DDetonationProfileProblem(
	std::ranges::common_range auto& u_sol,
	std::ranges::common_range auto& flux,
	std::ranges::common_range auto& u_s_sol,
	std::ranges::common_range auto& times,
	T t0, T dx, std::size_t n_size, T t_fin,
	T amp, T wn,
	auto&& timeStepFunction,
	T cfl = 0.4
) {
	/* Time Operator of the detonation problem: perform the time loop
	 * and solve it, storing the result in `u_sol` and the numerical
	 * flux needed for the calculation in `flux`.
	 */

	std::valarray<T> y1(0., std::ranges::size(u_sol));
	std::valarray<T> y2(0., std::ranges::size(u_sol));
	std::valarray<T> y3(0., std::ranges::size(u_sol));
	std::valarray<T> y4(0., std::ranges::size(u_sol));
	std::valarray<T> y5(0., std::ranges::size(u_sol));
//	std::array<std::reference_wrapper<std::valarray<T>>, 2> fluxes = {
//		std::ref(y2), std::ref(y3)
//	};
	std::valarray<T> dy1(0., std::ranges::size(u_sol));
	std::valarray<T> dy2(0., std::ranges::size(u_sol));
	std::valarray<T> dy3(0., std::ranges::size(u_sol));
	std::valarray<T> dy4(0., std::ranges::size(u_sol));
	std::valarray<T> dy5(0., std::ranges::size(u_sol));

	std::array<std::reference_wrapper<std::valarray<T>>, 10> fluxes = {
		std::ref(y1), std::ref(dy1),
		std::ref(y2), std::ref(dy2),
		std::ref(y3), std::ref(dy3),
		std::ref(y4), std::ref(dy4),
		std::ref(y5), std::ref(dy5)
	};

	T x0_lab = 0;
	std::size_t time_index = 0;

	timeOperator<T>(
		u_sol, flux, fluxes, t0, dx, n_size, t_fin,
		timeStepFunction,
		[amp, wn, &x0_lab,
		&time_index,
		&u_sol, &u_s_sol, &times, t0](
				const decltype(u_sol)& u, T dt) -> T {
			T u_s = u[std::ranges::size(u)-1-3];
			T D = calcDshock<T>(
							u[std::ranges::size(u)-1-3],
							amp, std::ref(x0_lab), wn);

			auto a0 = std::ranges::views::transform(
					[D](const auto& u_pt) -> T {
				// return (u_pt - D);
				return reactiveBurgersFluxDerivative(u_pt, D);
			});

			x0_lab += D * dt;

			if (std::ranges::size(u_s_sol) <= time_index) {
				u_s_sol.resize(2 * std::ranges::size(u_s_sol));
				times.resize(2 * std::ranges::size(times));
			}

			u_s_sol[time_index] = u_s;
			if (time_index > 0)
				times[time_index] = times[time_index - 1] + dt;
			else
			 	times[0] = t0;
			++ time_index;

			return std::ranges::max(std::as_const(u) | a0);
		},
		cfl,
		std::ref(x0_lab)
	);

	u_s_sol.resize(time_index);
	times.resize(time_index);
}


template <typename T>
std::valarray<T> solve1DDetonationProfileProblem(
	std::ranges::common_range auto& u_init,
	std::ranges::common_range auto& x,
	T alpha,
	T beta,
	T amp,
	T wn,
	std::ranges::common_range auto& u_s_sol,
	std::ranges::common_range auto& times,
	T epsilon=1e-6,
	// T q0, // Initial coordinate of the discontinuity
	T t0 = 0., T t_max = 3000, T l_min = -10., T l_max = 0.,
	std::size_t mesh_size = 201, T cfl = 0.4
) {
	T t = t0;
	T dx = (l_max - l_min) / (mesh_size/* - 1*/);
	const std::size_t number_of_ghost_points = 3;

	prepare1DDetonationProfileProblem<T>(
		std::span<T>{u_init},
		std::span<T>{std::begin(u_init)+number_of_ghost_points,
			std::end(u_init)-number_of_ghost_points},
		std::span<T>{std::begin(u_init), number_of_ghost_points},
		std::span<T>{std::end(u_init)-3, number_of_ghost_points},
		std::span<T>{std::begin(x), std::end(x)},
		dx, /* T q0, */ l_min, l_max,
		alpha, beta, amp, epsilon);

//	std::array<std::valarray<T>, 2> monotone_flux_components {
//		std::valarray<T>(std::ranges::size(u_init)),
//		std::valarray<T>(std::ranges::size(u_init))
//	};

	auto spaceOp = [alpha, beta, amp, wn, &x/*, &monotone_flux_components*/](
			std::span<T> const u,
			T t, T dx, const std::valarray<T>& max_eigenvalues,
			T n_size, T x0_lab) {
		return calcdSpaceDet<T>(
				u, std::span{x},
				t, dx, max_eigenvalues, n_size,
				[amp, wn](
						std::span<T> u,
						T t, const std::valarray<T>& lam,
						std::size_t n_size, T x0_lab,
						T eps = 1e-40, T p = 2.) {
					return calcFluxForKFR<T>(
						u, t, lam, x0_lab, amp, wn, eps, p);
				},
				[alpha, beta, amp](
						std::span<T> const u,
						std::valarray<T>& f,
						std::span<T> x) {
					std::transform(
							std::execution::par_unseq,
							std::begin(x), std::end(x),
							std::begin(f), std::begin(f),
							[&u, alpha, beta, amp](
									const auto x_el,
									const auto f_el) {
						T u_s = u[std::ranges::size(u) - 1 - 3];
						return f_el + reactiveBurgersSource<T>(
								x_el, u_s, alpha, beta, amp);
					});
				},
				x0_lab, 1e-40, 2.
		);
	};

//	std::valarray<T> ddflux(0., std::ranges::size(u_init));

//	auto secondDerivative = [alpha, beta, amp, wn, &x](
//			std::span<T> const u, std::span<T> const du,
//			T t, T dx, const std::valarray<T>& max_eigenvalues,
//			T n_size, T x0_lab) {
//		return calcSpatialOperatorTimeDerivativeForKFR<T>(
//			u, du, x,
//			t, max_eigenvalues, x0_lab,
//			amp, wn, alpha, beta);
//	};

	auto updateGhostPoints = [&x, number_of_ghost_points, dx](
			std::valarray<T>& u) {
		updateGhostPointsForDetonationProfile<T>(
			std::span<T>{std::begin(u)+number_of_ghost_points,
				std::end(u)-number_of_ghost_points},
			std::span<T>{std::begin(u), number_of_ghost_points},
			std::span<T>{std::end(u)-3, number_of_ghost_points},
			std::span<T>{std::begin(x), std::end(x)},
			dx
		);
	};

	std::valarray<T> flux(0., std::ranges::size(u_init));

	integrate1DDetonationProfileProblem<T>(u_init, flux,
		u_s_sol, times,
		t0, dx, mesh_size, t_max, amp, wn,
		[&spaceOp, /*&secondDerivative, */&updateGhostPoints/*, &ddflux*/](
//		[&spaceOp, &updateGhostPoints](
			std::valarray<T>& u,
			std::valarray<T>& dflux,
//			std::array<
//				std::reference_wrapper<std::valarray<T>
//			>, 2> fluxes,
			std::array<
				std::reference_wrapper<std::valarray<T>
			>, 10> fluxes,
			T t, T dt, T dx,
			const std::valarray<T>& lam,
			std::size_t n_size,
			T& x0_lab
		) {
//			advanceTimestepTVDRK3<T>(
//				u, dflux, fluxes[0].get(), fluxes[1].get(),
//				t, dt, dx, lam,
//				3, spaceOp, updateGhostPoints, x0_lab);
//			advanceTimestepSSPRK10_4<T>(
//				u, dflux, fluxes[0].get(),
//				t, dt, dx, lam,
//				n_size, spaceOp, updateGhostPoints, x0_lab);
//			advanceTimestepTDRK3_5<T>(
//				u, dflux, ddflux,
//				fluxes[0].get(), fluxes[1].get(),
//				fluxes[2].get(), fluxes[3].get(),
//				t, dt, dx, lam,
//				n_size, spaceOp, secondDerivative,
//				updateGhostPoints, x0_lab);
			advanceTimestepRK6_5<T>(
				u,
				dflux,
				fluxes[0].get(), fluxes[1].get(),
				fluxes[2].get(), fluxes[3].get(),
				fluxes[4].get(), fluxes[5].get(),
				fluxes[6].get(), fluxes[7].get(),
				fluxes[8].get(), fluxes[9].get(),
				t, dt, dx, lam,
				3, spaceOp, updateGhostPoints, x0_lab);
		}, cfl);

	return u_init;
}
