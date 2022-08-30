#ifndef SSPRK10_4_H
#define SSPRK10_4_H

#include <algorithm>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimestepSSPRK10_4(
		std::ranges::common_range auto& u,
		std::ranges::common_range auto& dflux,
		std::ranges::common_range auto& u_temp,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_size,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* Optimal Low-Storage 4th-Order 10-Stage Explicit
	 * Total Variation Diming / Diminishing
	 * (Strong Stability Preserving)
	 * Runge-Kutta Scheme (TVD RK4 / SSPRK(10,4))
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 * This method has been found in Ketcheson (2008).
	 * It has an SSP coefficient of C = 6 (which is, as
	 * mentioned, provably optimal). C_{eff} = 0.60.
	 * See 'Highly Efficient Strong Stability-Preserving
	 * Runge–Kutta Methods with Low-Storage Implementations'
	 * by David I. Ketcheson (SIAM J. Sci. Comput. Vol. 30,
	 * No. 4, pp. 2113–2136).
	 *
	 * (It is important to note that this method is
	 * not only strong stability preserving but
	 * also internally stable.)
	 */

	// std::slice Nint(3, nSize, 1);
	dflux.resize(std::ranges::size(u));
	u_temp = u;

	// ------------------------First Five Stages------------------------
	// L[u] = (-) dF[u]/dx

	for (const int k [[maybe_unused]] : std::ranges::iota_view{0, 5}) {
		dflux = calcdSpace(u, t, dx, max_eigenvalues,
			n_size, opts_args...);  // L1 = L[u^n]
		std::transform(
					std::execution::par_unseq,
					std::ranges::begin(u), std::ranges::end(u),
					std::ranges::begin(dflux),
					std::ranges::begin(u),
					[dt](const auto u, const auto df) {
			return u + dt * df / 6.;
		});  // u(1) = u^n + Δt L[u^n] / 6.
		updateGhostPoints(u);
	}

	// ----------------------------Interim------------------------------
	std::transform(
		std::execution::par_unseq,
		std::ranges::begin(u), std::ranges::end(u),
		std::ranges::begin(u_temp),
		std::ranges::begin(u_temp),
		[dt](const auto q1, const auto q2) {
			return 0.04 * q2 + 0.36 * q1;
	});

	std::transform(
		std::execution::par_unseq,
		std::ranges::begin(u), std::ranges::end(u),
		std::ranges::begin(u_temp),
		std::ranges::begin(u),
		[dt](const auto q1, const auto q2) {
			return 15. * q2 - 5. * q1;
	});

	// ------------------------Next Four Stages-------------------------
	for (const int k [[maybe_unused]] : std::ranges::iota_view{5, 9}) {
		dflux = calcdSpace(u, t, dx, max_eigenvalues,
			n_size, opts_args...);  // L1 = L[u^n]
		std::transform(
					std::execution::par_unseq,
					std::ranges::begin(u), std::ranges::end(u),
					std::ranges::begin(dflux),
					std::ranges::begin(u),
					[dt](const auto u, const auto df) {
			return u + dt * df / 6.;
		});  // u(1) = u^n + Δt L[u^n] / 6.
		updateGhostPoints(u);
	}

	// ------------------------The Final Stage--------------------------
	std::transform(
		std::execution::par_unseq,
		std::ranges::begin(u), std::ranges::end(u),
		std::ranges::begin(dflux),
		std::ranges::begin(u),
		[dt](const auto q1, const auto df) {
			return 0.6 * q1 + 0.1 * dt * df;
	});

	std::transform(
		std::execution::par_unseq,
		std::ranges::begin(u), std::ranges::end(u),
		std::ranges::begin(u_temp),
		std::ranges::begin(u),
		[dt](const auto q1, const auto q2) {
			return q1 + q2;
	});

	updateGhostPoints(u);

	// return U;
	// U = res;
}

#endif // SSPRK10_4_H
