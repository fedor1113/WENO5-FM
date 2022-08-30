#include <algorithm>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

#include "arithmeticwith.h"
// #include "weno5coefs.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void advanceTimestepTVDRK3(
		std::ranges::common_range auto& U,
		std::ranges::common_range auto& dflux,
		std::ranges::common_range auto& Y2,
		std::ranges::common_range auto& Y3,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_size,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	/* Optimal 3rd Order 3 Stage Explicit Total Variation Diming
	 * / Diminishing (Strong Stability Preserving)
	 * Runge-Kutta Scheme (TVD RK3 / SSPRK(3,3))
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 * This is generally known as the Shu–Osher method.
	 * It has 3 stages. SSP coefficient = 1.
	 * See Shu and Osher (1988).
	 *
	 * (It is important to note that this method is
	 * not only strong stability preserving but
	 * also internally stable.)
	 */

	// std::slice Nint(3, nSize, 1);
	dflux.resize(std::ranges::size(U));

	// ------------------------First Stage----------------------------
	// std::valarray<Vector4<T>> flux = calcFlux(U, t, lam);
	dflux = calcdSpace(U, t, dx, max_eigenvalues,
		n_size, opts_args...);  // L1 = L[u^n]
	// std::valarray<Vector4<T>> res(Vector4<T>::ZERO, U.size());
	// L[u] = (-) dF[u]/dx

	// Y2 = U + Vector4<T>(dt) * dflux;
	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(U), std::ranges::end(U),
				std::ranges::begin(dflux),
				std::ranges::begin(Y2),
				[dt](const auto u, const auto df) {
		return u + dt * df;
	});
	// u(1) = u^n + Δt L[u^n]

	updateGhostPoints(Y2);


	// ------------------------Second Stage---------------------------
	dflux = calcdSpace(Y2, t, dx, max_eigenvalues,
		n_size, opts_args...);  // L2 = L[u(1)]

	// Y3 = Vector4<T>(3.)*U + Vector4<T>(dt) * dflux + Y2;
	// Y3 *= Vector4<T>(0.25);
	auto iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(std::ranges::size(U))
	);
	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv), std::ranges::end(iv),
				[dt, &Y3, &U, &dflux, &Y2](std::size_t k) {
		Y3[k] = (3. * U[k] + dt * dflux[k] + Y2[k]) * 0.25;
	});
	// u(2) = 0.75 * u^n + 0.25 * u(1) + 0.25 * Δt L[u(1)]

	updateGhostPoints(Y3);


	// ------------------------Third Stage----------------------------
	dflux = calcdSpace(Y3, t, dx, max_eigenvalues,
		n_size, opts_args...);  // L3 = L[u(2)]


	// U += Vector4<T>(2.) * (Y3 + Vector4<T>(dt) * dflux);
	// U *= Vector4<T>(1./3.);
	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(Y3), std::ranges::end(Y3),
				std::ranges::begin(dflux),
				std::ranges::begin(Y3),
				[dt](const auto u, const auto df) {
		return 2. * (u + dt * df);
	});
	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(U), std::ranges::end(U),
				std::ranges::begin(Y3),
				std::ranges::begin(U),
				[dt](const auto u, const auto y) {
		return (u + y) * (1./3.);
	});
	// u^(n+1) = (1/3) * u^n + (2/3) * u(2) + (2/3) * Δt L[u(2)]

	updateGhostPoints(U);

	// return U;
	// U = res;
}