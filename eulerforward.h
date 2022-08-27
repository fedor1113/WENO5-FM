#ifndef EULERFORWARD_H
#define EULERFORWARD_H

#include <algorithm>
#include <cstddef>
#include <execution>
#include <ranges>

#include "arithmeticwith.h"


template<ArithmeticWith<numeric_val> T, typename... Args>
void EulerForward(
		std::ranges::common_range auto& U,
		std::ranges::common_range auto& dflux,
		T t, T dt, T dx,
		const std::ranges::common_range auto& max_eigenvalues,
		std::size_t n_size,
		auto&& calcdSpace,
		auto&& updateGhostPoints,
		Args... opts_args) {
	dflux = calcdSpace(U, t, dx, max_eigenvalues,
		n_size, opts_args...);  // L1 = L[u^n]
	// std::valarray<Vector4<T>> res(Vector4<T>::ZERO, U.size());
	// L[u] = (-) dF[u]/dx

	// Y2 = U + Vector4<T>(dt) * dflux;
	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(U), std::ranges::end(U),
				std::ranges::begin(dflux),
				std::ranges::begin(U),
				[dt](const auto u, const auto df) {
		return u + dt * df;
	});

	updateGhostPoints(U);
}

#endif // EULERFORWARD_H
