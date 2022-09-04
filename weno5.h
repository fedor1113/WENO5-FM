#ifndef WENO5_H
#define WENO5_H

// #include "weno5coefs.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <execution>
// #include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <ranges>
// #include <span>
// #include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

#include "arithmeticwith.h"

// template class Vector4<double>;

// using Vec4 = Vector4<double>;


template <ArithmeticWith<numeric_val> T, typename T1, typename T2>
std::valarray<T> vecMatDot(const T1& vec, const T2& mat) {
	/* Multiply a vector by a matrix (a sum product over
	 * the last axis of mat and vec).
	 */

	std::valarray<T> result(std::ranges::size(mat));

	for (std::size_t k = 0; k < std::ranges::size(mat)
			&& k < std::ranges::size(vec); ++ k)
		result[k] = std::inner_product(std::ranges::begin(mat[k]),
									   std::ranges::end(mat[k]),
									   std::ranges::begin(vec), 0.);

	return result;
}




template <ArithmeticWith<numeric_val> T>
//std::array<T, 3> smoothness_indicators(const T1& f_stencil) {
std::valarray<T> betaSmoothnessIndicators(
		const std::ranges::sized_range auto& f_stencil) {
	/* Return the WENO5 smoothness indicators of Jiang and Shu (1996)
	 * for each of the 3 substencils.
	 * That is the sum of the normalized squares of the scaled
	 * L2-norms of all the derivatives of 3 local interpolating
	 * polynomials in the sub-stencils of 5-node `f_stencil`.
	 *
	 * This allows for (2*3-1)=5th order accuracy from the 3rd
	 * order Eno schemes.
	 */

	// std::array<T, 3> res;
	std::valarray<T> betas(3);

	betas[0] = ((13./12.) * std::pow(
					f_stencil[0] - 2.*f_stencil[1] + f_stencil[2], 2)
				+ (1./4.) * std::pow(
					f_stencil[0] - 4.*f_stencil[1] + 3.*f_stencil[2], 2
					));

	betas[1] = ((13./12.) * std::pow(
					f_stencil[1] - 2.*f_stencil[2] + f_stencil[3], 2)
				+ (1./4.) * std::pow((f_stencil[1] - f_stencil[3]), 2));

	betas[2] = ((13./12.) * std::pow(
					f_stencil[2] - 2.*f_stencil[3] + f_stencil[4], 2)
				+ (1./4.) * std::pow(
					3.*f_stencil[2] - 4.*f_stencil[3] + f_stencil[4], 2
					));

	return betas;
}


//template <ArithmeticWith<numeric_val> T>
//std::valarray<T> betaSmoothnessIndicatorsMat(
//		const std::ranges::common_range auto& f_stencil,
//		std::array<std::array<const T, 6>, 6> _coefs[] = plus_coefs) {
//	/* Return the smoothness indicators beta_k, k=0,1,2
//	 * for each of the 3 substencils of `f_stencil`.
//	 */

//	std::valarray<T> res(3);

//	for (std::size_t k = 0; k < 3; ++ k) {
//		// for (std::size_t k = 0; k < half_size + 1; ++ k)
//		res[k] = std::inner_product(
//			std::ranges::begin(f_stencil),
//			std::ranges::end(f_stencil),
//			std::ranges::begin(vecMatDot(f_stencil, _coefs[k])),
//			0.
//		);
//	}


//	return res;
//}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> f3OrdReconstructionFromStencil(
		const std::ranges::sized_range auto& f_stencil) {
	/* 3rd order reconstructions of f(j) from all the 3 3-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *               (or [j+3, j+2, j+1, j+0, j-1, ...] for '-').
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 */
	//		{{2./6, -7./6, 11./6, 0., 0., 0.}},
	//		{{0., -1./6, 5./6, 2./6, 0., 0.}},
	//		{{0., 0., 2./6, 5./6, -1./6, 0.}}

	std::valarray<T> q_res(3);

	q_res[0] = (2. * f_stencil[0]
			  - 7. * f_stencil[1]
			 + 11. * f_stencil[2]) / 6.;

	q_res[1] = (-1. * f_stencil[1]
			   + 5. * f_stencil[2]
			   + 2. * f_stencil[3]) / 6.;

	q_res[2] = (2. * f_stencil[2]
			  + 5. * f_stencil[3]
			  - 1. * f_stencil[4]) / 6.;

	return q_res;
}


template <ArithmeticWith<numeric_val> T>
T henrickGMappingForLambda(T lambda_weno_weight,
						   T lambda_ideal = 1./3.) {
	/* The mapping function g by Henrick modified for symmetric
	 * lambda-weights by Zheng Hong, Zhengyin Ye and Kun Ye.
	 */

	T square_ideal = lambda_ideal * lambda_ideal;

	return lambda_weno_weight
			* (lambda_ideal
			   + square_ideal
			   - 3. * lambda_ideal * lambda_weno_weight
			   + lambda_weno_weight * lambda_weno_weight)
			/ (square_ideal
			   + lambda_weno_weight
			   * (1. - 2.*square_ideal));
}


//template <ArithmeticWith<numeric_val> T>
//T henrickGMapping(T omega_weno_weight,
//				  T d_ideal = 1./3.) {
//	/* The original mapping function g by Henrick et al. */

//	T d_square = d_ideal * d_ideal;
//	return omega_weno_weight
//			* (d_ideal
//			   + d_square
//			   - 3. * d_ideal * omega_weno_weight
//			   + omega_weno_weight * omega_weno_weight)
//			/ (d_square
//			   + omega_weno_weight
//			   * (1. - 2.*d_ideal));
//}


// const std::ranges::common_range auto&&
template <ArithmeticWith<numeric_val> T>
T alphaWENO5FMWeight(
		T beta_IS_coefficient,
		T epsilon = 1e-40,
		T p = 2.) {
	/* Compute appropriate alpha(α)-weights for the WENO5-FM scheme,
	 * by which the inverses of smoothness indicators are meant,
	 * so inverse beta(β) with the caveat of aritificially finite
	 * answers using the added epsilon-parameter to the beta-weights.
	 *
	 * `p` controls (increases) the amount of numerical dissipation
	 * (it's recommended to take it = r-1 for 2r-1 order schemes,
	 * so 2 for WENO5).
	 */

	return 1. / std::pow(epsilon + beta_IS_coefficient, p);
}


template <ArithmeticWith<numeric_val> T>
void lambdaWENO5FMWeights(
		const std::ranges::common_range auto&& alpha_weights,
		std::ranges::common_range auto&& res) {
	/* FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
	 * for WENO5-FM or WENO5-ZM
	 * due to Zheng Hong, Zhengyin Ye and Kun Ye:
	 * lambda_weights = alpha_weights / alpha_weights.sum();
	 */

	// return alpha_weights / alpha_weights.sum();
	T sum = std::reduce(
				/*std::execution::par_unseq,*/
				std::ranges::begin(alpha_weights),
				std::ranges::end(alpha_weights),
				0.);

	std::transform(
				/*std::execution::par_unseq,*/
				std::ranges::begin(alpha_weights),
				std::ranges::end(alpha_weights),
				std::ranges::begin(res),
				[sum](const auto alpha) {
		return alpha / sum;
	});
}


template <ArithmeticWith<numeric_val> T>
std::ranges::common_range auto omegaWENO5FMWeights(
		const std::ranges::common_range auto&& lambda_weights) {
	/* From Henrick et al.'s mappings of g(λ_k) for the improved
	 * symmetric normalized lambda-weights of Hong, Ye & Ye
	 * and linear weights d_k we get the new corrected resultant
	 * normalized WENO5-FM (WENO5-ZM) omega (ω_k-)weights for WENO5-FM
	 * (again due to Zheng Hong, Zhengyin Ye and Kun Ye).
	 */

	// The ideal weights (they generate the central upstream fifth-order
	// scheme for the 5-point stencil), which are in WENO usu. called
	// linear weights:
	std::valarray<T> d_lin_weights = {0.1, 0.6, 0.3};
	// From them WENO5-Z and WENO-M will calculate the non-linear
	// alpha and omega weights.

	// In WENO5-FM, further, we have one ideal value for λ
	// \overbar{lambda} = 1/3
	// T lambda_ideal = 1/3;
	// In the smooth region the smoothness indicators β_k ought to
	// be equal for all sub-stencils, and thus the weight's ideal
	// value must be unique.

	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
	// omega_weights = d_lin_weights * alpha_weights;
	// lambda_weights

	// And only to the λ-weights a mapping in the spirit of
	// Henrick et al. is applied:
	auto gMap = [](T x) -> T {
		return henrickGMappingForLambda(x);
	};

	std::valarray<T> alpha_weights(3);
	std::ranges::transform(
				lambda_weights,
				std::ranges::begin(alpha_weights),
				gMap);
	// α*-weights

	// From α*=g(λ_k) and d_k we get the new corrected resultant
	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
	// omega_weights = d_lin_weights * alpha_weights;
	std::valarray<T> omega_weights = d_lin_weights * alpha_weights;
	omega_weights /= omega_weights.sum();

	return omega_weights;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5JSReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);
//	beta_IS_coefs[0] = (13.0/12.0)*std::pow(
//				f_stencil[0]-2.0*f_stencil[1]+f_stencil[2], 2)
//			+ 0.25*std::pow(
//				f_stencil[0]-4.0*f_stencil[1]+3.0*f_stencil[2], 2);
//	beta_IS_coefs[1] = (13.0/12.0)*std::pow(
//				f_stencil[1]-2.0*f_stencil[2]+f_stencil[3], 2)
//			+ 0.25*std::pow(f_stencil[1]-f_stencil[3], 2);
//	beta_IS_coefs[2] = (13.0/12.0)*std::pow(
//				f_stencil[2]-2.0*f_stencil[3]+f_stencil[4], 2)
//			+ 0.25*std::pow(
//				3.0*f_stencil[2]-4.0*f_stencil[3]+f_stencil[4], 2);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {0.1, 0.6, 0.3};
	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
//	std::valarray<T> alpha_weights(3);
//	alpha_weights[0] = 1.0e-1/std::pow((eps+beta_IS_coefs[0]), 2);
//	alpha_weights[1] = 6.0e-1/std::pow((eps+beta_IS_coefs[1]), 2);
//	alpha_weights[2] = 3.0e-1/std::pow((eps+beta_IS_coefs[2]), 2);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();
//	std::valarray<T> omega_weights(3);
//	omega_weights[0] = alpha_weights[0]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);
//	omega_weights[1] = alpha_weights[1]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);
//	omega_weights[2] = alpha_weights[2]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<T> eno_reconstructed_f
			= f3OrdReconstructionFromStencil<T>(f_stencil);

//	 eno_reconstructed_f[0] = f_stencil[0]/3.0
//			 - 7.0/6.0*f_stencil[1] + 11.0/6.0*f_stencil[2];
//	 eno_reconstructed_f[1] =-f_stencil[1]/6.0
//			 + 5.0/6.0*f_stencil[2] + f_stencil[3]/3.0;
//	 eno_reconstructed_f[2] = f_stencil[2]/3.0
//			 + 5.0/6.0*f_stencil[3] - f_stencil[4]/6.0;

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5JSReconstructionKernelRev(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {0.3, 0.6, 0.1};
	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights / alpha_weights.sum();

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of f_{i+1/2}
	// via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<T> eno_reconstructed_f(3);
	eno_reconstructed_f[0] = (-1. * f_stencil[0]
						   + 5. * f_stencil[1]
						   + 2. * f_stencil[2]) / 6.;
	eno_reconstructed_f[1] = (2. * f_stencil[1]
						   + 5. * f_stencil[2]
						   - 1. * f_stencil[3]) / 6.;
	eno_reconstructed_f[2] = (11. * f_stencil[2]
						   - 7. * f_stencil[3]
						   + 2. * f_stencil[4]) / 6.;

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];


	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5MReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {0.1, 0.6, 0.3};
	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::transform(
				/*std::execution::par_unseq,*/
				std::ranges::begin(omega_weights),
				std::ranges::end(omega_weights),
				std::ranges::begin(d_lin_weights),
				std::ranges::begin(omega_weights),
				[](auto w, auto d) {
		return henrickGMappingForLambda(w, d);
	});  // we obtain a new non-normalized weight alpha*

	omega_weights = omega_weights / omega_weights.sum(); // normalize it

	std::valarray<T> eno_reconstructed_f
			= f3OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5FMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);

	std::array<T, 3> alpha_weights;
	std::ranges::transform(
				beta_IS_coefs,
				std::ranges::begin(alpha_weights),
				[eps, p](auto beta) {
		return alphaWENO5FMWeight(beta, eps, p);
	});

	std::array<T, 3> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	std::valarray<T> omega_weights = omegaWENO5FMWeights<T>(
				std::move(lambda_weights));

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<T> eno_reconstructed_f
			= f3OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


//template <ArithmeticWith<numeric_val> T>
//T computeFHatWENO5FMReconstructionKernelRev(std::span<T, 5> f_stencil,
//											T eps = 1e-40, T p = 2.) {
//	/* Calculate (reconstruct) one of the two split monotone numerical
//	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
//	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
//	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
//	 *                      ^    ^    ^    ^    ^    ^
//	 *                      0    1    2    3    4    |
//	 * in either case for convenience).
//	 *
//	 * I.e. this function implements the downwind reconstruction which
//	 * should be used for negative fluxes (with information propagated
//	 * from right to left) if the nodes are passed in order.
//	 */

//	std::valarray<T> beta_IS_coefs(3);

//	T f_hat = 0.;

//	// smoothness indicators of the stencil
//	// (measure how smooth u is in the stencil)
////	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

//	// The non-matrix variant seems to be faster(?)
//	// beta_IS_coefs = betaSmoothnessIndicators(f_stencil);

//	T f_prev2 = f_stencil[0];
//	T f_prev1 = f_stencil[1];
//	T f_curr0 = f_stencil[2];
//	T f_next1 = f_stencil[3];
//	T f_next2 = f_stencil[4];
////	T f_prev2 = f_stencil[4];
////	T f_prev1 = f_stencil[3];
////	T f_curr0 = f_stencil[2];
////	T f_next1 = f_stencil[1];
////	T f_next2 = f_stencil[0];

//	beta_IS_coefs[0] = ((13./12.) * std::pow(
//							f_prev2 - 2.*f_prev1 + f_curr0, 2)
//				+ (1./4.) * std::pow(
//							f_prev2 - 4.*f_prev1 + 3.*f_curr0, 2));

//	beta_IS_coefs[1] = ((13./12.) * std::pow(
//							f_prev1 - 2.*f_curr0 + f_next1, 2)
//				+ (1./4.) * std::pow((f_/usr/bin/prev1 - f_next1), 2));

//	beta_IS_coefs[2] = ((13./12.) * std::pow(
//							f_curr0 - 2.*f_next1 + f_next2, 2)
//				+ (1./4.) * std::pow(
//							3.*f_curr0 - 4.*f_next1 + f_next2, 2));

//	// non-linear non-scaled (α-)weights
//	std::valarray<T> d_lin_weights = {0.3, 0.6, 0.1};
//	std::valarray<T> alpha_weights = d_lin_weights
//			/ std::pow(eps + beta_IS_coefs, p);
//	// — we have no need of them in WENO-FM!

//	// Instead we use:
////	std::valarray<T> alpha_weights = alphaWENO5FMWeights(
////		std::move(beta_IS_coefs), eps, p
////	);

//	// scaled (normalized) non-linear (ω-)weights (ENO weights)
//	 std::valarray<T> omega_weights = alpha_weights
//			 / alpha_weights.sum();

//	// FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
//	// due to Zheng Hong, Zhengyin Ye and Kun Ye:
////	 lambda_weights = alpha_weights / alpha_weights.sum();
////	std::valarray<T> lambda_weights = lambdaWENO5FMWeights(
////		std::move(alpha_weights)
////	);

////	std::valarray<T> omega_weights = omegaWENO5FMWeights(
////		std::move(lambda_weights)
////	);

//	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
//	// f_{i+1/2} via linear combinations with WmNplus coefficients
//	// for each substencil which is then used to calculate
//	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
//	// using the nonlinear weights [ω]
////	f_hat = std::inner_product(
////		std::ranges::begin(omega_weights),
////		std::ranges::end(omega_weights),
////		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)),
////		0.
////	);

//	 std::valarray<
//		 T
//	 > eno_reconstructed_f(3);
//	 eno_reconstructed_f[0] = (-1. * f_stencil[0]
//							   + 5. * f_stencil[1]
//							   + 2. * f_stencil[2]) / 6.;

//	 eno_reconstructed_f[1] = (2. * f_stencil[1]
//							   + 5. * f_stencil[2]
//							   - 1. * f_stencil[3]) / 6.;

//	 eno_reconstructed_f[2] = (11. * f_stencil[2]
//							   - 7. * f_stencil[3]
//							   + 2. * f_stencil[4]) / 6.;

//	f_hat = omega_weights[0] * eno_reconstructed_f[0]
//			+ omega_weights[1] * eno_reconstructed_f[1]
//			+ omega_weights[2] * eno_reconstructed_f[2];

//	return f_hat;
//}


// FD WENO5FM (WENO5-FM) - method
// Reconstruction based on LF flux splitting + improved mapped WENO
// of 5th order
// (see Mapped weighted essentially non-oscillatory schemes:
// achieving optimal order near critical points, 2005 by Henrick et al.)
// and 'An improved WENO-Z scheme with symmetry-preserving mapping'
// by Zheng Hong, Zhengyin Ye and Kun Ye, 2020
template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO5FM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	/* Component-wise finite-difference WENO5FM (FD WENO5-FM) - space
	 * reconstruction method with the global Lax-Friedrichs (LF) flux
	 * splitting.
	 *
	 * Usually, componentwise reconstruction produces satisfactory
	 * results for schemes up to third-order accuracy, while characteristic
	 * reconstruction produces better nonoscillatory results for
	 * higher-order accuracy, albeit with an increased computational cost.
	 */

	const unsigned order = 5;
	const std::size_t stencil_size = order;
	const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = order / 2;  // 2

	// r = (order + 1) / 2 = 3
	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(
				numerical_flux) - n_ghost_cells - 1;
	auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	// [g      g      g      i      i      i      i      i      i      ...]
	// {0      1      2      3      4      5}     6      7      8      ...
	//  |             |      |      |
	// itr            j nGhostCells end()

	// WENO5 stencils

	// Coefficients WmN(+/-) before fluxes at the stencil nodes
	// to find component stencils
	// [q_k] = WmN(+/-) * [ f[j-2] ... f[j+2] ]
	// of the WENO interpolator.
	// So 3rd order approximation coefficients associated with
	// the substencils.

	// Calculation of f_hat, the numerical flux of u (whichever
	// is chosen), requires the approximation of u that uses at
	// the j-th cell of u-discretization a group of cell average
	// values on the left (`u_minus`) and on the right (`u_plus`).
	// So left- and right-biased approximations respectively.
	// `u_plus` represents the cells [j-2, j-1, j, j+1, j+2],
	// and `u_minus` represents the cells [j-1, j, j+1, j+2, j+3];
	// for convenience and uniformity we represent both using the
	// same combined structure of    [j-2, j-1, j, j+1, j+2, j+3].
	// std::valarray<T> u_plus(_actual_stencil_size);   // f_plus
	// std::valarray<T> u_minus(_actual_stencil_size);  // f_minus

	// For the purpose of linear stability (upwinding),
	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	// dfminus/du <= 0), is performed.
	// Lax-Friedrichs (LF) flux splitting (because it is the simplest
	// and smooth - we need the positive and negative fluxes to have
	// as many derivatives as the order of our finite-difference WENO)
	// is chosen here. `f_plus` uses a biased stencil with 1 point to
	// the left, while `f_minus` uses a biased stencil with 1 point to
	// the right.

	T fhatminus = 0.;
	T fhatplus = 0.;

	// So an LF flux	`numerical_flux`, f_hat(u_minus, u_plus),
	// a monotone numerical flux consistent with the physical one
	// (f_hat(u, u) = f(u)), will be construced at the end of
	// the loop below from `fhatminus` and `fhatplus` using
	// `f_minus` and `f_plus`.
	// N.B.! This can be replaced by an exact or approximate
	// Riemann solver (see Toro, 2009). Somehow...
	// Not every monotone flux can be writtenin the flux split form.
	// For example, the Godunov flux cannot.
	// std::valarray<T> numerical_flux(0., u.size());
	// f = std::valarray<T>(u.size());

	auto j_it_p = std::ranges::begin(f_plus);  // f_plus
	auto j_it_m = std::ranges::begin(f_minus);  // f_minus

//	std::ranges::transform(
//			monotone_flux_components[0],
//			monotone_flux_components[1],
//			std::ranges::begin(f), [](const auto fp, const auto fm) {

//	})
	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size);
	std::advance(j_it_m, mini - 1 + half_size + 1 - stencil_size);
	auto u_plus = std::ranges::views::counted(j_it_p, 6);
	auto u_minus = std::ranges::views::counted(j_it_m, 6)
					| std::ranges::views::reverse;

	for (std::size_t j : shifted_index_range) {
		j_it_p = std::ranges::begin(f_plus);  // f_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size);
		u_plus = std::ranges::views::counted(j_it_p, 6);

		j_it_m = std::ranges::begin(f_minus);  // f_minus
		std::advance(j_it_m, j + half_size + 1 - stencil_size);
		u_minus = std::ranges::views::counted(j_it_m, 6)
					| std::ranges::views::reverse;

		fhatplus = computeFHatWENO5FMReconstructionKernel<T>(
			std::ranges::views::counted(
						std::ranges::begin(u_plus), 5), eps, p
		);

		fhatminus = computeFHatWENO5FMReconstructionKernel<T>(
			std::ranges::subrange(
				std::ranges::begin(u_minus),
						std::ranges::end(u_minus) - 1), eps, p
		);
//		fhatminus = computeFHatWENO5JSReconstructionKernelRev<T>(
//			std::ranges::views::counted(
//						std::ranges::begin(u_minus)+1, 5), eps, p
//		);

		numerical_flux[j] = fhatplus + fhatminus;
	}

	// return numerical_flux;
	// std::cout << " done!" << "\n";
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFVWENO5FM(
		const std::ranges::common_range auto&& u,
		T t,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	/* Component-wise finite-volume WENO5FM (FV WENO5-FM) - space
	 * reconstruction method with the global Lax-Friedrichs (LF) flux
	 * splitting.
	 *
	 * Usually, componentwise reconstruction produces satisfactory
	 * results for schemes up to third-order accuracy, while characteristic
	 * reconstruction produces better nonoscillatory results for
	 * higher-order accuracy, albeit with an increased computational cost.
	 */

	const unsigned order = 5;
	const std::size_t stencil_size = order;
	// const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = order / 2;  // 2

	// r = (order + 1) / 2 = 3
	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;  // 3
	const std::size_t mini = n_ghost_cells;  // at least 3
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(u) - n_ghost_cells - 1;
	auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	// [g      g      g      i      i      i      i      i      i      ...]
	// {0      1      2      3      4      5}     6      7      8      ...
	//  |             |      |      |
	// itr            j nGhostCells end()

	// WENO5 stencils

	// Coefficients WmN(+/-) before fluxes at the stencil nodes
	// to find component stencils
	// [q_k] = WmN(+/-) * [ f[j-2] ... f[j+2] ]
	// of the WENO interpolator.
	// So 3rd order approximation coefficients associated with
	// the substencils.

	// Calculation of f_hat, the numerical flux of u (whichever
	// is chosen), requires the approximation of u that uses at
	// the j-th cell of u-discretization a group of cell average
	// values on the left (`u_minus`) and on the right (`u_plus`).
	// So left- and right-biased approximations respectively.
	// `u_plus` represents the cells [j-2, j-1, j, j+1, j+2],
	// and `u_minus` represents the cells [j-1, j, j+1, j+2, j+3];
	// for convenience and uniformity we represent both using the
	// same combined structure of    [j-2, j-1, j, j+1, j+2, j+3].
	// std::valarray<double> u_plus(_actual_stencil_size);   // f_plus
	// std::valarray<double> u_minus(_actual_stencil_size);  // f_minus

	// For the purpose of linear stability (upwinding),
	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	// dfminus/du <= 0), is performed.
	// Lax-Friedrichs (LF) flux splitting (because it is the simplest
	// and smooth - we need the positive and negative fluxes to have
	// as many derivatives as the order of our finite-difference WENO)
	// is chosen here. `f_plus` uses a biased stencil with 1 point to
	// the left, while `f_minus` uses a biased stencil with 1 point to
	// the right.

	// So an LF flux	`numerical_flux`, f_hat(u_minus, u_plus),
	// a monotone numerical flux consistent with the physical one
	// (f_hat(u, u) = f(u)), will be construced at the end of
	// the loop below from `fhatminus` and `fhatplus` using
	// `f_minus` and `f_plus`.
	// N.B.! This can be replaced by an exact or approximate
	// Riemann solver (see Toro, 2009). Somehow...
	// Not every monotone flux can be writtenin the flux split form.
	// For example, the Godunov flux cannot.
	// std::valarray<double> numerical_flux(0., u.size());
	// f = std::valarray<double>(u.size());

	auto j_it_p = std::ranges::begin(u);  // f_plus

//	std::ranges::transform(
//			monotone_flux_components[0],
//			monotone_flux_components[1],
//			std::ranges::begin(f), [](const auto fp, const auto fm) {

//	})
	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size);
	auto u_plus = std::ranges::views::counted(j_it_p, 6);
	auto u_minus = u_plus | std::ranges::views::reverse;

	for (std::size_t j : shifted_index_range) {
		j_it_p = std::ranges::begin(u);  // f_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size);
		u_plus = std::ranges::views::counted(j_it_p, 6);

		u_plus_rec[j] = computeFHatWENO5FMReconstructionKernel(
			std::ranges::views::counted(
						std::ranges::begin(u_plus), 5), eps, p
		);


		u_minus = u_plus | std::ranges::views::reverse;
		u_minus_rec[j] = computeFHatWENO5FMReconstructionKernel(
			std::ranges::subrange(
				std::ranges::begin(u_minus),
						std::ranges::end(u_minus) - 1), eps, p
		);
//		u_minus_rec[j] = computeFHatWENO5JSReconstructionKernelRev(
//			std::ranges::views::counted(
//						std::ranges::begin(u_minus)+1, 5), eps, p
//		);
	}

	// std::cout << " done!" << "\n";
}


template <ArithmeticWith<numeric_val> T>
void updateGhostPointsTransmissive(
		std::ranges::common_range auto&& U,
		std::size_t left_bound_size = 3,
		std::optional<std::size_t> right_bound_size = std::nullopt) {
	/* Update ghost points in U with transmissive (Neumann) b.c.s. */

	if (!right_bound_size)
		right_bound_size.emplace(left_bound_size);

	const std::size_t n_full_size = std::ranges::size(U);
	const std::size_t right_start_index = (n_full_size
										   - right_bound_size.value());
	auto left_boundary_start = std::ranges::begin(U);
	auto right_boundary_start = std::ranges::begin(U);
	std::advance(right_boundary_start, right_start_index);

	// Transmissive b.c.s
	std::for_each_n(/*std::execution::par_unseq,*/
					left_boundary_start,
					left_bound_size,
					[&U, left_bound_size](auto& n) {
		n = U[left_bound_size];
	});
	// U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];

	std::for_each_n(/*std::execution::par_unseq,*/
					right_boundary_start,
					right_bound_size.value(),
					[&U, right_start_index](auto& n) {
		n = U[right_start_index-1];
	});
	// U[maxi+1] = U[maxi]; U[maxi+2] = U[maxi]; U[maxi+3] = U[maxi];
}


template <ArithmeticWith<numeric_val> T, typename... Args>
void timeOperator(
	std::ranges::common_range auto& U,
	std::ranges::common_range auto& flux,
	std::ranges::common_range auto& intermediate_fluxes,
	T t0, T dx, std::size_t n_size, T fin_t,
	auto&& timeStepFunction,
	auto&& calcMaxWaveSpeed,
	T cfl = 0.4,
	Args... opts_args
) {
	/* Time Operator of some problem: perform the time loop
	 * and solve it, storing the result in `U` and the numerical flux
	 * needed for the calculation in `flux`.
	 */

	T t = t0;
	T dt = 0.;

	T cpu = calcMaxWaveSpeed(U, dt);
	std::valarray<T> lam = std::valarray(cpu, 4);

	dt = cfl * dx / lam[0];

	while (t < fin_t) {
		if (t + dt > fin_t)
			dt = fin_t - t;

		timeStepFunction(U, flux, intermediate_fluxes,
				t, dt, dx, lam, n_size, opts_args...);

		t += dt;

		cpu = calcMaxWaveSpeed(U, dt);
		lam = std::valarray(cpu, 4);

		dt = cfl * dx / lam[0];
	}
	std::cout << "t = " << t << "\n";

	// return U;
}

#endif // WENO5_H
