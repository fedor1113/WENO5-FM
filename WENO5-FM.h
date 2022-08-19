#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
// #include <type_traits>
#include <utility>
#include <valarray>
#include <vector>
#include <algorithm>

#include "_vector4.h"

// template class Vector4<double>;

// using Vec4 = Vector4<double>;


// const std::vector<std::vector<double>> WmAm0 {};
std::array<std::array<const double, 6>, 6> WmAm0 {{
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 20./6, -31./6, 11./6}},
	{{0., 0., 0., -31./6, 50./6, -19./6}},
	{{0., 0., 0., 11./6, -19./6, 8./6}}
}};


std::array<std::array<const double, 6>, 6> WmAm1 {{
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 8./6, -13./6, 5./6, 0.}},
	{{0., 0., -13./6, 26./6, -13./6, 0.}},
	{{0., 0., 5./6, -13./6, 8./6, 0.}},
	{{0., 0., 0., 0., 0., 0.}}
}};

std::array<std::array<const double, 6>, 6> WmAm2 {{
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 8./6, -19./6, 11./6, 0., 0.}},
	{{0., -19./6, 50./6, -31./6, 0., 0.}},
	{{0., 11./6, -31./6, 20./6, 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}}
}};

std::array<std::array<const double, 6>, 6> minus_coefs[] = {
	WmAm0, WmAm1, WmAm2
};


std::array<std::array<const double, 6>, 6> WmAp0 {{
	{{8./6, -19./6, 11./6, 0., 0., 0.}},
	{{-19./6, 50./6, -31./6, 0., 0., 0.}},
	{{11./6, -31./6, 20./6, 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}}
}};

std::array<std::array<const double, 6>, 6> WmAp1 {{
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 8./6, -13./6, 5./6, 0., 0.}},
	{{0., -13./6, 26./6, -13./6, 0., 0.}},
	{{0., 5./6, -13./6, 8./6, 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}}
}};

std::array<std::array<const double, 6>, 6> WmAp2 {{
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 0., 0., 0., 0.}},
	{{0., 0., 20./6, -31./6, 11./6, 0.}},
	{{0., 0., -31./6, 50./6, -19./6, 0.}},
	{{0., 0., 11./6, -19./6, 8./6, 0}},
	{{0., 0., 0., 0., 0., 0.}}
}};

std::array<std::array<const double, 6>, 6> plus_coefs[] = {
	WmAp0, WmAp1, WmAp2
};

//	std::array<std::array<const T, 6>, 3> WmNminus {{
//		{{0., 0., 0., 11./6, -7./6, 2./6}},
//		{{0., 0., 2./6, 5./6, -1./6, 0.}},
//		{{0., -1./6, 5./6, 2./6, 0., 0.}}
//	}};

//	std::array<std::array<const T, 6>, 3> WmNplus {{
//		{{2./6, -7./6, 11./6, 0., 0., 0.}},
//		{{0., -1./6, 5./6, 2./6, 0., 0.}},
//		{{0., 0., 2./6, 5./6, -1./6, 0.}}
//	}};


template <typename T>
T gete(T rho, T p, T gamma) {
	if (rho != 0.)
		return p / (gamma-1.) / rho;

	return 0.;
}


//template <typename T>
//T getp(T rho, T e) {
//	return (GAMMA-1.) * rho * e;
//}


//template <typename T>
//T eFromConservative(T rho, T j, T rhoE) {
//	return (rhoE - 0.5 * j*j / rho) / rho;
//}


template <typename T>
Vector4<T> calcPhysicalFlux(T rho, T u, T p, T last, T gamma) {
	if (rho == 0) return Vector4<T>::ZERO;


	T e = gete(rho, p, gamma);
	return Vector4<T>(rho * u, p + rho*u*u,
				u*(p + rho*(e + 0.5*u*u)),
				u*last);
}


template <typename T>
Vector4<T> primitiveToConservative(Vector4<T> u, T gamma = 1.4) {
	// Conservative variables
	return Vector4<T>(u[0], u[0] * u[1],
		u[2] / (gamma - 1.) + 0.5 * u[0] * u[1] * u[1], u[3]);
}


template <typename T>
Vector4<T> conservativeToPrimitive(Vector4<T> q, T gamma) {
	// Primitive variables
	T rho = q[0];
	T u = q[1] / rho;
	T E = q[2] / rho;
	T p = (gamma - 1.) * rho * (E - 0.5*u*u);
	T e = gete(rho, p, gamma);

	return Vector4<T>(rho, u, p, e);
}


template <typename T>
Vector4<T> calcPhysicalFluxFromConservativeVec(Vector4<T> u, T gamma) {
//	return calcPhysicalFlux(u[0],
//			u[1] / u[0],
//			getp(u[0], eFromConservative(u[0], u[1], u[2])));
	Vector4<T> prim = conservativeToPrimitive(u, gamma);

	return calcPhysicalFlux(prim[0], prim[1], prim[2], prim[3], gamma);
}


template <typename T>
std::valarray<T> operator + (const std::valarray<T>& arr,
							 auto some_range) {
	std::valarray<T> res(arr.size());
	std::transform(
				std::begin(arr), std::end(arr),
				std::ranges::begin(some_range),
				std::begin(res), std::plus<>{});

	return res;
}


template <typename T>
std::valarray<T> operator - (const std::valarray<T>& arr,
							 auto some_range) {
	std::valarray<T> res(arr.size());
	std::transform(
				std::begin(arr), std::end(arr),
				std::ranges::begin(some_range),
				std::begin(res), std::minus<>{});

	return res;
}


template <typename T>
std::valarray<T>& operator += (std::valarray<T>& arr, auto some_range) {
	std::transform(
				std::begin(arr), std::end(arr),
				std::ranges::begin(some_range),
				std::begin(arr), std::plus<>{});
	return arr;
}


template <typename T>
std::valarray<T>& operator -= (std::valarray<T>& arr, auto some_range) {
	std::transform(
				std::begin(arr), std::end(arr),
				std::ranges::begin(some_range),
				std::begin(arr), std::minus<>{});
	return arr;
}


//template <typename T>
//Vector4<T> primitiveToConservative(Vector4<T> u, T gamma = 1.4) {
//	// Conservative variables
//	double G = 1005.9 / 717.09;
//	return Vector4<T>(u[0], u[1] / u[0], u[2] / (G-1.0), u[3]);
//}


//template <typename T>
//Vector4<T> conservativeToPrimitive(Vector4<T> q, T gamma = 1.4) {
//	// Primitive variables
//	T rho = q[0];
//	T u = q[1] / rho;
//	// T E = q[2] / rho;


//	T s1 = q[3] / q[0];
//	T s2 = 1. - s1;
//	T v = q[1] / q[0];
//	T cp1 = 1005.9; T cp2 = 627.83;  // Cp values for air and sf6
//	T cv1 = 717.09; T cv2 = 566.95;  // Cv values for air and sf6
//	T gammaeff = (cp1*s1+cp2*s2) / (cv1*s1+cv2*s2);  // Calculate an effective gamma
//	T p = (q[2]-q[0]*std::pow(v, 2.0)/2.0)*(gammaeff-1.);  // Calculate pressure from ch10.pdf, eq 10.2

//	return Vector4<T>(rho, u, p, q[3] / q[0]);
//}


//template <typename T>
//Vector4<T> calcPhysicalFluxFromConservativeVec(Vector4<T> u,
//											   T gamma = 1.4) {
//	T s1 = 1.;  // u[3] / u[0];
//	T s2 = 1. - s1;
//	T v = u[1] / u[0];

//	T cp1 = 1005.9; T cp2 = 627.83;  // Cp values for air and sf6
//	T cv1 = 717.09; T cv2 = 566.95;  // Cv values for air and sf6
//	T gammaeff = (cp1*s1+cp2*s2) / (cv1*s1+cv2*s2);  // Calculate an effective gamma
//	T P = (u[2]-u[0]*std::pow(v, 2.0)/2.0)*(gammaeff-1.);  // Calculate pressure from ch10.pdf, eq 10.2
//	//	ret[0] = u[1]
//	//	ret[1] = rho*np.power(v,2.0)+P;
//	//	ret[2] = (u[2]+P)*v;
//	//	ret[3] = v*u[3];
//	//	return ret;

//	return Vector4<T>(u[1], u[0]*v*v + P, (u[2]+P)*v, u[3]*v);
//}


template <typename T, typename T1, typename T2>
std::valarray<T> vecMatDot(const T1& vec, const T2& mat) {
	/* Multiply a vector by a matrix (a sum product over
	 * the last axis of mat and vec).
	 */

	std::valarray<T> result(std::ranges::size(mat));

	for (std::size_t k = 0; k < std::ranges::size(mat)
			&& k < std::ranges::size(vec); ++ k)
		result[k] = std::inner_product(std::begin(mat[k]),
									   std::end(mat[k]),
									   std::begin(vec), 0.);

	return result;
}


template <typename T>
std::valarray<Vector4<T>> calcPhysFlux(
		const std::valarray<Vector4<T>>& u_arr, T gamma = 1.4) {
	/* Calculate physical fluxes of conservative variables
	 * in all points in the computation domain u_arr.
	 */

//	return u_arr.apply(
//				[](Vector4<T> v) {
//		return calcPhysicalFluxFromConservativeVec<T>(v);
//	});  // f_arr
	std::valarray<Vector4<T>> res(std::ranges::size(u_arr));

	std::transform(std::begin(u_arr), std::end(u_arr), std::begin(res),
				   [&gamma](Vector4<T> v) {
		return calcPhysicalFluxFromConservativeVec<T>(v, gamma);
	});

	return res;  // f_arr
}


//template <typename T>
//T calcMaxWaveSpeedD(const std::valarray<Vector4<T>>& u_arr) {
//	/* Calculate df/du. */

//	std::size_t k = 0;

//	T cp1 = 1005.9; T cp2 = 627.83; // Cp values for air and sf6
//	T cv1 = 717.09; T cv2 = 566.95; // Cv values for air and sf6

//	std::valarray<T> s1(u_arr.size());
//	for (k = 0; k < u_arr.size(); ++ k)
//		s1[k] = 1.; // s2 = W[4,:]/W[0,:]
//	std::valarray<T> s2 = 1. - s1;

//	std::valarray<T> arr_G = (s1*cp1+s2*cp2) / (s1*cv1+s2*cv2);
//	// std::valarray<T> arr_G(1.4, U.size());
//	std::valarray<T> a0(u_arr.size());
//	for (k = 0; k < u_arr.size(); ++ k)
//		a0[k] = arr_G[k] * u_arr[k][2] * (arr_G[k]-1.) / u_arr[k][0];
//	a0 = std::sqrt(a0);

//	return a0.max();
//}


template <typename T>
T calcSquareSoundSpeed(T rho, T rho_v, T rho_E, T gamma = 1.4) {
	return gamma * (gamma - 1.)
			* (rho_E - rho_v * rho_v * 0.5 / rho) / rho;
}


template <typename T>
auto calcMaxWaveSpeedD(const auto& u_arr, T gamma = 1.4) {
	/* Calculate |df/du| for 1D Euler eq'ns. */

	const std::size_t size = std::ranges::size(u_arr);
	std::valarray<T> a0(size);
	std::ranges::transform(std::as_const(u_arr),
						   std::begin(a0),
						   [gamma](const auto& u_arr_vec_pt) -> T {
		return (std::sqrt(std::abs(calcSquareSoundSpeed(
							u_arr_vec_pt[0],
							u_arr_vec_pt[1],
							u_arr_vec_pt[2], gamma)))
				+ std::abs(u_arr_vec_pt[1] / u_arr_vec_pt[0]));
	});
//	std::size_t k = 0;
//	for (k = 0; k < size; ++ k)
//		a0[k] = gamma * u_arr[k][2] * (gamma-1.) / u_arr[k][0];
	// a0 = std::sqrt(std::abs(a0));

	// std::ranges::for_each(a0, [](auto& n) { n = std::abs(n); });
	// std::ranges::for_each(a0, [](auto& n) { n = std::sqrt(n); });

//	std::valarray<T> a1(u_arr.size());
//	std::valarray<T> a2(u_arr.size());
//	for (k = 0; k < u_arr.size(); ++ k) {
//		a1[k] = a0[k] + u_arr[k][1]/u_arr[k][0];
//		a2[k] = a0[k] - u_arr[k][1]/u_arr[k][0];
//	}

//	return std::max({a0.max(), a1.max(), a2.max()});
	// return a1.max();
	// return a0.max();
	return *std::ranges::max_element(a0);
}


template <typename T>
//std::array<T, 3> smoothness_indicators(const T1& f_stencil) {
std::valarray<T> betaSmoothnessIndicators(std::span<T, 5> f_stencil) {
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
	std::valarray<T> res(3);
	T f_prev2 = f_stencil[0];
	T f_prev1 = f_stencil[1];
	T f_curr0 = f_stencil[2];
	T f_next1 = f_stencil[3];
	T f_next2 = f_stencil[4];

	T beta_0 = ((13./12.) * std::pow(f_prev2 - 2.*f_prev1 + f_curr0, 2)
				+ (1./4.) * std::pow(f_prev2 - 4.*f_prev1 + 3.*f_curr0, 2));

	T beta_1 = ((13./12.) * std::pow(f_prev1 - 2.*f_curr0 + f_next1, 2)
				+ (1./4.) * std::pow((f_prev1 - f_next1), 2));

	T beta_2 = ((13./12.) * std::pow(f_curr0 - 2.*f_next1 + f_next2, 2)
				+ (1./4.) * std::pow(3.*f_curr0 - 4.*f_next1 + f_next2, 2));

	res[0] = beta_0;
	res[1] = beta_1;
	res[2] = beta_2;

	return res;
}


template <typename T, typename T1>
std::valarray<T> betaSmoothnessIndicatorsMat(
	const T1& f_stencil,
	std::array<std::array<const T, 6>, 6> _coefs[] = plus_coefs
) {
	/* Return the smoothness indicators beta_k, k=0,1,2
	 * for each of the 3 substencils of `f_stencil`.
	 */

	std::valarray<T> res(3);

	for (std::size_t k = 0; k < 3; ++ k) {
		// for (std::size_t k = 0; k < half_size + 1; ++ k)
		res[k] = std::inner_product(
			std::begin(f_stencil),
			std::end(f_stencil),
			std::begin(vecMatDot(f_stencil, _coefs[k])),
			0.
		);
	}


	return res;
}


template <typename T>
std::valarray<T> f3OrdReconstructionFromStencil(
		std::span<T, 5> f_stencil) {
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


template <typename T>
T henrickGMappingForLambda(T lambda_weno_weight,
						   T lambda_ideal = 1./3.) {
	/* The mapping function g by Henrick modified for symmetric
	 * lambda-weights by Zheng Hong, Zhengyin Ye and Kun Ye.
	 */

	T square_ideal = lambda_ideal * lambda_ideal;
	return lambda_weno_weight * (lambda_ideal
								 + square_ideal
								 - 3. * lambda_ideal * lambda_weno_weight
								 + lambda_weno_weight * lambda_weno_weight)
							  / (square_ideal
								 + lambda_weno_weight
								   * (1. - 2.*square_ideal));
}


template <typename T>
std::valarray<T> alphaWENO5FMWeights(
	const std::valarray<T>&& beta_IS_coefs,
	T epsilon = 1e-40,
	T p = 2.
) {
	/* Compute appropriate alpha(α)-weights for the WENO5-FM scheme,
	 * by which the inverses of smoothness indicators are meant,
	 * so inverse beta(β) with the caveat of aritificially finite
	 * answers using the added epsilon-parameter to the beta-weights
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * `p` controls (increases) the amount of numerical dissipation
	 * (it's recommended to take it = r-1 for 2r-1 order schemes,
	 * so 2 for WENO5).
	 */

	return 1. / std::pow(epsilon + beta_IS_coefs, p);
}


template <typename T>
std::valarray<T> lambdaWENO5FMWeights(
	const std::valarray<T>&& alpha_weights
) {
	/* FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
	 * for WENO5-FM or WENO5-ZM
	 * due to Zheng Hong, Zhengyin Ye and Kun Ye:
	 * lambda_weights = alpha_weights / alpha_weights.sum();
	 */

	return alpha_weights / alpha_weights.sum();
}


template <typename T>
std::valarray<T> omegaWENO5FMWeights(
	const std::valarray<T>&& lambda_weights
) {
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

	std::valarray<T> alpha_weights = lambda_weights.apply(gMap);
	// α*-weights

	// From α*=g(λ_k) and d_k we get the new corrected resultant
	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
	// omega_weights = d_lin_weights * alpha_weights;
	std::valarray<T> omega_weights = d_lin_weights * alpha_weights;
	omega_weights /= omega_weights.sum();

	return omega_weights;
}


template <typename T>
T computeFHatWENO5JSReconstructionKernel(std::span<T, 5> f_stencil,
										 T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical fluxes
	 * `fhatplus` / `fhatminus` at a point j+0 for a given stencil
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
	// Computationally, we will only need to remember
	// some 3 weights at each WENO5 reconstruction step,
	// so this array will be more than enough, but,
	// for ease of understanding, we will create
	// a bunch of other arrays.
//	std::valarray<T> alpha_weights(3);
//	std::valarray<T> lambda_weights(3);
//	std::valarray<T> omega_weights(3);
//	std::valarray<T> eno_reconstructed_f(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {0.1, 0.6, 0.3};
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
//		std::begin(omega_weights), std::end(omega_weights),
//		std::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	 std::valarray<
		 T
	 > eno_reconstructed_f = f3OrdReconstructionFromStencil(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];


	return f_hat;
}


template <typename T>
T computeFHatWENO5FMReconstructionKernel(std::span<T, 5> f_stencil,
										 T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical fluxes
	 * `fhatplus` / `fhatminus` at a point j+0 for a given stencil
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
	// Computationally, we will only need to remember
	// some 3 weights at each WENO5 reconstruction step,
	// so this array will be more than enough, but,
	// for ease of understanding, we will create
	// a bunch of other arrays.
//	std::valarray<T> alpha_weights(3);
//	std::valarray<T> lambda_weights(3);
//	std::valarray<T> omega_weights(3);
//	std::valarray<T> eno_reconstructed_f(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators(f_stencil);

	// non-linear non-scaled (α-)weights
//	std::valarray<T> d_lin_weights = {0.1, 0.6, 0.3};
//	std::valarray<T> alpha_weights = d_lin_weights / std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// Instead we use:
	std::valarray<T> alpha_weights = alphaWENO5FMWeights(
		std::move(beta_IS_coefs), eps, p
	);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
//	 std::valarray<T> omega_weights = alpha_weights / alpha_weights.sum();

	// FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
	// due to Zheng Hong, Zhengyin Ye and Kun Ye:
//	 lambda_weights = alpha_weights / alpha_weights.sum();
	std::valarray<T> lambda_weights = lambdaWENO5FMWeights(
		std::move(alpha_weights)
	);

	std::valarray<T> omega_weights = omegaWENO5FMWeights(
		std::move(lambda_weights)
	);

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of f_{i+1/2}
	// via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::begin(omega_weights), std::end(omega_weights),
//		std::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	 std::valarray<
		 T
	 > eno_reconstructed_f = f3OrdReconstructionFromStencil(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];


	return f_hat;
}


//template <typename T>
//T computeFHatWENO5FMReconstructionKernelRev(std::span<T, 5> f_stencil,
//											T eps = 1e-40, T p = 2.) {
//	/* Calculate (reconstruct) one of the two split monotone numerical fluxes
//	 * `fhatplus` / `fhatminus` at a point j+0 for a given stencil
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

//	beta_IS_coefs[0] = ((13./12.) * std::pow(f_prev2 - 2.*f_prev1 + f_curr0, 2)
//				+ (1./4.) * std::pow(f_prev2 - 4.*f_prev1 + 3.*f_curr0, 2));

//	beta_IS_coefs[1] = ((13./12.) * std::pow(f_prev1 - 2.*f_curr0 + f_next1, 2)
//				+ (1./4.) * std::pow((f_prev1 - f_next1), 2));

//	beta_IS_coefs[2] = ((13./12.) * std::pow(f_curr0 - 2.*f_next1 + f_next2, 2)
//				+ (1./4.) * std::pow(3.*f_curr0 - 4.*f_next1 + f_next2, 2));

//	// non-linear non-scaled (α-)weights
//	std::valarray<T> d_lin_weights = {0.3, 0.6, 0.1};
//	std::valarray<T> alpha_weights = d_lin_weights / std::pow(eps + beta_IS_coefs, p);
//	// — we have no need of them in WENO-FM!

//	// Instead we use:
////	std::valarray<T> alpha_weights = alphaWENO5FMWeights(
////		std::move(beta_IS_coefs), eps, p
////	);

//	// scaled (normalized) non-linear (ω-)weights (ENO weights)
//	 std::valarray<T> omega_weights = alpha_weights / alpha_weights.sum();

//	// FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
//	// due to Zheng Hong, Zhengyin Ye and Kun Ye:
////	 lambda_weights = alpha_weights / alpha_weights.sum();
////	std::valarray<T> lambda_weights = lambdaWENO5FMWeights(
////		std::move(alpha_weights)
////	);

////	std::valarray<T> omega_weights = omegaWENO5FMWeights(
////		std::move(lambda_weights)
////	);

//	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of f_{i+1/2}
//	// via linear combinations with WmNplus coefficients
//	// for each substencil which is then used to calculate
//	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
//	// using the nonlinear weights [ω]
////	f_hat = std::inner_product(
////		std::begin(omega_weights), std::end(omega_weights),
////		std::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
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


template <typename T>
std::array<std::valarray<T>, 2> splitFluxAsLaxFriedrichs(
		const auto& u, const auto& f, T alpha) {
	/* Global Lax-Friedrichs (LF) flux splitting.
	 *
	 * For the purpose of linear stability (upwinding),
	 * a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	 * dfminus/du <= 0), needs to be performed in FD WENO5.
	 * LF flux splitting remains especially simple, while still
	 * retaining the necessary number of derivatives,
	 * so it's perfect for the job.
	 */

	std::array<std::valarray<T>, 2> monotone_lf_flux_components {
		std::valarray<T>(std::ranges::size(f)),
		std::valarray<T>(std::ranges::size(f))
	};

	// 	std::valarray<T> f_plus = 0.5 * (f + alpha * u);
	//	std::valarray<T> f_plus = f / alpha;
	//	f_plus += u;
	//	f_plus *= 0.5 * alpha;
	std::transform(std::begin(f), std::end(f),
		std::begin(u), std::begin(monotone_lf_flux_components[0]),
		[&alpha](T f_pt, T u_pt) {
		return 0.5 * (f_pt + alpha * u_pt);
	});

	//	std::valarray<T> f_minus = 0.5 * (f - alpha * u);
	//	std::valarray<T> f_minus = f / alpha;
	//	f_minus -= u;
	//	f_minus *= 0.5 * alpha;
	std::valarray<T> f_minus(std::ranges::size(f));
	std::transform(std::begin(f), std::end(f),
		std::begin(u), std::begin(monotone_lf_flux_components[1]),
		[&alpha](T f_pt, T u_pt) {
		return 0.5 * (f_pt - alpha * u_pt);
	});

	return monotone_lf_flux_components;
}


// FD WENO5FM (WENO5-FM) - method
// Reconstruction based on LF flux splitting + improved mapped WENO
// of 5th order
// (see Mapped weighted essentially non-oscillatory schemes:
// achieving optimal order near critical points, 2005 by Henrick et al.)
// and 'An improved WENO-Z scheme with symmetry-preserving mapping'
// by Zheng Hong, Zhengyin Ye and Kun Ye, 2020
template <typename T>
void calcHydroStageWENO5FM(const auto& u,
						   T t,
						   T lam,
						   auto& f,
						   std::size_t nSize,
						   T eps = 1e-40,
						   T p = 2.) {
	/* Component-wise finite-difference WENO5FM (WENO5-FM) - space
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
	const std::size_t nGhostCells = (stencil_size + 1) / 2;
	const std::size_t mini = nGhostCells;
	const std::size_t maxi = nGhostCells + nSize - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini + 2, maxi - 2};
	auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 2};
	// [g      g      g      i      i      i      i      i      i      ...]
	// {0      1      2      3      4      5}     6      7      8      ...
	//                |             |
	//                j        nGhostCells

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
	std::valarray<T> u_plus(_actual_stencil_size);   // f_plus
	std::valarray<T> u_minus(_actual_stencil_size);  // f_minus

	// For the purpose of linear stability (upwinding),
	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	// dfminus/du <= 0), is performed.
	// Lax-Friedrichs (LF) flux splitting (because it is the simplest
	// and smooth - we need the positive and negative fluxes to have
	// as many derivatives as the order of our finite-difference WENO)
	// is chosen here. `f_plus` uses a biased stencil with 1 point to
	// the left, while `f_minus` uses a biased stencil with 1 point to
	// the right.
	// T alpha = abs(U).max();
	// T alpha = std::sqrt(u.max().square());
	T alpha = std::abs(lam);  // α = max |df/du|

	T fhatminus = 0.;
	T fhatplus = 0.;

	std::array<
			std::valarray<T>, 2
			> monotone_flux_components = splitFluxAsLaxFriedrichs(
				u, f, alpha
				);

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

	auto j_it_p = std::begin(monotone_flux_components[0]);  // f_plus
	auto j_it_m = std::begin(monotone_flux_components[1]);  // f_minus

	for (auto j : shifted_index_range) {
	// for (std::size_t j = 10; j < 11; ++ j) {
		j_it_p = std::begin(monotone_flux_components[0]); // f_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size);
		std::copy_n(j_it_p, u_plus.size(), std::begin(u_plus));

		j_it_m = std::begin(monotone_flux_components[1]);  // f_minus
		std::advance(j_it_m, j + half_size + 1 - stencil_size);
		std::copy_n(j_it_m, u_minus.size(), std::begin(u_minus));

		fhatplus = computeFHatWENO5FMReconstructionKernel(
			std::span<T, 5>{std::begin(u_plus), 5}, eps, p
		);

		std::reverse(std::begin(u_minus), std::end(u_minus));
		fhatminus = computeFHatWENO5FMReconstructionKernel(
			std::span<T, 5>{std::begin(u_minus), 5}, eps, p
		);
//		fhatminus = computeFHatWENO5FMReconstructionKernelRev(
//			std::span<T, 5>{std::begin(u_minus)+1, 5}, eps, p
//		);

		f[j] = fhatplus + fhatminus;

//		for (auto n : f)
//			std::cout << n << " ";
//		std::cout << "\n";
	}

//	for (std::size_t j = maxi+1; j < u.size(); ++ j) {
//		// numerical_flux[j] = 0;
//		f[j] = 0;
//	}

//	for (std::size_t j = 0; j < 2; ++ j) {
//		// numerical_flux[j] = 0;
//		f[j] = 0;
//	}

	// return numerical_flux;
	// f = numerical_flux;
	// std::cout << " done!" << "\n";
}


template <typename T>
void updateGhostPoints(
		auto& U,
		std::size_t left_bound_size = 3,
		std::optional<std::size_t> right_bound_size = std::nullopt) {
	/* Update ghost points in U with transmissive (Neumann) b.c.s. */

	//	for (k = nGhostCells; k < nGhostCells + nSize; ++ k) {
	//			U[k] = ms[k-nGhostCells].W;
	//	}

	// const std::size_t n_full_size = U.size();
	// const std::size_t n_ghost_cells = 3;  // 3
	// const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_full_size - n_ghost_cells - 1;
	// ... n_full_size-4   n_full_size-3   n_full_size-2   n_full_size-1 ]
	// ...     maxi
	// const std::size_t n_size = n_full_size - 2 * n_ghost_cells;
	if (!right_bound_size)
		right_bound_size.emplace(left_bound_size);

	const std::size_t n_full_size = std::ranges::size(U);
	const std::size_t right_start_index = (n_full_size
										   - right_bound_size.value());
	auto left_boundary_start = std::begin(U);
	auto right_boundary_start = std::begin(U);
	std::advance(right_boundary_start, right_start_index);

	// Transmissive b.c.s
	std::ranges::for_each_n(left_boundary_start,
							left_bound_size,
							[&U, left_bound_size](auto& n) {
		n = U[left_bound_size];
	});
	// U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];

	std::ranges::for_each_n(right_boundary_start,
							right_bound_size.value(),
							[&U, right_start_index](auto& n) {
		n = U[right_start_index];
	});
	// U[maxi+1] = U[maxi]; U[maxi+2] = U[maxi]; U[maxi+3] = U[maxi];
}


template <typename T>
std::valarray<Vector4<T>> calcFluxComponentWise(
		const std::valarray<Vector4<T>>& U,
		T t, const std::valarray<T>& lam,
		std::size_t nSize,
		T eps = 1e-40,
		T p = 2.) {
	std::valarray<Vector4<T>> res = calcPhysFlux(U);
	std::valarray<std::valarray<T>> component_fs(
				std::valarray<T>(U.size()), 4);

	std::size_t k = 0;

	for (std::size_t j = 0; j < 4; ++ j)
		for (k = 0; k < U.size(); ++ k) {
			component_fs[j][k] = res[k][j];
		}

	auto getVector4Component = [](
			const Vector4<T>& x,
			std::size_t k) { return x[k]; };

	for (size_t k : std::ranges::iota_view{0, 4}) {
		auto kThVector4Component = [&k, &getVector4Component](
				const Vector4<T>& x) {
			return getVector4Component(x, k);
		};

		calcHydroStageWENO5FM<T>(
			std::ranges::views::transform(U, kThVector4Component),
			t, lam[k],
			component_fs[k], nSize, eps, p
		);
	}

	for (std::size_t j = 0; j < 4; ++ j)
		for (k = 0; k < U.size(); ++ k)
			res[k][j] = component_fs[j][k];

	return res;
}


//template <typename T>
//std::valarray<Vector4<T>> calcFluxComponentWise(
//		const std::valarray<Vector4<T>>& U,
//		T t, const std::valarray<T>& lam,
//		std::size_t nSize,
//		T eps = 1e-40,
//		T p = 2.) {
//	std::valarray<Vector4<T>> res = calcPhysFlux(U);
//	std::valarray<std::valarray<T>> components(
//				std::valarray<T>(U.size()), 4);
//	std::valarray<std::valarray<T>> component_fs(
//				std::valarray<T>(U.size()), 4);

//	std::size_t k = 0;

//	for (std::size_t j = 0; j < 4; ++ j)
//		for (k = 0; k < U.size(); ++ k) {
//			components[j][k] = U[k][j];
//			component_fs[j][k] = res[k][j];
//		}

//	k = 0;
//	// for (auto component : components) {
//	for (k = 0; k < 4; ++ k) {
//		calcHydroStageWENO5FM<T>(components[k], t, lam[k],
//								 component_fs[k], nSize, eps, p);
//		// ++ k;
//	}

//	for (std::size_t j = 0; j < 4; ++ j)
//		for (k = 0; k < U.size(); ++ k)
//			res[k][j] = component_fs[j][k];

//	return res;
//}


template <typename T>
std::valarray<Vector4<T>> calcdSpace(const std::valarray<Vector4<T>>& U,
							   T t,
							   T dx,
							   const std::valarray<T>& lam,
							   std::size_t nSize,
							   T eps = 1e-40,
							   T p = 2.) {
	std::valarray<Vector4<T>> dflux(Vector4<T>::ZERO, U.size());
	std::valarray<Vector4<T>> lf = calcFluxComponentWise<T>(
				U, t, lam, nSize, eps, p
				);

	const std::size_t ghost_point_number = 3;

	std::slice Nweno(ghost_point_number, nSize, 1);
	std::slice Nweno_shifted_by_neg_1(ghost_point_number-1, nSize, 1);
//	std::slice Nweno(1, U.size()-1, 1);
//	std::slice Nweno_shifted_by_neg_1(0, U.size()-1, 1);

	std::valarray<Vector4<T>> f_mn = lf[Nweno_shifted_by_neg_1];
	std::valarray<Vector4<T>> f_pl = lf[Nweno];


	dflux[Nweno] = -(f_pl - f_mn) / Vector4<T>(dx);
	// dflux[3] -= lf[0] / Vector4<T>(dx);
	// dflux[std::ranges::size(U)-1-ghost_point_number]
	// 	+= lf[std::ranges::size(U)-1-ghost_point_number] / Vector4<T>(dx);

	// dflux += source terms...

	return dflux;
}


template <typename T>
void advanceTimestepTVDRK3(
		auto& U,
		auto& dflux,
		auto& Y2,
		auto& Y3,
		T t, T dt, T dx,
		const auto& lam,
		std::size_t n_size
//		auto& calcdSpace,
//		auto& updateGhostPoints,
//		auto& calcFlux,
//		auto& calcSource,
		/*auto& opts*/) {
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
	dflux = calcdSpace<T>(U, t, dx, lam, n_size);  // L1 = L[u^n]
	// std::valarray<Vector4<T>> res(Vector4<T>::ZERO, U.size());
	// L[u] = (-) dF[u]/dx

	// Y2 = U + Vector4<T>(dt) * dflux;
	std::ranges::transform(U,
						   dflux,
						   std::begin(Y2),
						   [dt](const auto u, const auto df) {
		return u + dt * df;
	});
	// u(1) = u^n + Δt L[u^n]

	updateGhostPoints<T>(Y2);


	// ------------------------Second Stage---------------------------
	dflux = calcdSpace<T>(Y2, t, dx, lam, n_size);  // L2 = L[u(1)]

	// Y3 = Vector4<T>(3.)*U + Vector4<T>(dt) * dflux + Y2;
	// Y3 *= Vector4<T>(0.25);
	for (std::size_t k = 0; k < std::ranges::size(U); ++ k)
		Y3[k] = (3. * U[k] + dt * dflux[k] + Y2[k]) * 0.25;
	// u(2) = 0.75 * u^n + 0.25 * u(1) + 0.25 * Δt L[u(1)]

	updateGhostPoints<T>(Y3);


	// ------------------------Third Stage----------------------------
	dflux = calcdSpace<T>(Y3, t, dx, lam, n_size);  // L3 = L[u(2)]


	// U += Vector4<T>(2.) * (Y3 + Vector4<T>(dt) * dflux);
	// U *= Vector4<T>(1./3.);
	std::ranges::transform(Y3,
						   dflux,
						   std::begin(Y3),
						   [dt](const auto u, const auto df) {
		return 2. * (u + dt * df);
	});
	std::ranges::transform(U,
						   Y3,
						   std::begin(U),
						   [dt](const auto u, const auto y) {
		return (u + y) * (1./3.);
	});
	// u^(n+1) = (1/3) * u^n + (2/3) * u(2) + (2/3) * Δt L[u(2)]

	updateGhostPoints<T>(U);

	// return U;
	// U = res;
}


template <typename T>
std::valarray<Vector4<T>> integrateRiemannProblem(
	auto& U,
	auto& flux,
	auto& Y2,
	auto& Y3,
	T t0, T dx, std::size_t nSize,
	T fin_t, T cfl = 0.4
) {
	/* ... */

	T t = t0;

	T cpu = calcMaxWaveSpeedD<T>(U);
	std::valarray<T> lam = std::valarray(cpu, 4);

	T dt = cfl * dx / lam[0];

	while (t < fin_t) {
		if (t + dt > fin_t)
			dt = fin_t - t;

		advanceTimestepTVDRK3<T>(U, flux, Y2, Y3, t, dt, dx, lam, nSize);
		// std::cout << n*dt << "\n";
		// std::cout << t << "\n";

		t += dt;

		 cpu = calcMaxWaveSpeedD<T>(U);
		 lam = std::valarray(cpu, 4);

		dt = cfl * dx / lam[0];
	}
	std::cout << "t = " << t << "\n";

	return U;
}


template <typename T>
std::function<Vector4<T>(Vector4<T>, T)> primitiveToConservativeU
	= primitiveToConservative<T>;


template <typename T>
void prepareRiemannProblem(
	auto& u_init, auto& x, T gamma,
	T rho_left, T v_left, T p_left, T e_left, T rhoE_left,
	T rho_right, T v_right, T p_right, T e_rightR, T rhoE_right,
	T q0, T l_min, T l_max,
	std::function<Vector4<T>(Vector4<T>, T)> primitiveToConservativeU,
	std::size_t mesh_size
) {
	/* ... */

//	const T Runiv = 8.314;  // Universal gas constant, [J/K-mol]
//	T mu1 = 0.0289647; T mu2 = 0.14606;  // Mollecular weights of air and sf6 [kg/mol]
//	T R1 = Runiv/mu1; T R2 = Runiv/mu2;  // Specific gas constants
//	T cp1 = 1005.9; T cp2 = 627.83;  // Cp values for air and sf6
//	T cv1 = 717.09; T cv2 = 566.95;  // Cv values for air and sf6

	std::size_t computational_domain_size = mesh_size;
	const std::size_t n_ghost_points = 3;
	std::size_t full_mesh_size = computational_domain_size
		+ 2*n_ghost_points;
	T dx = (l_max - l_min) / mesh_size;  // [L]

	x = std::valarray<T>(0., full_mesh_size);
	u_init = std::valarray<Vector4<T>>(Vector4<T>::ZERO, full_mesh_size);

	std::size_t k = 0;
	for (k = 0; k < n_ghost_points; ++ k)
		x[k] = l_min - dx * (n_ghost_points - k);

	x[n_ghost_points] = l_min;
	for (k = n_ghost_points + 1; k < full_mesh_size; ++ k)
		x[k] = x[k-1] + dx;

	std::size_t x0_coord = 0;
	while (x[x0_coord] < q0)
		++ x0_coord;

//	T Tatm = 293.0;  // [K], approx 70 F
//	T Patm = 101300.0;  // [Pa], Atmospheric Pressure
//	T Rhoatm = Patm / (Tatm * R1);  // Density of the first gas at STP

// 	T G = gamma;
	Vector4<T> vec(rho_left, v_left, p_left, 0.);
	vec = primitiveToConservativeU(vec, gamma);
	for (k = 0; k <= x0_coord; ++ k) {
		u_init[k] = vec;
		// u_init[k][3] = 1. * u_init[k][0];
	}

	vec = primitiveToConservativeU(Vector4(rho_right, v_right, p_right, 0.), gamma);
	for (k = x0_coord + 1; k < full_mesh_size; ++ k) {
		u_init[k] = vec;
		// u_init[k][3] = 1. * u_init[k][0];
	}
}


template <typename T>
std::valarray<Vector4<T>> solve1DRiemannProblemForEulerEq(
	auto& u_init, auto& x, T gamma,
	T rho_left, T v_left, T p_left, T e_left, T rhoE_left,
	T rho_right, T v_right, T p_right, T e_rightR, T rhoE_right,
	T q0, // Initial coordinate of the discontinuity
	T t0, T t_max, T l_min, T l_max,
	std::function<Vector4<T>(Vector4<T>, T)> primitiveToConservativeU,
	std::size_t mesh_size, T cfl = 0.4
) {
	/* ... */

	prepareRiemannProblem<T>(
		u_init, x, gamma,
		rho_left, v_left, p_left, e_left, rhoE_left,
		rho_right, v_right, p_right, e_rightR, rhoE_right,
		q0, l_min, l_max, primitiveToConservativeU, mesh_size
	);

	std::valarray<Vector4<T>> flux(Vector4<T>::ZERO,
								  std::ranges::size(u_init));
	std::valarray<Vector4<T>> y2(Vector4<T>::ZERO,
								 std::ranges::size(u_init));
	std::valarray<Vector4<T>> y3(Vector4<T>::ZERO,
								 std::ranges::size(u_init));

	u_init = integrateRiemannProblem<T>(u_init, flux, y2, y3,
		t0, (l_max-l_min) / mesh_size, mesh_size, t_max, cfl);

	return u_init;
}


template <typename T>
std::valarray<Vector4<T>> solve1DRiemannProblemForEulerEq(
	auto& u_init, auto& x, T gamma,
	T rho_left, T v_left, T p_left,
	T rho_right, T v_right, T p_right,
	T q0, // Initial coordinate of the discontinuity
	T t0, T t_max, T l_min, T l_max,
	std::function<Vector4<T>(Vector4<T>, T)> primitiveToConservativeU,
	std::size_t mesh_size, T cfl = 0.4
) {
	double e_left = p_left / (gamma - 1.) / rho_left;
	double e_right = p_right / (gamma - 1.) / rho_right;
	double E_left = (e_left + (v_left*v_left)/2.);
	double E_right = (e_right + (v_right*v_right)/2.);

	return solve1DRiemannProblemForEulerEq(
		u_init, x, gamma,
		rho_left, v_left, p_left, e_left, E_left,
		rho_right, v_right, p_right, e_right, E_right,
		q0, t0, t_max, l_min, l_max,
		primitiveToConservativeU, mesh_size, cfl
	);
}
