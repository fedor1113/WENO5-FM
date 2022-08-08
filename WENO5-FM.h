#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
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

// const double GAMMA = 1.4;


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


//template <typename T>
//T gete(T rho, T p) {
//	if (rho != 0.)
//		return p / (GAMMA-1.) / rho;

//	return 0.;
//}


//template <typename T>
//T getp(T rho, T e) {
//	return (GAMMA-1.) * rho * e;
//}


//template <typename T>
//T eFromConservative(T rho, T j, T rhoE) {
//	return (rhoE - 0.5 * j*j / rho) / rho;
//}


//template <typename T>
//Vector4<T> calcPhysicalFlux(T rho, T u, T p, T last) {
//	if (rho == 0) return Vector4<T>::ZERO;


//	T e = gete(rho, p);
//	return Vector4<T>(rho * u, p + rho*u*u,
//				u*(p + rho*(e + 0.5*u*u)),
//				u*last);
//}


//template <typename T>
//Vector4<T> conservativeToPrimitive(Vector4<T> q) {
//	// Primitive variables
//	T rho = q[0];
//	T u = q[1] / rho;
//	T E = q[2] / rho;
//	T p = (GAMMA - 1.) * rho * (E - 0.5*u*u);

//	return Vector4<T>(rho, u, p, q[3]);
//}


//template <typename T>
//Vector4<T> calcPhysicalFluxFromConservativeVec(Vector4<T> u) {
////	return calcPhysicalFlux(u[0],
////			u[1] / u[0],
////			getp(u[0], eFromConservative(u[0], u[1], u[2])));
//	Vector4<T> prim = conservativeToPrimitive(u);

//	return calcPhysicalFlux(prim[0], prim[1], prim[2], prim[3]);
//}


template <typename T>
std::valarray<T> operator + (const std::valarray<T>& arr,
							 auto some_range) {
	std::valarray<T> res(arr.size());
	std::transform(std::begin(arr), std::end(arr),
				   std::ranges::begin(some_range),
				   std::begin(res), std::plus<>{});

	return res;
}


template <typename T>
std::valarray<T> operator - (const std::valarray<T>& arr,
							 auto some_range) {
	std::valarray<T> res(arr.size());
	std::transform(std::begin(arr), std::end(arr),
				   std::ranges::begin(some_range),
				   std::begin(res), std::minus<>{});

	return res;
}


template <typename T>
std::valarray<T>& operator += (std::valarray<T>& arr, auto some_range) {
	std::transform(std::begin(arr), std::end(arr),
				   std::ranges::begin(some_range),
				   std::begin(arr), std::plus<>{});
	return arr;
}


template <typename T>
std::valarray<T>& operator -= (std::valarray<T>& arr, auto some_range) {
	std::transform(std::begin(arr), std::end(arr),
				   std::ranges::begin(some_range),
				   std::begin(arr), std::minus<>{});
	return arr;
}


template <typename T>
Vector4<T> conservativeToPrimitive(Vector4<T> q) {
	// Primitive variables
	T rho = q[0];
	T u = q[1] / rho;
	// T E = q[2] / rho;


	T s1 = q[3] / q[0];
	T s2 = 1. - s1;
	T v = q[1] / q[0];
	T cp1 = 1005.9; T cp2 = 627.83;  // Cp values for air and sf6
	T cv1 = 717.09; T cv2 = 566.95;  // Cv values for air and sf6
	T gammaeff = (cp1*s1+cp2*s2) / (cv1*s1+cv2*s2);  // Calculate an effective gamma
	T p = (q[2]-q[0]*std::pow(v, 2.0)/2.0)*(gammaeff-1.);  // Calculate pressure from ch10.pdf, eq 10.2

	return Vector4<T>(rho, u, p, q[3] / q[0]);
}


template <typename T>
Vector4<T> calcPhysicalFluxFromConservativeVec(Vector4<T> u) {
	T s1 = u[3] / u[0];
	T s2 = 1. - s1;
	T v = u[1] / u[0];

	T cp1 = 1005.9; T cp2 = 627.83;  // Cp values for air and sf6
	T cv1 = 717.09; T cv2 = 566.95;  // Cv values for air and sf6
	T gammaeff = (cp1*s1+cp2*s2) / (cv1*s1+cv2*s2);  // Calculate an effective gamma
	T P = (u[2]-u[0]*std::pow(v, 2.0)/2.0)*(gammaeff-1.);  // Calculate pressure from ch10.pdf, eq 10.2
	//	ret[0] = u[1]
	//	ret[1] = rho*np.power(v,2.0)+P;
	//	ret[2] = (u[2]+P)*v;
	//	ret[3] = v*u[3];
	//	return ret;

	return Vector4<T>(u[1], u[0]*v*v + P, (u[2]+P)*v, u[3]*v);
}


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
std::valarray<Vector4<T>> calcPhysFlux(const std::valarray<Vector4<T>>& u_arr) {
	/* Calculate physical fluxes of conservative variables
	 * in all points in the computation domain u_arr.
	 */

	return u_arr.apply(
				calcPhysicalFluxFromConservativeVec
			);  // f_arr
}


template <typename T>
T calcMaxWaveSpeedD(const std::valarray<Vector4<T>>& u_arr) {
	/* Calculate df/du. */

	std::size_t k = 0;

	T cp1 = 1005.9; T cp2 = 627.83; // Cp values for air and sf6
	T cv1 = 717.09; T cv2 = 566.95; // Cv values for air and sf6

	std::valarray<T> s1(u_arr.size());
	for (k = 0; k < u_arr.size(); ++ k)
		s1[k] = u_arr[k][3] / u_arr[k][0]; // s2 = W[4,:]/W[0,:]
	std::valarray<T> s2 = 1 - s1;

	std::valarray<T> arr_G = (s1*cp1+s2*cp2) / (s1*cv1+s2*cv2);
	// std::valarray<T> arr_G(1.4, U.size());
	std::valarray<T> a0(u_arr.size());
	for (k = 0; k < u_arr.size(); ++ k)
		a0[k] = arr_G[k] * u_arr[k][2] * (arr_G[k]-1.) / u_arr[k][0];
	a0 = std::sqrt(a0);

	return a0.max();
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

	T beta_0 = ((13./12.) * pow(f_prev2 - 2.*f_prev1 + f_curr0, 2)
				+ (1./4.) * pow(f_prev2 - 4.*f_prev1 + 3.*f_curr0, 2));

	T beta_1 = ((13./12.) * pow((f_prev1 - 2.*f_curr0 + f_next1), 2)
				+ (1/4.) * pow((f_prev1 - f_next1), 2));

	T beta_2 = ((13./12.) * pow(f_curr0 - 2.*f_next1 + f_next2, 2)
				+ (1./4.) * pow(3*f_curr0 - 4.*f_next1 + f_next2, 2));

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
T computeFHatWENO5FMReconstructionKernel(std::span<T, 5> f_stencil,
										 T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical fluxes
	 * `fhatplus` / `fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
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
//	alpha_weights = d_lin_weights / std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// Instead we use:
	std::valarray<T> alpha_weights = alphaWENO5FMWeights(
		std::move(beta_IS_coefs), eps, p
	);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
//	 omega_weights = alpha_weights / alpha_weights.sum();

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


template <typename T>
std::array<std::valarray<T>, 2> splitFluxAsLaxFriedrichs(
		const auto& u, const auto& f, T alpha) {
	/* Global Lax-Friedrichs (LF) flux splitting.
	 *
	 * For the purpose of linear stability (upwinding),
	 * a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	 * dfminus/du <= 0), is performed.
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
			std::span<T, 5>{std::begin(u_plus), 5}
		);

		std::reverse(std::begin(u_minus), std::end(u_minus));
		fhatminus = computeFHatWENO5FMReconstructionKernel(
			std::span<T, 5>{std::begin(u_minus), 5}
		);

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
void update_ghost_points(std::valarray<Vector4<T>>& U/*,
										std::size_t mini,
										std::size_t maxi*/) {
	/* Update ghost points in U. */

	//	for (k = nGhostCells; k < nGhostCells + nSize; ++ k) {
	//			U[k] = ms[k-nGhostCells].W;
	//	}

	const std::size_t nFullSize = U.size();
	const std::size_t nGhostCells = 3;  // 3
	const std::size_t mini = nGhostCells;
	const std::size_t maxi = nFullSize - nGhostCells - 1;
	// ... nFullSize-4   nFullSize-3   nFullSize-2   nFullSize-1 ]
	// ...    maxi
	// const std::size_t nSize = nFullSize - 2 * nGhostCells;

	// Transmissive b.c.s
	U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];
	U[maxi+1] = U[maxi]; U[maxi+2] = U[maxi]; U[maxi+3] = U[maxi];
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
	// dflux[3] -= lf[3];
	// dflux[ghost_point_number-1-3] += lf[ghost_point_number-1-3];

	// dflux += source terms...

	return dflux;
}


template <typename T>
void advanceTimestepTVDRK3(
	std::valarray<Vector4<T>>& U,
	std::valarray<Vector4<T>>& dflux,
	std::valarray<Vector4<T>>& Y2,
	std::valarray<Vector4<T>>& Y3,
	T t, T dt, T dx,
	const std::valarray<T>& lam, std::size_t nSize
) {
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

	// ----------------------First Stage--------------------
	// std::valarray<Vector4<T>> flux = calcFlux(U, t, lam);
	dflux = calcdSpace<T>(U, t, dx, lam, nSize);  // L1 = L[u^n]
	// std::valarray<Vector4<T>> res(Vector4<T>::ZERO, U.size());
	// L[u] = (-) dF[u]/dx

	Y2 = U + Vector4<T>(dt) * dflux;
	// u(1) = u^n + Δt L[u^n]

	update_ghost_points<T>(Y2);


	// ----------------------Second Stage-------------------
	dflux = calcdSpace<T>(Y2, t, dx, lam, nSize);  // L2 = L[u(1)]

	Y3 = Vector4<T>(3.)*U + Vector4<T>(dt) * dflux + Y2;
	Y3 *= Vector4<T>(0.25);
	// u(2) = 0.75 * u^n + 0.25 * u(1) + 0.25 * Δt L[u(1)]

	update_ghost_points<T>(Y3);


	// ----------------------Third Stage--------------------
	dflux = calcdSpace<T>(Y3, t, dx, lam, nSize);  // L3 = L[u(2)]


	U += Vector4<T>(2.) * (Y3 + Vector4<T>(dt) * dflux);
	U *= Vector4<T>(1./3.);
	// u^(n+1) = (1/3) * u^n + (2/3) * u(2) + (2/3) * Δt L[u(2)]

	update_ghost_points<T>(U);

	// return U;
	// U = res;
}


template <typename T>
void integrate(
	std::valarray<Vector4<T>>& U,
	std::valarray<Vector4<T>>& flux,
	std::valarray<Vector4<T>>& Y2,
	std::valarray<Vector4<T>>& Y3,
	T t0, T dx, std::size_t nSize,
	T fin_t, T cfl = 0.4
) {
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
}
