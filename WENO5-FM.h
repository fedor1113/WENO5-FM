#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <ranges>
#include <span>
#include <type_traits>
#include <valarray>
#include <vector>
#include <algorithm>

#include "_vector4.h"

template class Vector4<double>;

using Vec4 = Vector4<double>;

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


//double gete(double rho, double p) {
//	if (rho != 0.)
//		return p / (GAMMA-1.) / rho;

//	return 0.;
//}


//double getp(double rho, double e) {
//	return (GAMMA-1.) * rho * e;
//}


//double eFromConservative(double rho, double j, double rhoE) {
//	return (rhoE - 0.5 * j*j / rho) / rho;
//}


//Vec4 calcPhysicalFlux(double rho, double u, double p, double last) {
//	if (rho == 0) return Vec4::ZERO;


//	double e = gete(rho, p);
//	return Vec4(rho * u, p + rho*u*u,
//				u*(p + rho*(e + 0.5*u*u)),
//				u*last);
//}


//Vec4 conservativeToPrimitive(Vec4 q) {
//	// Primitive variables
//	double rho = q[0];
//	double u = q[1] / rho;
//	double E = q[2] / rho;
//	double p = (GAMMA - 1.) * rho * (E - 0.5*u*u);

//	return Vec4(rho, u, p, q[3]);
//}


//Vec4 calcPhysicalFluxFromConservativeVec(Vec4 u) {
////	return calcPhysicalFlux(u[0],
////			u[1] / u[0],
////			getp(u[0], eFromConservative(u[0], u[1], u[2])));
//	Vec4 prim = conservativeToPrimitive(u);

//	return calcPhysicalFlux(prim[0], prim[1], prim[2], prim[3]);
//}


Vec4 conservativeToPrimitive(Vec4 q) {
	// Primitive variables
	double rho = q[0];
	double u = q[1] / rho;
	// double E = q[2] / rho;


	double s1 = q[3] / q[0];
	double s2 = 1. - s1;
	double v = q[1] / q[0];
	double cp1 = 1005.9; double cp2 = 627.83; // Cp values for air and sf6
	double cv1 = 717.09; double cv2 = 566.95; // Cv values for air and sf6
	double gammaeff = (cp1*s1+cp2*s2) / (cv1*s1+cv2*s2); // Calculate an effective gamma
	double p = (q[2]-q[0]*std::pow(v, 2.0)/2.0)*(gammaeff-1.); // Calculate pressure from ch10.pdf, eq 10.2

	return Vec4(rho, u, p, q[3] / q[0]);
}


Vec4 calcPhysicalFluxFromConservativeVec(Vec4 u) {
	double s1 = u[3] / u[0];
	double s2 = 1. - s1;
	double v = u[1] / u[0];

	double cp1 = 1005.9; double cp2 = 627.83; // Cp values for air and sf6
	double cv1 = 717.09; double cv2 = 566.95; // Cv values for air and sf6
	double gammaeff = (cp1*s1+cp2*s2) / (cv1*s1+cv2*s2); // Calculate an effective gamma
	double P = (u[2]-u[0]*std::pow(v, 2.0)/2.0)*(gammaeff-1.); // Calculate pressure from ch10.pdf, eq 10.2
	//		ret[0] = u[1]
	//		ret[1] = rho*np.power(v,2.0)+P;
	//		ret[2] = (u[2]+P)*v;
	//		ret[3] = v*u[3];
	//		return ret;

	return Vec4(u[1], u[0]*v*v + P, (u[2]+P)*v, u[3]*v);
}


template <typename T1, typename T2>
std::valarray<double> vecMatDot(const T1& vec, const T2& mat) {
	/* Multiply a vector by a matrix (a sum product over
	 * the last axis of mat and vec).
	 */

	std::valarray<double> result(std::ranges::size(mat));

	for (size_t k = 0; k < std::ranges::size(mat)
			&& k < std::ranges::size(vec); ++ k)
		result[k] = std::inner_product(std::begin(mat[k]),
									   std::end(mat[k]),
									   std::begin(vec), 0.);

	return result;
}


std::valarray<Vec4> calcPhysFlux(std::valarray<Vec4> u_arr) {
	/* Calculate physical fluxes of conservative variables
	 * in all points in the computation domain u_arr.
	 */

	// return std::valarray(Vec4::ZERO, 500); // TO DO
	std::valarray<Vec4> f_arr = u_arr.apply(
		calcPhysicalFluxFromConservativeVec
	);

	return f_arr;
}


double calcMaxWaveSpeedD(std::valarray<Vec4> u_arr) {
	/* Calculate df/du. */

	size_t k = 0;

	double cp1 = 1005.9; double cp2 = 627.83; // Cp values for air and sf6
	double cv1 = 717.09; double cv2 = 566.95; // Cv values for air and sf6

	std::valarray<double> s1(0., u_arr.size());
	for (k = 0; k < u_arr.size(); ++ k)
		s1[k] = u_arr[k][3] / u_arr[k][0]; // s2 = W[4,:]/W[0,:]
	std::valarray<double> s2 = 1 - s1;

	std::valarray<double> arr_G = (s1*cp1+s2*cp2) / (s1*cv1+s2*cv2);
	// std::valarray<double> arr_G(1.4, U.size());
	std::valarray<double> a0(0., u_arr.size());
	for (k = 0; k < u_arr.size(); ++ k)
		a0[k] = arr_G[k] * u_arr[k][2] * (arr_G[k]-1.) / u_arr[k][0];
	a0 = std::sqrt(a0);

	return a0.max();
}


//template <typename T>
//std::array<double, 3> smoothness_indicators(const T& f_stensil) {
std::valarray<double> betaSmoothnessIndicators(
		std::span<double, 5> f_stensil) {
	/* Return the WENO smoothness indicators of Jiang and Shu (1996)
	 * for each of the 3 substancils.
	 * That is the sum of the normalized squares of the scaled
	 * L2-norms of all the derivatives of 3 local interpolating
	 * polynomials in the sub-stencils of 5-node `f_stensil`.
	 *
	 * This allows for (2*3-1)=5th order accuracy from the 3rd
	 * order Eno schemes.
	 */

	// std::array<double, 3> res;
	std::valarray<double> res(0., 3);
	double f_prev2 = f_stensil[0];
	double f_prev1 = f_stensil[1];
	double f_curr0 = f_stensil[2];
	double f_next1 = f_stensil[3];
	double f_next2 = f_stensil[4];

	double beta_0 = ((13./12.) * pow(f_prev2 - 2.*f_prev1 + f_curr0, 2)
					 + (1./4.) * pow(f_prev2 - 4.*f_prev1 + 3.*f_curr0, 2));

	double beta_1 = ((13./12.) * pow((f_prev1 - 2.*f_curr0 + f_next1), 2)
					 + (1/4.) * pow((f_prev1 - f_next1), 2));

	double beta_2 = ((13./12.) * pow(f_curr0 - 2.*f_next1 + f_next2, 2)
					 + (1./4.) * pow(3*f_curr0 - 4.*f_next1 + f_next2, 2));
	res[0] = beta_0;
	res[1] = beta_1;
	res[2] = beta_2;

	return res;
}


template <typename T>
std::valarray<double> betaSmoothnessIndicatorsMat(
	const T& f_stensil,
	std::array<std::array<const double, 6>, 6> _coefs[] = plus_coefs
) {
	/*Return the smoothness indicators beta_k, k=0,1,2
	 * for each of the 3 substancils of f_stensil.
	 */

	std::valarray<double> res(0., 3);
	std::valarray<double> prod_buf(std::ranges::size(f_stensil));

	for (size_t k = 0; k < 3; ++ k) {
		// for (size_t k = 0; k < half_size + 1; ++ k)
		prod_buf = vecMatDot(f_stensil, _coefs[k]);
		res[k] = std::inner_product(std::begin(prod_buf),
									std::end(prod_buf),
									std::begin(f_stensil), 0.);
	}


	return res;
}


double henrickGMappingForLambda(double lambda_weno_weight,
								double lambda_ideal = 1./3.) {
	/* The mapping function g by Henrick modified for symmetric
	 * lambda-weights by Zheng Hong, Zhengyin Ye and Kun Ye.
	 */

	double square_ideal = lambda_ideal * lambda_ideal;
	return lambda_weno_weight * (lambda_ideal
								 + square_ideal
								 - 3. * lambda_ideal * lambda_weno_weight
								 + lambda_weno_weight * lambda_weno_weight)
							  / (square_ideal
								 + lambda_weno_weight
								   * (1. - 2.*square_ideal));
}


// FD WENO5FM (WENO5-FM) - method
// Reconstruction based on LF flux splitting + improved mapped WENO
// of 5th order
// (see Mapped weighted essentially non-oscillatory schemes:
// achieving optimal order near critical points, 2005 by Henrick et al.)
// and 'An improved WENO-Z scheme with symmetry-preserving mapping'
// by Zheng Hong, Zhengyin Ye and Kun Ye, 2020
void calcHydroStageWENO5FM(const std::valarray<double>& u,
						   double t,
						   double lam,
						   std::valarray<double>& f,
						   size_t nSize,
						   double eps = 1e-40,
						   double p = 2.) {
	/* Component-wise finite-difference WENO5FM (WENO5-FM) - space
	 * reconstruction method with Lax-Friedrichs (LF) flux splitting.
	 *
	 * Usually, componentwise reconstruction produces satisfactory
	 * results for schemes up to third-order accuracy, while characteristic
	 * reconstruction produces better nonoscillatory results for
	 * higher-order accuracy, albeit with an increased computational cost.
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

	const unsigned order = 5;
	const size_t stensil_size = order;
	const size_t _actual_stensil_size = stensil_size + 1;  // 6
	const size_t half_size = order / 2;  // 2

	const size_t nGhostCells = (stensil_size + 1) / 2;  // r = (order + 1) / 2 = 3
	const size_t mini = nGhostCells;
	const size_t maxi = nGhostCells + nSize - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini + 2, maxi - 2};
	auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 2};
	// [g      g      g      i      i      i      i      i      i      ...]
	// {0      1      2      3      4      5}     6      7      8      ...
	//                |             |
	//                j        nGhostCells

	// WENO5 stencils

	// Coefficients WmN(+/-) before fluxes at the stensil nodes
	// to find component stensils
	// [q_k] = WmN(+/-) * [ f[j-2] ... f[j+2] ]
	// of the WENO interpolator.
	std::array<std::array<const double, 6>, 3> WmNminus {{
		{{0., 0., 0., 11./6, -7./6, 2./6}},
		{{0., 0., 2./6, 5./6, -1./6, 0.}},
		{{0., -1./6, 5./6, 2./6, 0., 0.}}
	}};

	std::array<std::array<const double, 6>, 3> WmNplus {{
		{{2./6, -7./6, 11./6, 0., 0., 0.}},
		{{0., -1./6, 5./6, 2./6, 0., 0.}},
		{{0., 0., 2./6, 5./6, -1./6, 0.}}
	}};

	// Calculation of f_hat, the numerical flux of u (whichever
	// is chosen), requires the approximation of u that uses at
	// the j-th cell of u-discretization a group of cell average
	// values on the left (u_minus) and on the right (u_plus).
	// `u_plus` represents the cells [j-2, j-1, j, j+1, j+2],
	// and `u_minus` represents the cells [j-1, j, j+1, j+2, j+3];
	// for convenience and uniformity we represent both using the
	// same combined structure of    [j-2, j-1, j, j+1, j+2, j+3].
	std::valarray<double> u_plus(0., _actual_stensil_size);   // f_plus
	std::valarray<double> u_minus(0., _actual_stensil_size);  // f_minus

	std::valarray<double> betaISplus(0., half_size + 1);
	std::valarray<double> betaISminus(0., half_size + 1);

	std::valarray<double> alphaplus(0., half_size + 1);
	std::valarray<double> alphaminus(0., half_size + 1);

	std::valarray<double> lambdaplus(0., half_size + 1);
	std::valarray<double> lambdaminus(0., half_size + 1);

	std::valarray<double> omegaplus(0., half_size + 1);
	std::valarray<double> omegaminus(0., half_size + 1);

	double fhatminus = 0.;
	double fhatplus = 0.;

	// std::valarray<double> matvecprod(stensil_size);
	std::valarray<double> matvecprod(half_size + 1);

	// The ideal weights (they generate the central upstream fifth-order
	// scheme for the 5-point stencil):
	std::valarray<double> d_lin_weights = {0.1, 0.6, 0.3};
	// From them WENO5-Z and WENO-M will calculate the non-linear
	// alpha and omega weights.

	// In WENO5-FM, further, we have one ideal value for λ
	// \overbar{lambda} = 1/3
	// double lambda_ideal = 1/3;
	// In the smooth region the smoothness indicators β_k ought to be equal
	// for all sub-stencils, and thus the weight's ideal value must be
	// unique.


	// For the purpose of linear stability (upwinding),
	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	// dfminus/du <= 0), is performed.
	// Lax-Friedrichs (LF) flux splitting (because it is the simplest
	// and smooth - we need the positive and negative fluxes to have
	// as many derivatives as the order of our finite-difference WENO)
	// is chosen here. `f_plus` uses a biased stencil with 1 point to
	// the left, while `f_minus` uses a biased stencil with 1 point to
	// the right.
	// double alpha = abs(U).max();
	// double alpha = std::sqrt(u.max().square());
	double alpha = std::abs(lam);  // α = max |df/du|
	// double alpha = std::abs(u).max();
//	std::valarray<Vec4> f_minus = Vec4(0.5) * (calcPhysFlux(U)
//						- Vec4(alpha) * U);
//	std::valarray<Vec4> f_plus = Vec4(0.5) * (calcPhysFlux(U)
//						+ Vec4(alpha) * U);
	std::valarray<double> f_plus = 0.5 * (f + alpha * u);
	std::valarray<double> f_minus = 0.5 * (f - alpha * u);

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
	for (size_t j = 0; j < u.size(); ++ j)
		f[j] = 0;

	auto j_it_p = std::begin(f_plus);
	auto j_it_m = std::begin(f_minus);

	for (auto j : shifted_index_range) {
	// for (size_t j = 10; j < 11; ++ j) {
		j_it_p = std::begin(f_plus);
		j_it_m = std::begin(f_minus);
		std::advance(j_it_p, j + half_size + 1 - stensil_size);
		std::advance(j_it_m, j + half_size + 1 - stensil_size);
		std::copy_n(j_it_p, u_plus.size(), std::begin(u_plus));
		std::copy_n(j_it_m, u_minus.size(), std::begin(u_minus));
//		std::copy_n(std::begin(f_plus)+j + half_size + 1 - stensil_size, u_plus.size(), std::begin(u_plus));
//		std::copy_n(std::begin(f_minus)+j + half_size + 1 - stensil_size, u_minus.size(), std::begin(u_minus));
//		std::cout << j + half_size + 1 - stensil_size << " " << j + half_size + 1 - stensil_size + u_plus.size() << "\n";
//		for (auto n : u_plus)//	for (size_t j = maxi+1; j < u.size(); ++ j) {
		//		// numerical_flux[j] = 0;
		//		f[j] = 0;
		//	}

		//	for (size_t j = 0; j < 2; ++ j) {
		//		// numerical_flux[j] = 0;
		//		f[j] = 0;
		//	}
//			std::cout << n << " ";
//		std::cout << "\n";
//		std::cout << j + half_size + 1 - stensil_size << " " << j + half_size + 1 - stensil_size + u_minus.size() << "\n";
//		for (auto n : u_minus)
//			std::cout << n << " ";
//		std::cout << "\n";

		// smoothness indicators of the stensil
		// (measure how smooth u is in the stensil)
//		betaISplus = betaSmoothnessIndicatorsMat(u_plus, plus_coefs);
//		betaISminus = betaSmoothnessIndicatorsMat(u_minus, minus_coefs);

		// The non-matrix variant seems to be faster(?)
		betaISplus = betaSmoothnessIndicators(
					std::span<double, 5>{std::begin(u_plus), 5});

		std::reverse(std::begin(u_minus), std::end(u_minus));
		betaISminus = betaSmoothnessIndicators(
					std::span<double, 5>{std::begin(u_minus), 5});
		// std::reverse(std::begin(u_minus), std::end(u_minus));

		// non-linear non-scaled (α-)weights
		//alphaplus = d_lin_weights / std::pow(eps + betaISplus, p);
		//alphaminus = d_lin_weights / std::pow(eps + betaISminus, p);
		// — we have no need of them in WENO-FM!

		// Instead we use:
		alphaplus = 1. / std::pow(eps + betaISplus, p);
		alphaminus = 1. / std::pow(eps + betaISminus, p);

		// scaled (normalized) non-linear (ω-)weights (ENO weights)
		// omegaplus = alphaplus / alphaplus.sum();
		// omegaminus = alphaminus / alphaminus.sum();

		// FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
		// due to Zheng Hong, Zhengyin Ye and Kun Ye:
		lambdaplus = alphaplus / alphaplus.sum();
		lambdaminus = alphaminus / alphaminus.sum();

		// And only to them a mapping in the spirit of Henrick et al.
		// is applied:
		auto gMap = [](double x) -> double {
			return henrickGMappingForLambda(x);
		};

		alphaplus = lambdaplus.apply(gMap);
		alphaminus = lambdaminus.apply(gMap);

		// From g(λ_k) and d_k we get the new corrected resultant
		// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
		omegaplus = d_lin_weights * alphaplus;
		omegaplus /= omegaplus.sum();

		omegaminus = d_lin_weights * alphaminus;
		omegaminus /= omegaminus.sum();


		// matvecprod stores an estimate of f_{i+1/2} for each substencil
		// which is then used to calculate
		// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
		matvecprod = vecMatDot(u_plus, WmNplus);
		fhatplus = std::inner_product(std::begin(omegaplus), std::end(omegaplus),
									  std::begin(matvecprod), 0.);
		// matvecprod = vecMatDot(u_minus, WmNminus);
		matvecprod = vecMatDot(u_minus, WmNplus);
		fhatminus = std::inner_product(std::begin(omegaminus), std::end(omegaminus),
									   std::begin(matvecprod), 0.);
		// std::reverse(std::begin(u_plus), std::end(u_minus));

		// numerical_flux[j] = fhatplus + fhatminus;
		f[j] = fhatplus + fhatminus;

//		for (auto n : f)
//			std::cout << n << " ";
//		std::cout << "\n";
	}

//	for (size_t j = maxi+1; j < u.size(); ++ j) {
//		// numerical_flux[j] = 0;
//		f[j] = 0;
//	}

//	for (size_t j = 0; j < 2; ++ j) {
//		// numerical_flux[j] = 0;
//		f[j] = 0;
//	}

	// return numerical_flux;
	// f = numerical_flux;
	// std::cout << " done!" << "\n";
}


void update_ghost_points(std::valarray<Vec4>& U/*,
										size_t mini,
										size_t maxi*/) {
	/* Update ghost points in U. */

	//	for (k = nGhostCells; k < nGhostCells + nSize; ++ k) {
	//			U[k] = ms[k-nGhostCells].W;
	//	}

	const size_t nFullSize = U.size();
	const size_t nGhostCells = 3;  // 3
	const size_t mini = nGhostCells;
	const size_t maxi = nFullSize - nGhostCells - 1;
	// ... nFullSize-4   nFullSize-3   nFullSize-2   nFullSize-1 ]
	// ...    maxi
	// const size_t nSize = nFullSize - 2 * nGhostCells;

	// Transmissive b.c.s
	U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];
	U[maxi+1] = U[maxi]; U[maxi+2] = U[maxi]; U[maxi+3] = U[maxi];
}


std::valarray<Vec4> calcFluxComponentWise(std::valarray<Vec4>& U,
							 double t, std::valarray<double>& lam,
							 size_t nSize,
							 double eps = 1e-40,
							 double p = 2.) {
	std::valarray<Vec4> res = calcPhysFlux(U);
	std::valarray<std::valarray<double>> components(
				std::valarray(0., U.size()), 4);
	std::valarray<std::valarray<double>> component_fs(
				std::valarray(0., U.size()), 4);

	size_t k = 0;

	for (size_t j = 0; j < 4; ++ j)
		for (k = 0; k < U.size(); ++ k) {
			components[j][k] = U[k][j];
			component_fs[j][k] = res[k][j];
		}

	k = 0;
	// for (auto component : components) {
	for (k = 0; k < 4; ++ k) {
		calcHydroStageWENO5FM(components[k], t, lam[k],
							  component_fs[k], nSize, eps, p);
		// ++ k;
	}

	for (size_t j = 0; j < 4; ++ j)
		for (k = 0; k < U.size(); ++ k)
			res[k][j] = component_fs[j][k];

	return res;
}


std::valarray<Vec4> calcdSpace(std::valarray<Vec4>& U,
							   double t,
							   double dx,
							   std::valarray<double>& lam,
							   size_t nSize,
							   double eps = 1e-40,
							   double p = 2.) {
	std::valarray<Vec4> dflux(Vec4::ZERO, U.size());
	std::valarray<Vec4> lf = calcFluxComponentWise(
				U, t, lam, nSize, eps, p
				);

	const size_t ghost_point_number = 3;

	std::slice Nweno(ghost_point_number, nSize, 1);
	std::slice Nweno_shifted_by_neg_1(ghost_point_number-1, nSize, 1);
//	std::slice Nweno(1, U.size()-1, 1);
//	std::slice Nweno_shifted_by_neg_1(0, U.size()-1, 1);

	std::valarray<Vec4> f_mn = lf[Nweno_shifted_by_neg_1];
	std::valarray<Vec4> f_pl = lf[Nweno];

	dflux[Nweno] = -(f_pl - f_mn) / Vec4(dx);
	// dflux[3] -= lf[3];
	// dflux[ghost_point_number-1-3] += lf[ghost_point_number-1-3];

	// dflux += source terms...

	return dflux;
}


void advanceTimestepTVDRK3(
	std::valarray<Vec4>& U,
	std::valarray<Vec4>& dflux,
	std::valarray<Vec4>& Y2,
	std::valarray<Vec4>& Y3,
	double t, double dt, double dx,
	std::valarray<double>& lam, size_t nSize
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
	// std::valarray<Vec4> flux = calcFlux(U, t, lam);
	dflux = calcdSpace(U, t, dx, lam, nSize);  // L1 = L[u^n]
	// std::valarray<Vec4> res(Vec4::ZERO, U.size());
	// L[u] = (-) dF[u]/dx

	Y2 = U + Vec4(dt) * dflux;
	// u(1) = u^n + Δt L[u^n]

	update_ghost_points(Y2);


	// ----------------------Second Stage-------------------
	dflux = calcdSpace(Y2, t, dx, lam, nSize);  // L2 = L[u(1)]

	Y3 = Vec4(3.)*U + Vec4(dt) * dflux + Y2;
	Y3 *= Vec4(0.25);
	// u(2) = 0.75 * u^n + 0.25 * u(1) + 0.25 * Δt L[u(1)]

	update_ghost_points(Y3);


	// ----------------------Third Stage--------------------
	dflux = calcdSpace(Y3, t, dx, lam, nSize);  // L3 = L[u(2)]


	U += Vec4(2.) * (Y3 + Vec4(dt) * dflux);
	U *= Vec4(1./3.);
	// u^(n+1) = (1/3) * u^n + (2/3) * u(2) + (2/3) * Δt L[u(2)]

	update_ghost_points(U);

	// return U;
	// U = res;
}


void integrate(
	std::valarray<Vec4>& U,
	std::valarray<Vec4>& flux,
	std::valarray<Vec4>& Y2,
	std::valarray<Vec4>& Y3,
	double t0, double dx, size_t nSize,
	double fin_t, double cfl = 0.4
) {
	double t = t0;

	double cpu = calcMaxWaveSpeedD(U);
	std::valarray<double> lam = std::valarray(cpu, 4);

	double dt = cfl * dx / lam[0];

	while (t < fin_t) {
		if (t + dt > fin_t)
			dt = fin_t - t;

		advanceTimestepTVDRK3(U, flux, Y2, Y3, t, dt, dx, lam, nSize);
		// std::cout << n*dt << "\n";
		// std::cout << t << "\n";

		t += dt;

		 cpu = calcMaxWaveSpeedD(U);
		 lam = std::valarray(cpu, 4);

		dt = cfl * dx / lam[0];
	}
	std::cout << "t = " << t << "\n";
}
