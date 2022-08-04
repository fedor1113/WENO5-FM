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
	 * the last axis of mat and vec). */

	std::valarray<double> result(std::ranges::size(mat));

	for (size_t k = 0; k < std::ranges::size(mat)
			&& k < std::ranges::size(vec); ++ k)
		result[k] = std::inner_product(std::begin(mat[k]),
									   std::end(mat[k]),
									   std::begin(vec), 0.);

	return result;
}


//Vect4 PrimitiveToConservative(Vect4 U) {
//	/*(rho, u, p, 0) -> (rho, j, e)*/
//	double E = eos.gete(res.ro, res.p) + .5*res.u*res.u;
//	Vector4 FGodunov = Vector4(res.ro*res.u, res.ro*res.u*res.u + res.p, res.u*(res.ro*E+res.p), 0.);
//	Vect4 res(U[0], U[0]*U[1], e_from_caloric_eos(U[0], U[2]), U[3])
//}


std::valarray<Vec4> calcPhysFlux(std::valarray<Vec4> u_arr) {
	/* Calculate physical fluxes of conservative variables
	 * in all points in the computation domain u_arr. */

	// return std::valarray(Vec4::ZERO, 500); // TO DO
	std::valarray<Vec4> f_arr = u_arr.apply(
		calcPhysicalFluxFromConservativeVec
	);

	return f_arr;
}

//template <typename T>
//std::array<double, 3> smoothness_indicators(const T& f_stensil) {
std::valarray<double> betaSmoothnessIndicators(
		std::span<double, 5> f_stensil) {
	/* Return the WENO smoothness indicators of Jiang and Shu (1996)
	 * for each of the 3 substancils.*/

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
	 * for each of the 3 substancils of f_stensil.*/

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
	 * lambda-weights by Zheng Hong, Zhengyin Ye and Kun Ye. */
	double square_ideal = lambda_ideal * lambda_ideal;
	return lambda_weno_weight * (lambda_ideal
								 + square_ideal
								 - 3. * lambda_ideal * lambda_weno_weight
								 + lambda_weno_weight * lambda_weno_weight)
							  / (square_ideal
								 + lambda_weno_weight * (1.-2.*square_ideal));
}


// WENO5ZM (WENO5-FM) - method
// Reconstruction based on LF flux + improved mapped WENO reconstruction of 5th order
// (see Mapped weighted essentially non-oscillatory schemes:
// achieving optimal order near critical points, 2005 by Henrick et al.)
// and 'An improved WENO-Z scheme with symmetry-preserving mapping'
// by Zheng Hong, Zhengyin Ye and Kun Ye, 2020


std::valarray<double> calcHydroStageWENO5ZM(const std::valarray<double>& u,
											double t,
											double lam,
											const std::valarray<double>& f,
											size_t nSize,
											double eps = 1e-40,
											double p = 2.) {
	/* Component-wise finite-difference WENO5-ZM (WENO5-FM) - space
	 * reconstruction method with Lax-Friedrichs (LF) flux splitting.
	 *
	 * Usually, componentwise reconstruction produces satisfactory
	 * results for schemes up to third-order accuracy, while characteristic
	 * reconstruction produces better nonoscillatory results for
	 * higher-order accuracy, albeit with an increased computational cost.
	*/

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-M);
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
	const size_t _actual_stensil_size = stensil_size + 1; // 6
	const size_t half_size = order / 2;  // 2
	// std::cout << "calcHydroStageWENO5M3G(): ";
	// EOSOld &eos = task.getEOS();
	// const double gamma = eos.getGamma();
	// double E = 0.;
	// int k = 0;
	// int nDimCounter = 0;
	// const int nSize = 500; // ms.getSize();
	// const size_t nSize = u.size();
	// const double h = ms[1].x-ms[0].x;

	// Vec4 L = Vec4::ZERO, R = Vec4::ZERO, D = Vec4::ZERO, V = Vec4::ZERO;
	// Special arrays keeping data (and flows) in extra fictitious cells (ghost points).
	const size_t nGhostCells = (stensil_size + 1) / 2;  // r = (order + 1) / 2 = 3
	const size_t mini = nGhostCells;
	const size_t maxi = nGhostCells + nSize;
	auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};

	// std::valarray<Vec4> U(Vec4::ZERO, nGhostCells + nSize + nGhostCells);
	// std::valarray<double> u(0., nGhostCells + nSize + nGhostCells);
	// std::valarray<Vec4> Up(Vec4::ZERO, nGhostCells + nSize + nGhostCells);
	// std::valarray<Vec4> Um(Vec4::ZERO, nGhostCells + nSize + nGhostCells);
	// std::valarray<Vec4> F(Vec4::ZERO, nGhostCells + nSize + nGhostCells);

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
	std::valarray<double> u_plus(_actual_stensil_size);  // f_plus
	std::valarray<double> u_minus(_actual_stensil_size); // f_minus

	std::valarray<double> betaISplus(half_size + 1);
	std::valarray<double> betaISminus(half_size + 1);

	std::valarray<double> alphaplus(half_size + 1);
	std::valarray<double> alphaminus(half_size + 1);

	std::valarray<double> lambdaplus(half_size + 1);
	std::valarray<double> lambdaminus(half_size + 1);

	std::valarray<double> omegaplus(half_size + 1);
	std::valarray<double> omegaminus(half_size + 1);

	double fhatminus = 0.;
	double fhatplus = 0.;

	// std::valarray<double> matvecprod(stensil_size);
	std::valarray<double> matvecprod(half_size + 1);

	// The ideal weights (they generate the central upstream fifth-order
	// scheme for the 5-point stencil):
	std::valarray<double> d_lin_weights = {0.1, 0.6, 0.3};
	// From them WENO5-Z and WENO-M will calculate the non-linear
	// alpha and omega weights.

	// In WENO5-ZM, further, we have one ideal value for λ
	// \overbar{lambda} = 1/3
	// double lambda_ideal = 1/3;
	// In the smooth region the smoothness indicators β_k ought to be equal
	// for all sub-stencils, and thus the weight's ideal value must be
	// unique.


	// For the purpose of linear stability (upwinding),
	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	// dfminus/du <= 0), is performed.
	// Lax-Friedrichs (LF) flux splitting (because it is the simplest)
	// is chosen here (for high-order WENO-schemes the exact choice of
	// flux splitting seems to be largely unimportant for the precision
	// of the scheme). `f_plus` uses a biased stencil with 1 point to
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
	// Riemann solver (see Toro, 2009).
	std::valarray<double> numerical_flux(u.size());

	for (auto j : shifted_index_range) {
		auto j_it_p = std::begin(f_plus);
		auto j_it_m = std::begin(f_minus);
		std::advance(j_it_p, j + half_size + 1 - stensil_size);
		std::advance(j_it_m, j + half_size + 1 - stensil_size);
		std::copy_n(j_it_p, u_plus.size(), std::begin(u_plus));
		std::copy_n(j_it_m, u_minus.size(), std::begin(u_minus));

		// smoothness indicators of the stensil
		// (measure how smooth u is in the stensil)
		betaISplus = betaSmoothnessIndicatorsMat(u_plus, plus_coefs);
		betaISminus = betaSmoothnessIndicatorsMat(u_minus, minus_coefs);

		// The non-matrix variant seems to be faster but is currently
		// non-usable: it contains some mistake that adds noticeable
		// artefacts...
//		auto rit = std::ranges::reverse_view{u_minus}.begin();
//		auto fit = rit.base()-1;
//		betaISplus = betaSmoothnessIndicators(std::span<double, 5>{std::begin(u_plus), 5});
//		betaISminus = betaSmoothnessIndicators(std::span<double, 5>{fit, 5});

		// non-linear non-scaled (α-)weights
		//alphaplus = d_lin_weights / std::pow(eps + betaISplus, p);
		//alphaminus = d_lin_weights / std::pow(eps + betaISminus, p);
		// — we have no need of them in WENO-ZM!

		// Instead we use:
		alphaplus = 1. / std::pow(eps + betaISplus, p);
		alphaminus = 1. / std::pow(eps + betaISminus, p);

		// scaled (normalized) non-linear (ω-)weights (ENO weights)
		// omegaplus = alphaplus / alphaplus.sum();
		// omegaminus = alphaminus / alphaminus.sum();

		// ZM-improved scaled (normalized) symmetric (λ-)weights
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
		// normalized WENO5-ZM (WENO5-FM) (ω_k-)weights:
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
		matvecprod = vecMatDot(u_minus, WmNminus);
		fhatminus = std::inner_product(std::begin(omegaminus), std::end(omegaminus),
									   std::begin(matvecprod), 0.);

		numerical_flux[j] = fhatplus + fhatminus;
	}

	return numerical_flux;
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
	const size_t maxi = nFullSize - nGhostCells;
	// const size_t nSize = nFullSize - 2 * nGhostCells;

	// Transmissive b.c.s
	U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];
	U[maxi] = U[maxi-1]; U[maxi+1] = U[maxi-1]; U[maxi+2] = U[maxi-1];
}


//void update_ghost_points2(std::valarray<Vec4>& U/*,
//										 size_t mini,
//										 size_t maxi*/) {
//	// Left hand cells will use a forward difference
//	// From: https://en.wikipedia.org/wiki/Finite_difference_coefficient#Forward_finite_difference
//	//	−137/60	5	−5	10/3	−5/4	1/5
//	//  dLdx = ((-137.0/60.0*y[:,2]+5.0*y[:,3]-5.0*y[:,4]+10.0/3.0*y[:,5]-5.0/4.0*y[:,6]+1.0/5.0*y[:,7]))/dx = 0
//	// Therefore, y[:,2] = -60.0/137.0*(-1.0/5.0*y[:,7]+5.0/4.0*y[:,6]-10.0/3.0*y[:,5]+5.0*y[:,4]-5.0*y[:,3])

//	const size_t nFullSize = U.size();
//	const size_t nGhostCells = 3;  // 3
//	const size_t mini = nGhostCells;
//	const size_t maxi = nFullSize - nGhostCells;

//	std::valarray<Vec4> u_new = U;

//	// for (size_t k = 0; k < nFullSize; ++ k) {
//	u_new[2] = (Vec4(-60./137.) * (Vec4(-1./5.) * U[7]
//			+ Vec4(5./4.) * U[6]-Vec4(10./3.) * U[5]+Vec4(5.) * U[4]
//			- Vec4(5.) * U[3]));
//	u_new[2][3] = U[2][3];
//	u_new[2][1] = 0.;
//	u_new[1] = u_new[2];
//	u_new[0] = u_new[2];
//	// Same idea for the backwards differences, except the signs of the coefficients are flipped
//	U[nFullSize-3] = (Vec4(60./137.)
//			* (Vec4(1./5.) * U[nFullSize-8] - Vec4(5./4.) * U[nFullSize-7]
//			+ Vec4(10./3.) * U[nFullSize-6] - Vec4(5.) * U[nFullSize-5]
//			+ Vec4(5.) * U[nFullSize-4]));
//	u_new[nFullSize-3][3] = U[nFullSize-3][3];
//	u_new[nFullSize-3][1] = 0.;
//	u_new[nFullSize-2] = u_new[nFullSize-3];
//	u_new[nFullSize-1] = u_new[nFullSize-3];
//	// Note that this just does a "zero-th order" extrapolation in to the ghost cells,
//	// and so does not conserve mass when characteristics start hitting the boundaries
//	// http://physics.princeton.edu/~fpretori/Burgers/Boundary.html <- useful link for extrapolation

//	U = u_new;
//}

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
		component_fs[k] = calcHydroStageWENO5ZM(components[k], t, lam[k],
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
							   double eps = 1e-20,
							   double p = 2.) {
	std::valarray<Vec4> dflux(Vec4::ZERO, U.size());
	std::valarray<Vec4> lf = calcFluxComponentWise(
				U, t, lam, nSize, eps, p
				);

	std::slice Nweno(3, nSize, 1);
	std::slice Nweno_shifted_by_neg_1(3-1, nSize, 1);

	std::valarray<Vec4> f_mn = lf[Nweno_shifted_by_neg_1];
	std::valarray<Vec4> f_pl = lf[Nweno];

	dflux[Nweno] = -(f_pl - f_mn) / Vec4(dx);
	// dflux[3+1] -= lf[3+1] / Vec4(dx);
	// dflux[nSize-1] += lf[nSize-1] / Vec4(dx);

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
	 * du/dt = L[u].
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
	dflux = calcdSpace(Y3, t, dx, lam, nSize); // L3 = L[u(2)]

	U *= Vec4(1./3.);
	U += Vec4(2./3.) * (Y3 + Vec4(dt) * dflux);
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
	size_t k = 0;

	double cp1 = 1005.9; double cp2 = 627.83; // Cp values for air and sf6
	double cv1 = 717.09; double cv2 = 566.95; // Cv values for air and sf6

	std::valarray<double> s1 = std::valarray(0., U.size());
	for (k = 0; k < U.size(); ++ k)
		s1[k] = U[k][3] / U[k][0]; // s2 = W[4,:]/W[0,:]
	std::valarray<double> s2 = 1 - s1;

	double t = t0;
	std::valarray<double> arr_G = (s1*cp1+s2*cp2) / (s1*cv1+s2*cv2);
	// std::valarray<double> arr_G(1.4, U.size());

	std::valarray<double> a0 = std::valarray(0., U.size());
	for (k = 0; k < U.size(); ++ k)
		a0[k] = arr_G[k] * U[k][2] * (arr_G[k]-1.) / U[k][0];
	a0 = std::sqrt(a0);
	std::valarray<double> lam = std::valarray(a0.max(), 4);

	double dt = cfl * dx / lam[0];

	while (t < fin_t) {
		if (t + dt > fin_t)
			dt = fin_t - t;

		advanceTimestepTVDRK3(U, flux, Y2, Y3, t, dt, dx, lam, nSize);
		// std::cout << n*dt << "\n";
		// std::cout << t << "\n";

		t += dt;

		for (k = 0; k < U.size(); ++ k)
			a0[k] = arr_G[k] * U[k][2] * (arr_G[k]-1.) / U[k][0];
		a0 = std::sqrt(a0);
		lam = std::valarray(a0.max(), 4);

		dt = cfl * dx / lam[0];
	}
	std::cout << "t = " << t << "\n";
}
