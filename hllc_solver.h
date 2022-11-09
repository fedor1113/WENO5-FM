#ifndef HLLC_SOLVER_H
#define HLLC_SOLVER_H

#include <algorithm>
#include <cstdlib>
#include <execution>
#include <iostream>
#include <ranges>

#include "arithmeticwith.h"

#include "euler1d.h"

//template <ArithmeticWith<numeric_val> T>
//T getP(T rho, T e, T gamma);

//template <ArithmeticWith<numeric_val> T>
//T calcSquareSoundSpeed(T rho, T rho_v, T rho_E, T gamma/* = 1.4*/);


//template <ArithmeticWith<numeric_val> T>
//Vector4<T> calcHLLCFlux(
//		T roL, T rouL, T roEL,
//		T roR, T rouR, T roER, T gamma = 1.4) {
//	T _ro = 0., _u = 0., _e = 0., _p = 0.;
//	T uL = rouL/roL, uR = rouR/roR;
//	T EL = roEL/roL, ER = roER/roR;
//	T eL = roEL/roL - .5*uL*uL, eR = roER/roR - .5*uR*uR;
//	if (roL == 0. || roR == 0.) {
//		std::cout << "Error: CSolver::calcHLLFluxEOSIdeal(): vacuum is present." << "\n";
//		std::exit(EXIT_FAILURE);
//	}
//	T pL = getP(roL, eL, gamma), pR = getP(roR, eR, gamma);
//	T HL = EL + pL/roL, HR = ER + pR/roR;
//	T cL = std::sqrt(calcSquareSoundSpeed(roL, rouL, roEL, gamma));
//	T cR = std::sqrt(calcSquareSoundSpeed(roR, rouR, roER, gamma));
//	Vector4<T> UL = Vector4<T>(roL, rouL, roEL, 0.), UR = Vector4<T>(roR, rouR, roER, 0.);
//	_ro = roL, _u = rouL/roL, _e = roEL/roL-.5*_u*_u, _p = getP(_ro, _e, gamma);
//	Vector4<T> FL = Vector4<T>(rouL, _p + _ro*_u*_u, _u*(_p + roEL), 0.);
//	_ro = roR, _u = rouR/roR, _e = roER/roR-.5*_u*_u, _p = getP(_ro, _e, gamma);
//	Vector4<T> FR = Vector4<T>(rouR, _p + _ro*_u*_u, _u*(_p + roER), 0.);
//	// Step 1 -- pressure estimate from primitive-variable Riemann solver (PVRS)
//	T roAv = .5*(roL + roR), cAv = .5*(cL + cR);
//	T pPVRS = .5*(pL + pR) - .5*(uR - uL)*roAv*cAv;
//	//double pStar = max(0., pPVRS);
//	T pStar = std::max(static_cast<T>(-15.e9), pPVRS);


//	// Step 2 -- wave speed estimates
//	// Uncomment for Roe averaging. Works quite well for Mie-Gruneisen EOS
//	// but fails on Toro test 1 at rarefaction wave (imposes the discontinuity)
//	T _roAv = sqrt(roL*roR);
//	T _uAv = (uL*sqrt(roL) + uR*sqrt(roR))/(sqrt(roL)+sqrt(roR));
//	T _HAv = (HL*sqrt(roL) + HR*sqrt(roR))/(sqrt(roL)+sqrt(roR));
//	T _EAv = (EL*sqrt(roL) + ER*sqrt(roR))/(sqrt(roL)+sqrt(roR));
//	T _pAv = (_HAv - _EAv) * _roAv;
//	T _cAv = std::sqrt(
//				calcSquareSoundSpeed(_roAv, _roAv * _uAv, _roAv * _EAv, gamma));

//	T SL = _uAv-_cAv, SR = _uAv+_cAv;

//	T SStar = (pR - pL + roL*uL*(SL - uL)
//					- roR*uR*(SR - uR))/(roL*(SL - uL) - roR*(SR - uR));
//	// Vector4 Uhll = (SR*UR - SL*UL + FL - FR)/(SR - SL);
//	Vector4<T> Fhllc = Vector4<T>::ZERO;
//	Vector4<T> D = Vector4<T>(0., 1., SStar, 0.);
//	Vector4<T> UStarL = (SL*UL - FL + pStar*D)/(SL - SStar);
//	Vector4<T> UStarR = (SR*UR - FR + pStar*D)/(SR - SStar);
//	Vector4<T> FStarL = FL + SL*(UStarL - UL), FStarR = FR + SR*(UStarR - UR);
//	if (0. <= SL)
//		Fhllc = FL;
//	else if (SL <= 0. && 0 <= SStar)
//		Fhllc = FStarL;
//	else if (SStar <= 0. && 0 <= SR)
//		Fhllc = FStarR;
//	else
//		Fhllc = FR;

//	return Fhllc;
//}


//template <ArithmeticWith<numeric_val> T>
//void calcHLLCFlux(
//		const std::ranges::common_range auto& u_p,
//		const std::ranges::common_range auto& u_m,
//		std::ranges::common_range auto& res_f,
//		T gamma = 1.4) {
//	/* Global Lax-Friedrichs (LF) numerical flux. LF flux is
//	 * attractive because it is extremely simple and easy to
//	 * compute, but it is also probably the most dissipative
//	 * of all the well-known fluxes, smearing the discontinuities
//	 * quite noticeably for lower-order reconstruction methods.
//	 *
//	 * For high-order reconstructions this drawback becomes much
//	 * less noticeable.
//	 *
//	 * P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations
//	 * and Their Numerical Computation, Commun. Pure and Applied
//	 * Mathematics, 7, 159-193, 1954.
//	 */

//	std::transform(std::ranges::begin(u_p)+5, std::ranges::end(u_p)-5,
//				std::ranges::begin(u_m)+5,
//				std::ranges::begin(res_f)+5,
//				[gamma](auto u_l, auto u_r) {
//		return calcHLLCFlux(u_l[0], u_l[1], u_l[2],
//							u_r[0], u_r[1], u_r[2], gamma);
//	});
//}

#endif // HLLC_SOLVER_H
