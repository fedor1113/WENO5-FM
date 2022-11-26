#ifndef ROE_H
#define ROE_H

//#include <algorithm>
//#include <cmath>
//#include <cstdlib>
//#include <execution>
//#include <iostream>
//#include <ranges>

//#include "arithmeticwith.h"

//#include "_vector4.h"
//#include "eos.h"

//template <ArithmeticWith<numeric_val> T>
//T getP(T rho, T e, T gamma);


//template <ArithmeticWith<numeric_val> T>
//T calcSquareSoundSpeed(T rho, T rho_v, T rho_E, T gamma/* = 1.4*/);


//template <ArithmeticWith<numeric_val> T>
//Vector4<T> calcRoeFlux_(
//		T rhoL, T rhouL, T rhoEL, T rhoR, T rhouR, T rhoER, T gamma) {
//	T uL = rhouL / rhoL;
//	T eL = rhoEL / rhoL - 0.5 * uL * uL;
//	T uR = rhouR / rhoR;
//	T eR = rhoER / rhoR - 0.5 * uR * uR;

//	T pL = getP(rhoL, eL, gamma);
//	T pR = getP(rhoR, eR, gamma);
//	T cL = std::sqrt(calcSquareSoundSpeed(rhoL, rhouL, rhoEL, gamma));
//	T cR = std::sqrt(calcSquareSoundSpeed(rhoR, rhouR, rhoER, gamma));

//	T EL = rhoEL / rhoL;
//	T ER = rhoER / rhoR;
//	T HL = (rhoEL + pL) / rhoL;
//	T HR = (rhoER + pR) / rhoR;

//	T rhoAv = std::sqrt(rhoL * rhoR);
//	T sqrhoL = std::sqrt(rhoL);
//	T sqrhoR = std::sqrt(rhoR);
//	T uAv = (sqrhoL * uL + sqrhoR * uR) / (sqrhoL + sqrhoR);
//	T HAv = (sqrhoL * HL + sqrhoR * HR) / (sqrhoL + sqrhoR);
//	T EAv = (sqrhoL * EL + sqrhoR * ER) / (sqrhoL + sqrhoR);
//	T pAv = (HAv - EAv) * rhoAv;
//	T cAv = std::sqrt(
//				calcSquareSoundSpeed(
//					rhoAv, rhoAv * uAv, rhoAv * EAv, gamma));
//	// Ideal version to compare:
//	/*double _gam = 5./3.;
//	cAv = std::sqrt((HAv-.5*uAv*uAv)*(_gam-1.));*/
//	Vector4<T> Lambda = Vector4<T>(uAv - cAv, uAv, uAv + cAv, 0.);
//	Vector4<T> K0 = Vector4<T>(1., uAv-cAv, HAv-uAv*cAv, 0.);
//	Vector4<T> K1 = Vector4<T>(1., uAv,     uAv*uAv/2.,  0.);
//	Vector4<T> K2 = Vector4<T>(1., uAv+cAv, HAv+uAv*cAv, 0.);
//	Vector4<T> Alpha = Vector4<T>::ZERO;
//	// For ideal gas:
//	//Alpha[1] = (gamma-1.)/cAv/cAv * (
//	//			(HAv-uAv*uAv)*(rhoR-rhoL)
//	//			+ uAv*(rhouR-rhouL) - (rhoER-rhoEL));
//	// Using cAv2 = (gamma-1) * (HAv-.5*uAv2), (gamma-1)/cAv2 = 1/(HAv-.5uAv2)
//	Alpha[1] = 1./(HAv-.5*uAv*uAv) * (
//				(HAv-uAv*uAv)*(rhoR-rhoL)
//				+ uAv*(rhouR-rhouL)
//				- (rhoER-rhoEL));
//	Alpha[0] = 1./2./cAv * (
//				(rhoR-rhoL)*(uAv+cAv)
//				- (rhouR-rhouL)
//				- cAv*Alpha[1]);
//	Alpha[2] = (rhoR - rhoL) - Alpha[0] - Alpha[1];
//	Vector4<T> FL = Vector4<T>(
//				rhoL * uL, pL + rhoL * uL * uL, uL * (pL + rhoEL), 0.);
//	Vector4<T> FR = Vector4<T>(
//				rhoR * uR, pR + rhoR * uR * uR, uR * (pR + rhoER), 0.);
//	T delta = 0.;
//	delta = std::max(delta, 4. * ((uR - cR) - (uL - cL)));
//	if (std::abs(Lambda[0]) < 0.5 * delta)
//		Lambda[0] = Lambda[0] * Lambda[0] / delta + .25 * delta;
//	delta = std::max(static_cast<T>(0.), 4. * ((uR + cR) - (uL + cL)));
//	if (std::abs(Lambda[2]) < 0.5 * delta)
//		Lambda[2] = Lambda[2] * Lambda[2] / delta + .25 * delta;
//	Vector4<T> FRoe = 0.5 * (FL + FR) - 0.5 * (
//				Alpha[0] * std::abs(Lambda[0]) * K0
//				+ Alpha[1] * std::abs(Lambda[1]) * K1
//				+ Alpha[2] * std::abs(Lambda[2]) * K2);

//	return FRoe;
//}


////template <ArithmeticWith<numeric_val> T>
////Vector4<T> calcMarquinaFlux_(Vector4<T> left, Vector4<T> right, T gamma) {
////	Vector4<T> char_u_left = projectOntoCharacteristics(
////				left, left, gamma);
////	Vector4<T> char_f_left = projectOntoCharacteristics(
////				left,
////				calcPhysicalFluxFromConservativeVec<T>(right, gamma),
////				gamma);
////	Vector4<T> char_u_right = projectOntoCharacteristics(
////				right, right, gamma);
////	Vector4<T> char_f_right = projectOntoCharacteristics(
////				right,
////				calcPhysicalFluxFromConservativeVec<T>(right, gamma),
////				gamma);
////}


//template <ArithmeticWith<numeric_val> T>
//void calcRoeFlux(
//		const std::ranges::common_range auto& u_p,
//		const std::ranges::common_range auto& u_m,
//		std::ranges::common_range auto& res_f,
//		T gamma) {
//	std::transform(std::ranges::begin(u_p), std::ranges::end(u_p),
//				std::ranges::begin(u_m),
//				std::ranges::begin(res_f),
//				[gamma](auto u_left, auto u_right) {
//		return calcRoeFlux_(
//					u_left[0], u_left[1], u_left[2],
//					u_right[0], u_right[1], u_right[2],
//					gamma);
//	});
//}

#endif // ROE_H
