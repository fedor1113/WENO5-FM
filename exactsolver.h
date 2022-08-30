/*
-----------------------------------------------------------------------------
Adapted from fsla: https://github.com/vadimvshepelev/fsla
-----------------------------------------------------------------------------
*/


#ifndef EXACTSOLVER_H
#define EXACTSOLVER_H

#include <cmath>

#include "_vector4.h"
#include "arithmeticwith.h"
#include "euler1d.h"


// User-defined data-type for structure of Riemann problem solution
enum RPSolutionType {
	SWRW, RWSW, SWSW, RWRW,
	VacRW, RWVac, RWVacRW, nowaves
};


// Structure for the Riemann problem solution vector of primitive
// variables -- (roL, roR, v ,p)
template <ArithmeticWith<numeric_val> T>
struct RPSolutionPrimitive {
	RPSolutionPrimitive()
		: roL(0.), roR(0.), v(0.), p(0.), type(RPSolutionType::nowaves) {}
	T roL, roR, v, p;
	RPSolutionType type;
};


template <ArithmeticWith<numeric_val> T>
struct CVectorPrimitive {
	T ro;
	T v;
	T p;
};


enum RPWaveConfig {nothing, swrw, rwsw, swsw, rwrw, vacrw, rwvac, rwvacrw};


// Structure for the Riemann problem solution vector of primitive
// variables -- (roL, roR, u ,p)
template <ArithmeticWith<numeric_val> T>
struct RPValues {
	RPValues()
		: roL(0.), roR(0.), u(0.), p(0.), type(RPWaveConfig::nothing) {}
	T roL, roR, u, p;
	RPWaveConfig type;
};


template <ArithmeticWith<numeric_val> T>
struct C1DVectorPrimitive {
	T ro;
	T u;
	T p;
};


//template <ArithmeticWith<numeric_val> T>
//T gete(T rho, T p, T gamma = 1.4) {
//	if (rho != 0.)
//		return p / (gamma - 1.) / rho;

//	return 0.;
//}


//template <ArithmeticWith<numeric_val> T>
//T getP(T rho, T e, T gamma = 1.4) {
//	return (gamma - 1.) * rho * e;
//}


template <ArithmeticWith<numeric_val> T>
T gete(T, T, T);


template <ArithmeticWith<numeric_val> T>
T getP(T, T, T);


template <ArithmeticWith<numeric_val> T>
T getc(T rho, T p, T gamma = 1.4) {
	if (rho != 0.)
		return std::sqrt(gamma * p / rho);

	return 0.;
}


template <ArithmeticWith<numeric_val> T>
T fL(T p, T roL, T uL, T pL) {
	const T gamma = getc(1., 1.) * getc(1., 1.);
	T f = 0.;
	if (p > pL) {
		T AL = 2. / (gamma + 1.) / roL;
		T BL = (gamma - 1.) / (gamma + 1.) * pL;
		f = (p - pL) * std::sqrt(AL / (p + BL));
		return f;
	} else {
		T cL = std::sqrt(gamma * pL / roL);
		f = 2. * cL / (gamma - 1.) * (
					(pow(p / pL, (gamma - 1.) / 2. / gamma)) - 1.);
		return f;
	}
}


template <ArithmeticWith<numeric_val> T>
T dfLdp(T p, T roL, T vL, T pL) {
	const double gamma = getc(1., 1.) * getc(1., 1.);
	double dfdp = 0.;
	if (p > pL) {
		T AL = 2. / (gamma + 1.) / roL;
		T BL = (gamma - 1.) / (gamma + 1.) * pL;
		dfdp = std::sqrt(
					AL / (p + BL)) * (1. - (p - pL) / 2. /(p + BL));
		return dfdp;
	}
	else {
		T cL = std::sqrt(gamma * pL / roL);
		dfdp = cL / pL / gamma * std::pow(
					p / pL, -(gamma + 1) / 2. / gamma);
		return dfdp;
	}
}


template <ArithmeticWith<numeric_val> T>
T fR(T p, T roR, T vR, T pR) {
	const T gamma = getc(1., 1.) * getc(1., 1.);
	T f = 0.;
	if (p > pR) {
		T AR = 2. / (gamma + 1.) / roR;
		T BR = (gamma - 1.)/(gamma + 1.) * pR;
		f = (p - pR) * std::sqrt(AR / (p + BR));
		return f;
	} else {
		T cR = std::sqrt(gamma * pR / roR);
		f = 2. * cR / (gamma - 1.) * (
					(std::pow(p / pR, (gamma - 1.) / 2. / gamma)) - 1.);
		return f;
	}
}


template <ArithmeticWith<numeric_val> T>
T dfRdp(T p, T roR, T vR, T pR) {
	const T gamma = getc(1., 1.) * getc(1.,1.);
	T dfdp = 0.;
	if (p > pR) {
		T AR = 2. / (gamma + 1.) / roR;
		T BR = (gamma - 1.) / (gamma + 1.) * pR;
		dfdp = std::sqrt(AR / (p + BR)) * (1. - (p-pR) / 2. / (p+BR));
		return dfdp;
	} else {
		T cR = std::sqrt(gamma * pR / roR);
		dfdp = cR / pR / gamma * std::pow(
					p / pR, -(gamma + 1) / 2. / gamma);
		return dfdp;
	}
}


// Function finds resulting ro, u, p and wave structure of Riemann problem
template <ArithmeticWith<numeric_val> T>
RPValues<T> calcValues(T roL, T uL, T pL,
					   T roR, T uR, T pR) {
	// Решаем нелинейное уравнение относительно давления
	// методом касательных Ньютона
	RPValues<T> res;
	T p = 0., pPrev = 0.;
	const T TOL = 1.e-6;
	const T gamma = getc(1., 1.) * getc(1., 1.);
	int itCounter = 0;
	T cL = 0., cR = 0.;
	if (roL != 0.) cL = getc(roL, pL);
	if (roR != 0.) cR = getc(roR, pR);;
	// Пытаюсь определить возможную конфигурацию решения,
	// чтобы вернее выставить начальное приближение
	// Похоже, итерации нужны только в случаях "УВ+УВ" и "УВ + ВР",
	// т.к. в случае ВР+ВР и ВР+вакуум есть
	// аналитические решения для идеального газа
	//
	// Также вызывает вопрос последний тест Торо,
	// где полученное решение отличается от его решения
	// во втором знаке после запятой
	if (roL == roR && uL == uR && pL == pR) {
		res.type = RPWaveConfig::rwrw;
		res.roL  = roL;
		res.roR  = roL;
		res.p    = pL;
		res.u	 = uL;
		return res;
	}
	if (roL == 0.) {
		res.type = RPWaveConfig::vacrw;
		res.roL  = 0.;
		res.roR  = roR;
		res.p	 = 0.;
		res.u    = uR - 2. * cR / (gamma - 1.);
		return res;
	}
	if (roR == 0.) {
		res.type = RPWaveConfig::rwvac;
		res.roL  = roL;
		res.roR  = 0.;
		res.p	 = 0.;
		res.u    = uL + 2. * cL / (gamma - 1.);
		return res;
	}
	if (2. * cL / (gamma - 1.) + 2. * cR / (gamma - 1.)
			< std::abs(uL - uR)){
		res.type  = RPWaveConfig::rwvacrw;
		res.roL	  = 0.;
		res.roR   = 0.;
		res.u     = 0.;
		res.p	  = 0.;
		return res;
	}

	double fLmin = fL(pL, roL, uL, pL)
			+ fR(pL, roR, uR, pR) + uR - uL;
	double fRMax = fL(pR, roL, uL, pL)
			+ fR(pR, roR, uR, pR) + uR - uL;
	// Начальное приближение
	//p = 0.5*(pL+pR);
	p = pL / 2.;
	do {
		pPrev = p;
		p = pPrev - (fL(pPrev, roL, uL, pL)
					 + fR(pPrev, roR, uR, pR) + uR - uL)
				/ (dfLdp(pPrev, roL, uL, pL)
				   + dfRdp(pPrev, roR, uR, pR));
		if (p <= 0.)
			p = TOL;
		++ itCounter;
	} while (std::abs(2. * (p - pPrev) / (p + pPrev)) > TOL);
	res.p   = p;
	res.u   = 0.5 * (uL + uR)
			+ 0.5 * (fR(p, roR, uR, pR) - fL(p, roL, uL, pL));
	if (p < pL && p > pR) {
		res.type = RPWaveConfig::rwsw;
		res.roL  = roL * std::pow(res.p / pL, 1. / gamma);
		res.roR  = roR * (res.p / pR + (gamma - 1.) / (gamma + 1.))
				/ ((gamma - 1.) / (gamma + 1.) * res.p / pR + 1.);
	} else if (p <= pL && p <= pR) {
		res.type = RPWaveConfig::swsw;
		res.roL  = roL * std::pow(res.p / pL, 1. / gamma);
		res.roR  = roR * std::pow(res.p / pR, 1. / gamma);
	} else if (p > pL && p < pR) {
		res.type = RPWaveConfig::swsw;
		res.roL  = roL*(res.p / pL + (gamma - 1.) / (gamma + 1.))
				/((gamma - 1.) / (gamma + 1.) * res.p / pL + 1.);
		res.roR  = roR * std::pow(res.p / pR, 1. / gamma);
	} else {
		res.type = RPWaveConfig::swsw;
		res.roL  = roL * (res.p / pL + (gamma - 1.) / (gamma+1.))
				/ ((gamma-1.) / (gamma+1.) * res.p/pL + 1.);
		res.roR  = roR * (res.p / pR + (gamma - 1.) / (gamma+1.))
				/ ((gamma-1.) / (gamma+1.) * res.p/pR + 1.);
	}

	return res;
}


template <ArithmeticWith<numeric_val> T>
C1DVectorPrimitive<T> calcSolution(
		T roL, T uL, T pL,
		T roR, T uR, T pR,
		T x, T t) {
	RPValues<T> res = calcValues(roL, uL, pL, roR, uR, pR);
	// V = (ro, v, p)T
	C1DVectorPrimitive<T> V;
	T xi = x / t;
	const T gamma = getc(1., 1.) * getc(1., 1.);
	T cL = 0., cR = 0.;
	if (roL != 0.) cL = getc(roL, pL);
	if (roR != 0.) cR = getc(roR, pR);;
	T xiFront=0., xiHead=0., xiTail=0.,
			xiHeadL=0., xiTailL=0.,
			xiHeadR=0., xiTailR=0.;
	// Если вакуум
	if (res.type == RPWaveConfig::vacrw) {
		xiHead = uR + cR;
		xiTail = uR - 2. * cR / (gamma - 1.);
		if (xi <= xiTail) {
			V.ro = 0.;
			V.u  = uR - 2. * cR / (gamma - 1.);
			V.p  = 0.;
		} else if (xi < xiHead) {
			V.ro = roR * std::pow(
						2. / (gamma + 1.) - (gamma - 1.) / (gamma + 1.)
							/ cR * (uR - xi),
						2. / (gamma - 1.));
			V.u  = 2. / (gamma + 1) * (-cR + (gamma - 1.) / 2. * uR + xi);
			V.p  = pR * std::pow(
						2. / (gamma + 1.) - (gamma - 1.)
							/ (gamma + 1.) / cR * (uR - xi),
						2. * gamma / (gamma - 1.));
		} else {
			V.ro = roR;
			V.u  = uR;
			V.p  = pR;
		}
		return V;
	}
	if (res.type == RPWaveConfig::rwvac) {
		xiHead = uL - cL;
		xiTail = uL + 2. * cL / (gamma - 1.);
		if(xi >= xiTail) {
			V.ro = 0.;
			V.u  = 0.;
			V.p  = 0.;
		} else if (xi > xiHead) {
			V.ro = roL * std::pow(
						2. / (gamma + 1.) + (gamma - 1.)
							/ (gamma + 1.) / cL * (uL - xi),
						2. / (gamma - 1.));
			V.u  = 2./(gamma+1)*(cL + (gamma-1.)/2.*uL + xi);
			V.p  = pL * std::pow(
						2. / (gamma + 1.) + (gamma - 1.)
							/ (gamma + 1.) / cL * (uL - xi),
						2. * gamma / (gamma - 1.));
		} else {
			V.ro = roL;
			V.u = uL;
			V.p = pL;
		}
		return V;
	}
	if (res.type ==  RPWaveConfig::rwvacrw) {
		xiHeadL = uL - cL;
		xiTailL = uL + 2. * cL / (gamma - 1.);
		xiHeadR = uR + cR;
		xiTailR = uR - 2. * cR / (gamma - 1.);
		if (xi <= xiHeadL) {
			V.ro = roL;
			V.u  = uL;
			V.p  = pL;
		} else if (xi < xiTailL) {
			V.ro = roL * std::pow(
						2. / (gamma + 1.) + (gamma - 1.)
							/ (gamma + 1.) / cL * (uL - xi),
						2. / (gamma - 1.));
			V.u  = 2. / (gamma + 1.) * (cL + (gamma - 1.) / 2. * uL + xi);
			V.p  = pL * std::pow(
						2. / (gamma + 1.) + (gamma - 1.)
							/ (gamma + 1.) / cL * (uL - xi),
						2. * gamma / (gamma - 1.));
		} else if (xi <= xiTailR) {
			V.ro = 0.;
			V.u  = 0.;
			V.p  = 0.;
		} else if (xi < xiHeadR) {
			V.ro = roR * std::pow(
						2. / (gamma + 1.) - (gamma - 1.)
							/ (gamma + 1.) / cR * (uR - xi),
						2. / (gamma - 1.));
			V.u  = 2. / (gamma + 1.) * (-cR + (gamma - 1.) / 2. * uR + xi);
			V.p  = pR * std::pow(
						2. / (gamma + 1.) - (gamma - 1.)
							/ (gamma + 1.) / cR * (uR - xi),
						2. * gamma / (gamma - 1.));
		} else {
			V.ro = roR;
			V.u  = uR;
			V.p  = pR;
		}
		return V;
	}
	T cLLocal = getc(res.roL, res.p),
			cRLocal = getc(res.roR, res.p);
	// Если не вакуум. Пусть точка слева от контактного разрыва
	// (xiContact = res.u)
	if (xi < res.u) {
		if (res.type == RPWaveConfig::swsw
				|| res.type == RPWaveConfig::swrw) {
			xiFront = uL - cL * std::sqrt(
						(gamma + 1.) / 2. / gamma * res.p / pL
						+ (gamma - 1.) / 2. / gamma);
			if (xi < xiFront) {
				V.ro = roL;
				V.u  = uL;
				V.p  = pL;
			} else {
				V.ro = res.roL;
				V.u = res.u;
				V.p = res.p;
			}
		} else if (res.type == RPWaveConfig::rwsw
				   || res.type == RPWaveConfig::rwrw) {
			xiHead = uL-cL;
			xiTail = res.u-cLLocal;
			if (xi <= xiHead) {
				V.ro = roL;
				V.u  = uL;
				V.p  = pL;
			} else if (xi >= xiTail) {
				V.ro = res.roL;
				V.u  = res.u;
				V.p  = res.p;
			} else {
				V.ro = roL * std::pow(
							2. / (gamma + 1.) + (gamma - 1.)
								/ (gamma + 1.) / cL * (uL - xi),
							2. / (gamma - 1.));
				V.u  = 2. / (gamma + 1)
						* (cL + (gamma - 1.) / 2. * uL + xi);
				V.p  = pL * std::pow(
							2. / (gamma + 1.) + (gamma - 1.)
								/ (gamma + 1.) / cL * (uL - xi),
							2. * gamma / (gamma - 1.));
			}
		}
	//Пусть точка справа от контактного разрыва (xiContact = res.v)
	} else {
		if (res.type ==  RPWaveConfig::rwsw
				|| res.type ==  RPWaveConfig::swsw) {
			xiFront = uR + cR * std::sqrt(
						(gamma + 1.) / 2. / gamma * res.p / pR
							+ (gamma - 1.) / 2. / gamma);
			if (xi > xiFront) {
				V.ro = roR;
				V.u  = uR;
				V.p  = pR;
			} else {
				V.ro = res.roR;
				V.u  = res.u;
				V.p  = res.p;
			}
		} else if (res.type == RPWaveConfig::rwrw
				   || res.type == RPWaveConfig::swrw) {
			xiHead = uR + cR;
			xiTail = res.u + cRLocal;
			if (xi >= xiHead) {
				V.ro = roR;
				V.u  = uR;
				V.p  = pR;
			} else if (xi <= xiTail) {
				V.ro = res.roR;
				V.u  = res.u;
				V.p  = res.p;
			} else {
				V.ro = roR * std::pow(
							2. / (gamma + 1.) - (gamma - 1.)
								/ (gamma + 1.) / cR * (uR - xi),
							2. / (gamma - 1.));
				V.u  = 2. / (gamma + 1)
						* (-cR + (gamma - 1.) / 2. * uR + xi);
				V.p  = pR * std::pow(
							2. / (gamma + 1.) - (gamma - 1.)
								/ (gamma + 1.) / cR * (uR - xi),
							2. * gamma / (gamma - 1.));
			}
		}
	}
	return V;
}


// Godunov-exact solver flux
template <ArithmeticWith<numeric_val> T>
Vector4<T> calcExactFlux(T roL, T rouL, T roEL,
						 T roR, T rouR, T roER, T gamma = 1.4) {
	T uL = 0., uR = 0., pL = 0., pR = 0.;
	if (roL != 0.) {
		uL = rouL / roL;
		pL = getP(roL, roEL/roL - .5*uL*uL, gamma);
	}

	if (roR != 0.) {
		uR = rouR / roR;
		pR = getP(roR, roER/roR - .5*uR*uR, gamma);
	}

	C1DVectorPrimitive<T> res = calcSolution<T>(
				roL, uL, pL, roR, uR, pR, 0., .01);
	T E = gete(res.ro, res.p, gamma) + .5*res.u*res.u;
	Vector4<T> FGodunov = Vector4(
				res.ro*res.u,
				res.ro*res.u*res.u + res.p,
				res.u*(res.ro*E+res.p),
				0.);

	return FGodunov;
}

#endif // EXACTSOLVER_H
