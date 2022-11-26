#ifndef EOS_MG_H
#define EOS_MG_H

#include <cmath>

#include "Eigen/Dense"

#include "_vector4.h"
#include "arithmeticwith.h"


constexpr const numeric_val rho0 = 2750.;
constexpr const numeric_val p0 = 560.964e9;
constexpr const numeric_val a = 1.12657;
constexpr const numeric_val b =0.975511;
constexpr const numeric_val ph = 15.e9;


template <ArithmeticWith<numeric_val> T>
T get_enthalpy(T rho, T rho_v, T rho_E) {
	const T E = rho_E / rho;
	const T u = rho_v / rho;
	const T e = E - .5 * u * u;
	const T p = getP(rho, e);
	const T h = E + p / rho;

	return h;
}


template <ArithmeticWith<numeric_val> T>
/*constexpr*/ T G([[maybe_unused]] T x) {
	return static_cast<T>(2.);
}


template <ArithmeticWith<numeric_val> T>
T pCold(T rho) {
	T x = rho / rho0;
	T pc = 0.;
	if(x >= 1.)
		pc = p0 * x * (std::pow(x, a) - std::pow(x, b));
	else {
		T n = p0 * (a - b) / ph;
		pc = ph * (std::pow(x, n) - 1.);
	}
	return pc;
}


template <ArithmeticWith<numeric_val> T>
T eCold(T rho) {  // [J/kg]
	T x = rho / rho0;
	T ec = 0.;
	if (x >= 1.)
		ec = 2.03987e8 * (std::pow(x, a) / a - std::pow(x, b) / b);
	else {
		T n = p0 * (a - b) / ph;
		ec = 2.03987e8 * (a - b) * (1. / n * (std::pow(x, n - 1.) / (n - 1.) + 1. / x) - 1. / (n - 1.) - 1. / a / b);
	}
	return ec;
}


template <ArithmeticWith<numeric_val> T>
T gete(T rho, T p) {
	if (rho != 0.) {
		T x = rho / rho0;
		T _ec = eCold(rho);
		T _pc = pCold(rho);

		return _ec + 1. / rho / G(x) * (p - _pc);
	}

	return 0.;
}


template <ArithmeticWith<numeric_val> T>
T getP(T rho, T e) {
	T x = rho / rho0;

	return pCold(rho) + G(x) * rho * (e-eCold(rho));
}


//template <ArithmeticWith<numeric_val> T>
//T eFromConservative(T rho, T j, T rhoE) {
//	return (rhoE - 0.5 * j*j / rho) / rho;
//}


template <ArithmeticWith<numeric_val> T>
T pColdPrime(T rho) {
	T x = rho / rho0;
	T pcx = 0.;
	if (x >= 1.)
		pcx = p0 * ((a + 1.) * std::pow(x, a) - (b + 1.) * std::pow(x, b));
	else {
		T n = p0 * (a - b) / ph;
		pcx = p0 * (a - b) * std::pow(x, n - 1.);
	}

	return pcx;
}


template <ArithmeticWith<numeric_val> T>
T getdpdrho(T rho, T e) {
	const T x = rho / rho0;
	/*constexpr*/ const T _G = G(x);
	const T _pc = pCold(rho);
	const T _pcx = pColdPrime(rho);
	const T _ec = eCold(rho);

	return 1. / rho0 * (_pcx + _G * rho0 * (e - _ec) - _G / x * _pc);
}


template <ArithmeticWith<numeric_val> T>
T getdpde(T rho, T e) {
	/*constexpr*/ const T _G = G(rho / rho0);

	return rho * _G;
}


template <ArithmeticWith<numeric_val> T>
T getc(T rho, T p) {
	if (rho == static_cast<T>(0.))  return static_cast<T>(0.);

	T _e = gete(rho, p);
	T _prho = getdpdrho(rho, _e);
	T _pe = getdpde(rho, _e);
	T c2 = p * _pe / rho / rho + _prho;

	return std::sqrt(std::abs(c2));
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> calcPhysicalFlux(T rho, T u, T p, T last) {
	/* Calculate for a Vector4<T> of conserved variables for
	 * the 1D Euler equations its corresponding flux.
	 */

	if (rho == static_cast<T>(0.)) return Vector4<T>::ZERO;

	const T e = gete(rho, p);
	// const T e = FEOSMieGruneisenAl<T>::gete(rho, p);
	return Vector4<T>(rho * u, p + rho*u*u,
				u*(p + rho*(e + 0.5*u*u)),
				u*last);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> primitiveToConservative(Vector4<T> u) {
	// Conservative variables
	const T e = gete(u[0], u[2]);
	// const T rho_E = u[2] / (gamma - 1.) + 0.5 * u[0] * u[1] * u[1];
	const T rho_E = u[0] * (e + 0.5 * u[1] * u[1]);

	return Vector4<T>(u[0], u[0] * u[1],
		rho_E, /*u[3]*/e);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> conservativeToPrimitive(Vector4<T> q) {
	// Primitive variables
	const T rho = q[0];
	const T u = q[1] / rho;
	const T E = q[2] / rho;
	const T e = E - .5 * u * u;
	const T p = getP(rho, e);

	return Vector4<T>(rho, u, p, e);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> calcPhysicalFluxFromConservativeVec(
		Vector4<T> u) {
	/* Calculate for a Vector4<T> of conserved variables
	 * for the 1D Euler equations (rho, j=rho*v, rhoE=rho(e+v^2/2), smth)
	 * its corresponding flux.
	 */

//	return calcPhysicalFlux(u[0],
//			u[1] / u[0],
//			getP(u[0], eFromConservative(u[0], u[1], u[2])));
	Vector4<T> prim = conservativeToPrimitive(u);

	return calcPhysicalFlux(prim[0], prim[1], prim[2], prim[3]);
}


template <ArithmeticWith<numeric_val> T>
T calcSquareSoundSpeed(T rho, T rho_v, T rho_E) {
	/* Compute the square of sound speed. */

	const T e = (rho_E - .5 * rho_v * rho_v / rho) / rho;
	const T c = getc(rho, getP(rho, e));

	return c * c;
}


template <ArithmeticWith<numeric_val> T>
T calcMaxWaveSpeedDAtPoint(Vector4<T> u_vec_pt) {
	/* Calculate |df/du| for 1D Euler eq'ns. */

	return (std::sqrt(calcSquareSoundSpeed(
				u_vec_pt[0],
				u_vec_pt[1],
				u_vec_pt[2]))
					+ std::abs(u_vec_pt[1] / u_vec_pt[0]));
}


template <ArithmeticWith<numeric_val> T>
Eigen::Matrix<T, 3, 3> EigenLeft1DEulerEigenMatrix(
		Vector4<T> vec) {
	// assert(gamma > 1.);

	// T gamma_m = gamma - 1.;

	T u = vec[1];
	if (vec[0] != static_cast<T>(0.))
		u /= vec[0];

//	double phi_square = 0.5 * gamma_m * u * u;
	T e = vec[2] / vec[0] - static_cast<T>(.5) * u * u;
	T p = getP(vec[0],e);
	T c = getc(vec[0], p);

	T c_s_square = c * c;
	T c_s = std::abs(std::sqrt(c_s_square));

	T uc = u * c_s;
	T h = (vec[3] + p) / vec[0];  // total specific enthalpy
//	if (vec[0] != 0.)
//		h = (vec[3] + p) / vec[0];

//	Eigen::Matrix<T, 3, 3> l_mat {
//		{1. - phi_square / c_s_square,
//					gamma_m * u / c_s_square, -gamma_m / c_s_square},
//		{phi_square - uc, +c_s - gamma_m * u, gamma_m              },
//		{phi_square + uc, -c_s - gamma_m * u, gamma_m              }
//	};
//	Eigen::Matrix<double, 3, 3> l_mat {
//		{1.,      1.,          1.     },
//		{u - c_s, u + 0.,      u + c_s},
//		{H - uc,  0.5 * u * u, H + uc }
//	};

	T b = 0.;
	if (vec[0] != static_cast<T>(0.))
		b = getdpde(vec[0], e) / vec[0];

	Eigen::Matrix<T, 3, 3> l_mat {
		{1.,                           1.,        1.},
		{u - c_s,                  u + 0.,   u + c_s},
		{h -  uc,      h - c_s_square / b,   h +  uc},
	};

	return l_mat;
}


template <ArithmeticWith<numeric_val> T>
Eigen::Matrix<T, 3, 3> EigenRight1DEulerEigenMatrix(
		Vector4<T> vec) {
	// assert(gamma > 1.);

	// T gamma_m = gamma - 1.;

	T u = vec[1];
	if (vec[0] != static_cast<T>(0.))
		u /= vec[0];

	//	double phi_square = 0.5 * gamma_m * u * u;
	T e = vec[2] / vec[0] - static_cast<T>(.5) * u * u;
	T p = getP(vec[0],e);
	T c = getc(vec[0], p);

	T c_s_square = c * c;
	T beta = 1.;
	if (c_s_square != 0.)
		beta /= (2. * c_s_square);
	//	else
	//		beta = 0.;
	T c_s = std::sqrt(c_s_square);

	T uc = u * c_s;
	T h = (vec[3] + p) / vec[0];  // total specific enthalpy

//	Eigen::Matrix<T, 3, 3> r_mat {
//		{1.,                   beta,             beta            },
//		{u,                    beta * (u + c_s), beta * (u - c_s)},
//		{phi_square / gamma_m, beta * (H + uc),  beta * (H - uc) },
//	};

//	Eigen::Matrix<T, 3, 3> r_mat {
//		{H + c_s * (u - c_s) / gamma_m,       -(u + c_s / gamma_m), 1.},
//		{-2. * H + 4. * c_s_square / gamma_m, 2. * u,              -2.},
//		{H - c_s * (u + c_s) / gamma_m,       -u + c_s / gamma_m,   1.},
//	};
//	r_mat *= gamma_m * beta;

//	Eigen::Matrix<T, 3, 3> r_mat {
//		{gamma_m * H + c_s * (u - c_s),       -(gamma_m * u + c_s), 1. * gamma_m},
//		{-2. * gamma_m * H + 4. * c_s_square, 2. * gamma_m * u,              -2. * gamma_m},
//		{gamma_m * H - c_s * (u + c_s),       -gamma_m * u + c_s,   1. * gamma_m},
//	};

//	r_mat *= beta;
	T dpde = getdpde(vec[0], e);
	T b = 0.;
	if (vec[0] != static_cast<T>(0.))
		b = dpde / vec[0];

	// harder to compute for MG
	// double dpdrho = eos.getdpdrho(vec[0], e);
	T dpdrho = c_s_square - p * b / vec[0];

	T theta = u * u
			- vec[3] / vec[0]
			+ vec[0] * dpdrho / dpde;

	Eigen::Matrix<T, 3, 3> r_mat {
		{theta  +  uc / b, -(u + c_s / b),    1.},
		{2. * (h - u * u),         2. * u,   -2.},
		{theta  -  uc / b,   -u + c_s / b,    1.},
	};
	r_mat *= b * beta;

	return r_mat;
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> projectOntoCharacteristics(
		Vector4<T> conservative_variables, Vector4<T> vec) {
	return Vector4<T>(EigenRight1DEulerEigenMatrix<T>(
				conservative_variables)
			* Eigen::Matrix<T, 3, 1>{vec[0], vec[1], vec[2]});
//	return Vector4<T>((Eigen::Matrix<T, 1, 3>{vec[0], vec[1], vec[2]}
//			* EigenRight1DEulerEigenMatrix(
//				conservative_variables, gamma)).transpose());
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> projectCharacteristicVariablesBackOntoConserved(
		Vector4<T> conservative_variables, Vector4<T> vec) {
	return Vector4<T>(EigenLeft1DEulerEigenMatrix<T>(
				conservative_variables)
			* Eigen::Matrix<T, 3, 1>{vec[0], vec[1], vec[2]});
//	return Vector4<T>((Eigen::Matrix<T, 1, 3>{vec[0], vec[1], vec[2]}
//			* EigenLeft1DEulerEigenMatrix(
//				conservative_variables, gamma)).transpose());
}

#endif // EOS_MG_H
