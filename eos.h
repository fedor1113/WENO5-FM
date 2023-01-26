#ifndef EOS_H
#define EOS_H

#include <cmath>

#include "Eigen/Dense"

#include "_vector4.h"
#include "arithmeticwith.h"

constexpr const numeric_val DEFAULT_GAMMA = 1.4;
// constexpr const numeric_val gamma = 1.4;

template <ArithmeticWith<numeric_val> T>
T get_enthalpy(T rho, T rho_v, T rho_E, T gamma = DEFAULT_GAMMA) {
	const T E = rho_E / rho;
	const T u = rho_v / rho;
	const T h = E + rho * (gamma - 1.) * (E - 0.5 * u * u);
//	const T e = E - .5 * u * u;
//	const T p = FEOSMieGruneisenAl<T>::getp(rho, e);
//	const T h = E + p;

	return h;
}


template <ArithmeticWith<numeric_val> T>
T gete(T rho, T p, T gamma = DEFAULT_GAMMA) {
	// return FEOSMieGruneisenAl<T>::gete(rho, p);

	if (rho != 0.)
		return p / (gamma - 1.) / rho;

	return 0.;
}


template <ArithmeticWith<numeric_val> T>
T getP(T rho, T e, T gamma = DEFAULT_GAMMA) {
	// return FEOSMieGruneisenAl<T>::getp(rho, e);

	return (gamma - 1.) * rho * e;
}


//template <ArithmeticWith<numeric_val> T>
//T eFromConservative(T rho, T j, T rhoE) {
//	return (rhoE - 0.5 * j*j / rho) / rho;
//}


template <ArithmeticWith<numeric_val> T>
Vector4<T> calcPhysicalFlux(T rho, T u, T p, T last,
							T gamma = DEFAULT_GAMMA) {
	/* Calculate for a Vector4<T> of conserved variables for
	 * the 1D Euler equations its corresponding flux.
	 */

	if (rho == 0) return Vector4<T>::ZERO;

	const T e = gete(rho, p, gamma);
	// const T e = FEOSMieGruneisenAl<T>::gete(rho, p);
	return Vector4<T>(rho * u, p + rho*u*u,
				u*(p + rho*(e + 0.5*u*u)),
				u*last);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> primitiveToConservative(Vector4<T> u/*, T gamma = DEFAULT_GAMMA*/) {
	// Conservative variables
	T gamma = DEFAULT_GAMMA;
	const T e = gete(u[0], u[2], gamma);
	// const T rho_E = u[2] / (gamma - 1.) + 0.5 * u[0] * u[1] * u[1];
	const T rho_E = u[0] * (e + 0.5 * u[1] * u[1]);

	return Vector4<T>(u[0], u[0] * u[1],
		rho_E, /*u[3]*/e);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> conservativeToPrimitive(Vector4<T> q, T gamma = DEFAULT_GAMMA) {
	// Primitive variables
	const T rho = q[0];
	const T u = q[1] / rho;
	const T E = q[2] / rho;
	// const T e = E - .5 * u * u;
	// const T p = FEOSMieGruneisenAl<T>::getp(rho, e);
	T p = (gamma - 1.) * rho * (E - 0.5*u*u);
	T e = gete(rho, p, gamma);

	return Vector4<T>(rho, u, p, e);
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> calcPhysicalFluxFromConservativeVec(
		Vector4<T> u,
		T gamma = DEFAULT_GAMMA) {
	/* Calculate for a Vector4<T> of conserved variables
	 * for the 1D Euler equations (rho, j=rho*v, rhoE=rho(e+v^2/2), smth)
	 * its corresponding flux.
	 */

//	return calcPhysicalFlux(u[0],
//			u[1] / u[0],
//			getP(u[0], eFromConservative(u[0], u[1], u[2])));
	Vector4<T> prim = conservativeToPrimitive(u, gamma);

	return calcPhysicalFlux(prim[0], prim[1], prim[2], prim[3],
		gamma);
}


template <ArithmeticWith<numeric_val> T>
T calcSquareSoundSpeed(T rho, T rho_v, T rho_E, T gamma = DEFAULT_GAMMA) {
	/* Compute the square of sound speed. */

	if (rho == static_cast<T>(0.)) return static_cast<T>(0.);

	// const T e = (rho_E - .5 * rho_v * rho_v / rho) / rho;
	// const T p = FEOSMieGruneisenAl<T>::getp(rho, e);
	const T p = (gamma - 1.) * (rho_E - rho_v * rho_v * 0.5 / rho);

	return std::abs(gamma * p / rho);
}


template <ArithmeticWith<numeric_val> T>
T calcMaxWaveSpeedDAtPoint(Vector4<T> u_vec_pt, T gamma = DEFAULT_GAMMA) {
	/* Calculate |df/du| for 1D Euler eq'ns. */

	return (std::sqrt(calcSquareSoundSpeed(
				u_vec_pt[0],
				u_vec_pt[1],
				u_vec_pt[2], gamma))
					+ std::abs(u_vec_pt[1] / u_vec_pt[0]));
}


template <ArithmeticWith<numeric_val> T>
Eigen::Matrix<T, 3, 3> EigenLeft1DEulerEigenMatrix(
		Vector4<T> vec, T gamma = DEFAULT_GAMMA) {
	assert(gamma > 1.);

	T gamma_m = gamma - 1.;

	T c_s_square = calcSquareSoundSpeed<T>(vec[0], vec[1], vec[2], gamma);
	T c_s = std::abs(std::sqrt(c_s_square));

	T u = vec[1];
	if (vec[0] != static_cast<T>(0.))
		u /= vec[0];

	T phi_square = 0.5 * gamma_m * u * u;
	T uc = u * c_s;
//	T h = get_enthalpy<T>(vec[0], vec[1], vec[2], gamma);
	T H = 0.5 * u * u + c_s_square / gamma_m;

//	Eigen::Matrix<T, 3, 3> l_mat {
//		{1. - phi_square / c_s_square,
//					gamma_m * u / c_s_square, -gamma_m / c_s_square},
//		{phi_square - uc, +c_s - gamma_m * u, gamma_m              },
//		{phi_square + uc, -c_s - gamma_m * u, gamma_m              }
//	};
	Eigen::Matrix<T, 3, 3> l_mat {
		{1.,      1.,          1.     },
		{u - c_s, u + 0.,      u + c_s},
		{H - uc,  0.5 * u * u, H + uc }
	};

	return l_mat;
}


template <ArithmeticWith<numeric_val> T>
Eigen::Matrix<T, 3, 3> EigenRight1DEulerEigenMatrix(
		Vector4<T> vec, T gamma = DEFAULT_GAMMA) {
	assert(gamma > 1.);

	T gamma_m = gamma - 1.;

	T c_s_square = calcSquareSoundSpeed<T>(vec[0], vec[1], vec[2], gamma);
	T beta = 1.;
	if (c_s_square != 0.)
		beta /= (2. * c_s_square);
	else
		beta = 0.;
	T c_s = std::sqrt(c_s_square);

	T u = vec[1];
	if (vec[0] != static_cast<T>(0.))
		u /= vec[0];

	T phi_square = 0.5 * gamma_m * u * u;
	T uc = u * c_s;
//	T h = get_enthalpy<T>(vec[0], vec[1], vec[2], gamma);
	T H = 0.5 * u * u + c_s_square / gamma_m;

//	Eigen::Matrix<T, 3, 3> r_mat {
//		{1.,                   beta,             beta            },
//		{u,                    beta * (u + c_s), beta * (u - c_s)},
//		{phi_square / gamma_m, beta * (H + uc),  beta * (H - uc) },
//	};

	Eigen::Matrix<T, 3, 3> r_mat {
		{H + c_s * (u - c_s) / gamma_m,       -(u + c_s / gamma_m), 1.},
		{-2. * H + 4. * c_s_square / gamma_m, 2. * u,              -2.},
		{H - c_s * (u + c_s) / gamma_m,       -u + c_s / gamma_m,   1.},
	};
	r_mat *= gamma_m * beta;

//	Eigen::Matrix<T, 3, 3> r_mat {
//		{2. * phi_square + c_s * u,       -(gamma_m * u + c_s), gamma_m},
//		{2. * c_s_square - 4. * phi_square, 2. * gamma_m * u, -2. * gamma_m},
//		{2. * phi_square - c_s * u,       -gamma_m * u + c_s,   gamma_m},
//	};
//	r_mat *= beta;

//	Eigen::Matrix<T, 3, 3> r_mat {
//		{gamma_m * H + c_s * (u - c_s),       -(gamma_m * u + c_s), 1. * gamma_m},
//		{-2. * gamma_m * H + 4. * c_s_square, 2. * gamma_m * u,              -2. * gamma_m},
//		{gamma_m * H - c_s * (u + c_s),       -gamma_m * u + c_s,   1. * gamma_m},
//	};

//	r_mat *= beta;

	return r_mat;
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> projectOntoCharacteristics(
		Vector4<T> conservative_variables, Vector4<T> vec,
		T gamma = DEFAULT_GAMMA) {
	if (vec[0] == static_cast<T>(0.))
			return Vector4<T>::ZERO;

	return Vector4<T>(EigenRight1DEulerEigenMatrix<T>(
				conservative_variables, gamma)
			* Eigen::Matrix<T, 3, 1>{vec[0], vec[1], vec[2]});
//	return Vector4<T>((Eigen::Matrix<T, 1, 3>{vec[0], vec[1], vec[2]}
//			* EigenRight1DEulerEigenMatrix(
//				conservative_variables, gamma)).transpose());
}


template <ArithmeticWith<numeric_val> T>
Vector4<T> projectCharacteristicVariablesBackOntoConserved(
		Vector4<T> conservative_variables, Vector4<T> vec,
		T gamma = DEFAULT_GAMMA) {
	if (vec[0] == static_cast<T>(0.))
			return Vector4<T>::ZERO;

	return Vector4<T>(EigenLeft1DEulerEigenMatrix<T>(
				conservative_variables, gamma)
			* Eigen::Matrix<T, 3, 1>{vec[0], vec[1], vec[2]});
//	return Vector4<T>((Eigen::Matrix<T, 1, 3>{vec[0], vec[1], vec[2]}
//			* EigenLeft1DEulerEigenMatrix(
//				conservative_variables, gamma)).transpose());
}

#endif // EOS_H
