#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization("unroll-loops")
#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#include <valarray>

// #include "vectorfieldcomponentview.h"
#include "WENO5-FM.h"

using numeric_val = double;
template class Vector4<numeric_val>;
using Vec4 = Vector4<numeric_val>;

//class CTestToro {
//public:
//	numeric_val gamma;
//	numeric_val roL, vL, pL, eL, EL,
//		   roR, vR, pR, eR, ER;
//	numeric_val q0;   // Начальное значение координаты разрыва
//	numeric_val tMax; // Время, до которого ведется счет
//	CTestToro(numeric_val _gamma,
//			  numeric_val _roL, numeric_val _uL, numeric_val _vL, numeric_val _wL, numeric_val _pL,
//			  numeric_val _roR, numeric_val _uR, numeric_val _vR, numeric_val _wR, numeric_val _pR, numeric_val _q0, numeric_val _tMax) :
//			  gamma(_gamma), roL(_roL), vL(_vL), pL(_pL), roR(_roR), vR(_vR), pR(_pR), q0(_q0), tMax(_tMax) {
//				  eL = pL/(gamma-1.)/roL;
//				  EL = eL + (vL*vL)/2.;
//				  eR = pR/(gamma-1.)/roR;
//				  ER = eR + (vR*vR)/2.;
//			  }
//};


int main() {
	 std::ios::sync_with_stdio(false);
	// std::cin.tie(nullptr);

	std::size_t k = 0;
	const numeric_val Runiv = 8.314;  // Universal gas constant, [J/K-mol]
	numeric_val mu1 = 0.0289647; numeric_val mu2 = 0.14606;  // Mollecular weights of air and sf6 [kg/mol]
	numeric_val R1 = Runiv/mu1; numeric_val R2 = Runiv/mu2;  // Specific gas constants
	numeric_val cp1 = 1005.9; numeric_val cp2 = 627.83;  // Cp values for air and sf6
	numeric_val cv1 = 717.09; numeric_val cv2 = 566.95;  // Cv values for air and sf6

	std::size_t Nx = 1000;
	const std::size_t N_ghost_points = 3;
	std::size_t N_full = Nx + 2*N_ghost_points;
	// Nx = Nx + 6;  // Add in ghost cells
	numeric_val L = 1.;  // [L]

	numeric_val cfl = 0.3;
	numeric_val t = 0;
	numeric_val tfinal = 0.15;  // [T]

	std::valarray<Vec4> u_init(Vec4::ZERO, N_full);
	std::valarray<Vec4> flux(Vec4::ZERO, N_full);
	std::valarray<Vec4> Y2(Vec4::ZERO, N_full);
	std::valarray<Vec4> Y3(Vec4::ZERO, N_full);
	numeric_val Tatm = 293.0;  // [K], approx 70 F
	numeric_val Patm = 101300.0;  // [Pa], Atmospheric Pressure
	numeric_val Rhoatm = Patm / (Tatm * R1);  // Density of the first gas at STP

	numeric_val c1, c2, R, G;

	c1 = 1.; c2 = 1. - c1;
	R = c1*R1 + c2*R2;
	G = (c1*cp1+c2*cp2) / (c1*cv1+c2*cv2);
//	std::valarray<numeric_val> s1(static_cast<numeric_val>(0.), N_full);
//	for (k = 0; k < N_full; ++ k)
//		s1[k] = u_init[k][3] / u_init[k][0];  // s2 = W[4,:]/W[0,:]
//	std::valarray<numeric_val> s2 = 1 - s1;

	// std::size_t Nt = std::ceil(tfinal / dt);
	// std::valarray<int> Ts = std::valarray(0, Nt+2);
	// for (std::size_t k = 0; k <= Nt+1; ++ k) Ts[k] = k;

	// for (auto n : Ts) {

	std::valarray<Vec4> u_res(Vec4::ZERO, N_full);
	u_res[0][3] = 1;
	std::valarray<numeric_val> x(0., N_full);
//	for (auto n : VectorFieldComponentView<std::valarray<Vec4>>(u_res, 3))
//		std::cout << n << "\n";

//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, G,
//		8., 0., 8.,
//		1., 0., 1., 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);

	solve1DRiemannProblemForEulerEq<numeric_val>(
		u_res, x, 1.4,
		1., 0., 1.,
		0.125, 0., 0.1, 0.5,
		t, tfinal, 0., L,
		primitiveToConservativeU<numeric_val>, Nx, cfl
	);

	std::ofstream outfile;

	outfile.open("res500.dat");

	k = 0;
	if (outfile.is_open())
		for (auto u : u_res) {
			Vec4 q = conservativeToPrimitive(u, 1.4);
			// std::cout << u << "\n";

//			std::cout << x[k]
//					  << " " << q[0]
//					  << " " << q[1]
//					  << " " << q[2]
//					  << " " << q[3] << "\n";
//			G = (u[3]/u[0]*cp1+(1-u[3]/u[0])*cp2)/(u[3]/u[0]*cv1+(1-u[3]/u[0])*cv2);
//			outfile << x[k]
//				   << " " << q[0]
//				   << " " << u[1]/u[0]
//				   << " " << (u[2]-u[1]*(u[1]/u[0])/2)*(G-1)
//				   << " " << u[3]/u[0] << "\n";
			outfile << x[k]
				   << " " << q[0]
				   << " " << q[1]
				   << " " << q[2]
				   << " " << q[3] << "\n";

			++ k;
		}

	outfile.close();

	std::cout << "Done!" << "\n";

	return 0;
}
