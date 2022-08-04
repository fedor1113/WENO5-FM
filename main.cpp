#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#include <valarray>

#include "WENO5-ZM.h"

//class CTestToro {
//public:
//	double gamma;
//	double roL, vL, pL, eL, EL,
//		   roR, vR, pR, eR, ER;
//	double q0;   // Начальное значение координаты разрыва
//	double tMax; // Время, до которого ведется счет
//	CTestToro(double _gamma,
//			  double _roL, double _uL, double _vL, double _wL, double _pL,
//			  double _roR, double _uR, double _vR, double _wR, double _pR, double _q0, double _tMax) :
//			  gamma(_gamma), roL(_roL), vL(_vL), pL(_pL), roR(_roR), vR(_vR), pR(_pR), q0(_q0), tMax(_tMax) {
//				  eL = pL/(gamma-1.)/roL;
//				  EL = eL + (vL*vL)/2.;
//				  eR = pR/(gamma-1.)/roR;
//				  ER = eR + (vR*vR)/2.;
//			  }
//};


int main() {
	// std::ios::sync_with_stdio(false);
	// std::cin.tie(nullptr);

	// std::cout << "Hello World!" << "\n";
	const double Runiv = 8.314; // Universal gas constant, [J/K-mol]
	double mu1 = 0.0289647; double mu2 = 0.14606; // Mollecular weights of air and sf6 [kg/mol]
	double R1 = Runiv/mu1; double R2 = Runiv/mu2; // Specific gas constants
	double cp1 = 1005.9; double cp2 = 627.83; // Cp values for air and sf6
	double cv1 = 717.09; double cv2 = 566.95; // Cv values for air and sf6

	size_t Nx = 200;
	size_t N_full = Nx + 2*3;
	double L = 1.; // [L]
	double dx = L / Nx; // [L]
	// Nx = Nx + 6; // Add in ghost cells
	std::valarray<double> x(0., N_full);

	for (size_t k = 1; k < N_full; ++ k)
		x[k] = x[k-1] + dx;

	double cfl = 0.55;
	double t = 0;
	double tfinal = 0.2; // [T]

	std::valarray<Vec4> u_init(Vec4::ZERO, N_full);
	std::valarray<Vec4> flux(Vec4::ZERO, N_full);
	std::valarray<Vec4> Y2(N_full);
	std::valarray<Vec4> Y3(N_full);
	double Tatm = 293.0; // [K], approx 70 F
	double Patm = 101300.0; // [Pa], Atmospheric Pressure
	double Rhoatm = Patm / (Tatm * R1); // Density of the first gas at STP
	size_t k = 0;

	double P, c1, c2, R, G, rho;
	for (k = 0; k < size_t(N_full/2) + 1; ++ k) {
		P = 8. * Patm;
		// P = 1.;
		c1 = 1.; c2 = 1. - c1;
		R = c1*R1 + c2*R2;
		G = (c1*cp1+c2*cp2) / (c1*cv1+c2*cv2);
		// G = 1.4;
		rho = P / (R * Tatm);
		// rho = 1.;
		u_init[k][0] = rho / Rhoatm; // rho
		// u_init[k][0] = rho; // rho
		u_init[k][1] = 0.; // j = rho*u
		u_init[k][2] = (P/Patm) / (G-1.0); // rho*E
		// u_init[k][2] = P / (G-1.0); // rho*E
		u_init[k][3] = c1 * u_init[k][0];
	}

	for (k = size_t(N_full/2) + 1; k < N_full; ++ k) {
		P = 1. * Patm;
		// P = .1;
		c1 = 1.; c2 = 1. - c1;
		R = c1*R1 + c2*R2;
		G = (c1*cp1+c2*cp2) / (c1*cv1+c2*cv2);
		// G = 1.4;
		rho = P / (R * Tatm);
		// rho = 1.;
		// rho = 0.125;
		u_init[k][0] = rho / Rhoatm; // rho
		// u_init[k][0] = rho; // rho
		u_init[k][1] = 0.; // j = rho*u
		u_init[k][2] = (P/Patm) / (G-1.0); // rho*E
		// u_init[k][2] = P / (G-1.0); // rho*E
		u_init[k][3] = c1 * u_init[k][0];
	}

	std::valarray<double> s1 = std::valarray(0., N_full);
	for (k = 0; k < N_full; ++ k)
		s1[k] = u_init[k][3] / u_init[k][0]; // s2 = W[4,:]/W[0,:]
	std::valarray<double> s2 = 1 - s1;

	// size_t Nt = std::ceil(tfinal / dt);
	// std::valarray<int> Ts = std::valarray(0, Nt+2);
	// for (size_t k = 0; k <= Nt+1; ++ k) Ts[k] = k;

	// for (auto n : Ts) {
	integrate(u_init, flux, Y2, Y3, t, dx, Nx, tfinal, cfl);

	std::ofstream outfile;

	outfile.open("res.dat");

	k = 0;
	if (outfile.is_open())
		for (auto u : u_init) {
			Vec4 q = conservativeToPrimitive(u);
			// std::cout << u << "\n";

//			std::cout << x[k]
//					  << " " << q[0]
//					  << " " << q[1]
//					  << " " << q[2]
//					  << " " << q[3] << "\n";
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
