#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimize("unroll-loops")

#include <algorithm>
#include <cstdio>
#include <fstream>
// #include <initializer_list>
#include <iostream>
#include <numeric>
#include <valarray>
#include <vector>

// #include "vectorfieldcomponentview.h"
// #include "euler1d.h"
// #include "kfr1d.h"
#include "inviscidburgers1d.h"

// template class Vector4<numeric_val>;
// using Vec4 = Vector4<numeric_val>;


int main(int argc, char **argv) {
	std::ios::sync_with_stdio(false);
	// std::cin.tie(nullptr);

	std::size_t k = 0;

	std::size_t Nx = 201;
	const std::size_t N_ghost_points = 3;
	std::size_t N_full = Nx + 2*N_ghost_points;
	// Nx = Nx + 6;  // Add in ghost cells
	numeric_val L = 1.;  // [L]

	numeric_val cfl = 0.4;
	numeric_val t = 0.;
	numeric_val tfinal = 0.0;  // [T]

	// std::valarray<Vec4> u_init(Vec4::ZERO, N_full);
	// std::valarray<Vec4> flux(Vec4::ZERO, N_full);
	// std::valarray<Vec4> Y2(Vec4::ZERO, N_full);
	// std::valarray<Vec4> Y3(Vec4::ZERO, N_full);

	// numeric_val gamma = 1.4;

//	std::valarray<numeric_val> s1(static_cast<numeric_val>(0.), N_full);
//	for (k = 0; k < N_full; ++ k)
//		s1[k] = u_init[k][3] / u_init[k][0];  // s2 = W[4,:]/W[0,:]
//	std::valarray<numeric_val> s2 = 1 - s1;

	// std::size_t Nt = std::ceil(tfinal / dt);
	// std::valarray<int> Ts = std::valarray(0, Nt+2);
	// for (std::size_t k = 0; k <= Nt+1; ++ k) Ts[k] = k;

	// for (auto n : Ts) {

	// std::valarray<Vec4> u_res(Vec4::ZERO, N_full);
//	for (auto n : VectorFieldComponentView<std::valarray<Vec4>>(u_res, 3))
//		std::cout << n << "\n";


	// tfinal = 0.2;
	// solve1DRiemannProblemForEulerEq<numeric_val>(
	// 	u_res, x, gamma,
	// 	1., 0., 8.,
	// 	1., 0., 1., 0.5,
	// 	t, tfinal, 0., L,
	// 	primitiveToConservativeU<numeric_val>, Nx, cfl
	// );
	// Nx = 11; N_full = Nx + 2*N_ghost_points;
	// L = 1.;
	// cfl = 0.9;
	// tfinal = 0.5;
	// solve1DInviscidBurgersProblem<numeric_val>(
	// 	u_res, x, 0., 0.5, -1., 1.,  Nx, cfl
	// );

//	tfinal = 0.15;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		1., 0., 1.,
//		0.125, 0., 0.1, 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Sod's problem (expansion-contact-shock)

//	Nx = 200; N_full = Nx + 2*N_ghost_points;
//	L = 1.;
//	cfl = 0.55;
//	tfinal = 0.2;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		1., 0.75, 1.,
//		0.125, 0., 0.1, 0.3,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Modified Sod's problem (expansion-contact-shock). Toro-1

//	tfinal = 0.12;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		.445, .698, 3.528,
//		.5, 0., .571, 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Lax's problem (expansion-contact-shock)

//	tfinal = 2.5 * std::pow(10., -6);
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		1., 0., std::pow(10, 10),
//		.125, 0., .1, 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Strong shock tube problem (expansion-contact-shock)

//	tfinal = 0.09;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		3.857, 0.92, 10.3333,
//		1., 3.55, 1., 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Mach 3 shock test (expansion-contact-shock)

//	tfinal = 1.75 * std::pow(10., -4);
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		10., 2000., 500.,
//		20., 0., 500., 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // High Mach flow test (shock-contact-shock)

//	Nx = 200; N_full = Nx + 2 * N_ghost_points;
//	L = 1.;
//	cfl = 0.55;
//	tfinal = 0.15;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		1., -2., 0.4,
//		1., 2., 0.4, 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Two symmetric rarefaction waves (expansion-contact-expansion). Toro-2

//	Nx = 200; N_full = Nx + 2 * N_ghost_points;
//	L = 1.;
//	cfl = 0.55;
//	tfinal = 0.012;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		1., 0., 1000.,
//		1., 0., .01, 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Toro-3

//	tfinal = 0.02;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		3., 0., 1000.,
//		2., 0., .1, 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Left Expansion and right strong shock

//	Nx = 200; N_full = Nx + 2*N_ghost_points;
//	L = 1.;
//	cfl = 0.55;
//	tfinal = 0.035;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		5.99924, 19.5975, 460.894,
//		5.99242, -6.19633, 46.0950, 0.4,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Toro-4

//	Nx = 200; N_full = Nx + 2*N_ghost_points;
//	L = 1.;
//	cfl = 0.55;
//	tfinal = 0.012;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		1., -19.59745, 1000.,
//		1., -19.59745, .01, 0.8,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Toro-5

//	tfinal = 2;
//	solve1DRiemannProblemForEulerEq<numeric_val>(
//		u_res, x, gamma,
//		27./7., 4.*std::sqrt(35)/9., 31./3.,
//		1+0.2*std::sin(5*x), -5., 5., 0.5,
//		t, tfinal, 0., L,
//		primitiveToConservativeU<numeric_val>, Nx, cfl
//	);  // Shu-Osher test - not a Riemann problem, so no dice

	for (std::size_t j
//		 : {51, 101, 151, 201, 251, 501, 751, 801,
//			1001, 1101, 1251, 1501,
//			10001}) {
		 : {101}) {
		Nx = j; N_full = Nx + 2 * N_ghost_points;
		L = 10.;
		// L = 1.;
		cfl = 0.8;
		tfinal = 3000.;
		// tfinal = 1.5;


		std::valarray<numeric_val> u_res(0., N_full);
//		std::size_t guessed_length = static_cast<std::size_t>(
//					(tfinal - t) / (cfl * L / Nx / 10)
//					);
//		std::vector<numeric_val> u_s(guessed_length, 0.);
//		std::vector<numeric_val> times(u_s.size(), 0.);
		std::valarray<numeric_val> x(0., N_full);

//		solve1DDetonationProfileProblem<numeric_val>(
//					u_res, x,
//					4.5, 0.1, 0.05, 0.01,
//					u_s, times,
//					1e-40,
//					0., tfinal, -L, 0.,
//					Nx, cfl
//					);
//		solve1DInviscidBurgersProblem<numeric_val>(
//			u_res, x, 0., 0.5, -1., 1.,  Nx, cfl
//		);
//		solve1DInviscidBurgersProblem<numeric_val>(
//			u_res, x,
//			0., 0.5/std::numbers::pi_v<numeric_val>, -1., 1.,
//			Nx, cfl
//		);
		tfinal = 1.;
		L = 1.;
		solve1DInviscidBurgersProblem<numeric_val>(
			u_res, x,
			0., tfinal, -L, L,
			Nx, cfl, 1
		);  // Toro-1

		std::ofstream outfile;

		outfile.open("res_n_" + std::to_string(j) + ".dat");

		k = 0;
		if (outfile.is_open()) {
			// outfile << "TITLE=\"Riemann Problem 1D slice t="
			// 		<< tfinal << "\"" << "\n";
			// outfile << "VARIABLES=\"x\",\"rho\",\"u\",\"p\",\"e\"" << "\n";
			// outfile << "VARIABLES=\"x\",\"rho\",\"u\",\"p\"" << "\n";
			outfile << "TITLE=\"Detonation Problem 1D slice t="
					<< tfinal << "\"" << "\n";
			outfile << "VARIABLES=\"x\",\"u\"" << "\n";
			outfile << "ZONE T=\"Numerical\", I="
					<< N_full << ", F=POINT" << "\n";

			for (auto u : u_res) {
				// Vec4 q = conservativeToPrimitive(u, gamma);
				// std::cout << u << "\n";

				//			std::cout << x[k]
				//					  << " " << q[0]
				//					  << " " << q[1]
				//					  << " " << q[2]
				//					  << " " << q[3] << "\n";
				//			G = (u[3]/u[0]*cp1+(1-u[3]/u[0])*cp2)/(u[3]/u[0]*cv1+(1-u[3]/u[0])*cv2);
				//			outfile << x[k]
				//					<< " " << q[0]
				//					<< " " << u[1]/u[0]
				//					<< " " << (u[2]-u[1]*(u[1]/u[0])/2)*(G-1)
				//					<< " " << u[3]/u[0] << "\n";
				// outfile << x[k]
				// 		<< " " << q[0]
				// 		<< " " << q[1]
				// 		<< " " << q[2]
				// 		<< " " << q[3] << "\n";
				outfile << x[k] << " " << u << "\n";

				++ k;
			}
		}

		outfile.close();


//		outfile.open("u_res_n_" + std::to_string(j) + ".dat");

//		k = 0;
//		if (outfile.is_open()) {
//			// outfile << "TITLE=\"Riemann Problem 1D slice t="
//			// 		<< tfinal << "\"" << "\n";
//			// outfile << "VARIABLES=\"x\",\"rho\",\"u\",\"p\",\"e\"" << "\n";
//			// outfile << "VARIABLES=\"x\",\"rho\",\"u\",\"p\"" << "\n";
//			outfile << "TITLE=\"Detonation Problem u_s time evolution"
//					<< tfinal << "\"" << "\n";
//			outfile << "VARIABLES=\"x\",\"u\"" << "\n";
//			outfile << "ZONE T=\"Numerical\", I="
//					<< times.size() << ", F=POINT" << "\n";

//			for (std::size_t k = 0; k < times.size(); ++ k) {
//				outfile << times[k] << " " << u_s[k] << "\n";
//			}
//		}

//		outfile.close();


		std::cout << "Done!" << "\n";
	}
	return 0;
}
