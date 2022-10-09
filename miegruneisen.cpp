#include "miegruneisen.h"

#include <cmath>


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::getG(T ro) {
//	T x = ro/ro0;
//	const T a0 = 2.95, a1 = 2.408, a2 = 12.151;
//	const T M = 18.;     // [g/mole]
//	const T CVLiq = 4150.;   // [J/kg/K]
//	const T CVGas = 1430.;	  // [J/kg/K]
//	const T R = 8.31;        // [J/mole/K]
//	// Здесь используем только жидкую фазу! Можно будет дифференцировать при усложнении
//	const T CV = CVLiq;
//	T G = R/CV/M*(a0 + (1.-a0)*std::exp(-std::pow(x/.5273, 1.7)) + a1*std::exp(-std::pow(x/1.0904, -3.5)) + a2*std::exp(-std::pow(x/1.3927, -5.)));
//	return G;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::getp0(T ro) {
//	const T A = .6726e9; // [Pa]
//	const T b = 11.55;
//	const T K = 1.15e9;  // [Pa]
//	const T beta = .3333;
//	const T xi = .85;
//	// Cold component, Born-Meyer potential
//	T x = ro/ro0;
//	T p0 =  x*(A*std::pow(x, -beta) * std::exp(b*(1.-std::pow(x, -beta))) - K*(std::pow(x, xi)));
//	return p0;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::gete0(T ro) {
//	const T A = .6726e9; // [Pa]
//	const T b = 11.55;
//	const T K = 1.15e9;  // [Pa]
//	const T beta = .3333;
//	const T xi = .85;
//	// Cold component, Born-Meyer potential
//	T x = ro/ro0;
//	T e0 = 1./ro0*(A/beta/b * std::exp(b*(1.-std::pow(x, -beta))) - K/xi*std::pow(x, xi));
//	return e0;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::getp(T ro, T e) {
//	T x = ro/ro0;
//	T p0 = getp0(ro);
//	T e0 = gete0(ro);
//	T G = getG(ro);
//	T p = p0 + ro*G*(e-e0);
//	return p;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::gete(T ro, T p) {
//	T x = ro/ro0;
//	T p0 = getp0(ro);
//	T e0 = gete0(ro);
//	T G = getG(ro);
//	T e = e0 + 1./ro/G*(p-p0);
//	return e;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::getc(T ro, T p) {
//	T x = ro/ro0;
//	T G = getG(ro);
//	T p0 = getp0(ro);
//	T p0Prime = getp0Prime(ro);
//	T GPrime = getGPrime(ro);
//	T c2 = 1/ro0*(p0Prime + (GPrime/G + 1./x + G/x)*(p-p0));
//	return std::sqrt(c2);
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::getGPrime(T ro) {
//	T x = ro/ro0;
//	const T a0 = 2.95, a1 = 2.408, a2 = 12.151;
//	const T M = 18.;     // [g/mole]
//	const T CVLiq = 4150.;   // [J/kg/K]
//	const T CVGas = 1430.;	  // [J/kg/K]
//	const T R = 8.31;        // [J/mole/K]
//	// Здесь используем только жидкую фазу! Можно будет дифференцировать при усложнении
//	const T CV = CVLiq;
//	/* T GPrime = R/CV/M*((1.-a0)*exp(-pow(x/.5273, 1.7))*pow(x/.5273, -2.7)*(-1.7/.5273) +
//								a1 *exp(-pow(x/1.0904, -3.5))*pow(x/1.0904, -4.5)*(3.5/1.0904) +
//								a2 *exp(-pow(x/1.3927, -5.0))*pow(x/1.3927, -6.0)*(5.0/1.3927)); */
//	T GPrime = .035412*std::exp(-5.23948*std::pow(x, -5.))*std::pow(x, -6)
//				  + .00126927*std::exp(-1.35379*std::pow(x, -3.5))*std::pow(x, -4.5)
//				  + .00109463*std::exp(-2.96826*std::pow(x, 1.7))*std::pow(x, .7);
//	return GPrime;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::getp0Prime(T ro) {
//	const T A = .7626e9; // [Pa]
//	const T b = 11.55;
//	const T K = 1.15e9;  // [Pa]
//	const T beta = .3333;
//	const T xi = .85;
//	T x = ro/ro0;
//	T p0Prime = /*A*pow(x, -beta+1.) * exp(b*(1.-pow(x, -beta))) - K*(pow(x, xi+1.)) +
//					 A*(-beta+1.)*pow(x, -beta)*exp(b*(1.-pow(x, -beta))) + A*pow(x, -beta+1.) * exp(b*(1.-pow(x, -beta))) *
//		x*(         A*pow(x, -beta+1.) * exp(b*(1.-pow(x, -beta))) - K*(pow(x, xi+1.)));
//		*/
//		/*A*pow(x, -beta)*exp(b*(1.-pow(x, -beta))) - K*pow(x, xi) +
//		A*beta*exp(b*(1.-pow(x, -beta)))*(-pow(x, -beta) + pow(x, -2.*beta)) - K*xi*pow(x, xi);*/

//		std::exp(-11.55/std::pow(x,.3333))*(2.68705e14/std::pow(x, .6666) + 4.65359e13/std::pow(x,.3333)) - 2.1275e9*std::pow(x,.85);
//	return p0Prime;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisen<T>::getKSPrime(T ro, T e) {
//	T p = getp(ro, e), c = getc(ro, p), KS = ro*c*c, G = getG(ro), x = ro/ro0;
//	// Partial derivatives
//	// G'(ro)
//	T a0 = 2.95, a1 = 2.408, a2 = 12.151;
//	T b0 = .5273, b1 = 1.0904, b2 = 1.3927;
//	T c0 = 1.7, c1 = -3.5, c2 = -5.;
//	T dGdx = - c0*(1.-a0)/b0*std::exp(-std::pow(x/b0, c0))*std::pow(x/b0, c0-1.) - a1*c1/b1*std::exp(-std::pow(x/b1, c1))*std::pow(x/b1, c1-1.) - a2*c2/b2*std::exp(-std::pow(x/b2, c2))*std::pow(x/b2, c2-1.);
//	T dGdro = ro0*dGdx;
//	// G"(ro)
//	T d2Gdx2 = - c0*(1.-a0)/b0/b0*std::exp(-std::pow(x/b0, c0))*(-c0*std::pow(x/b0, 2.*c0-2.)+(c0-1.)*std::pow(x/b0, c0-2.))
//						 - a1*c1/b1/b1*std::exp(-std::pow(x/b1, c1))*(-c1*std::pow(x/b1, 2.*c1-2.)+(c1-1.)*std::pow(x/b1, c1-2.))
//						 - a2*c2/b2/b2*std::exp(-std::pow(x/b2, c2))*(-c2*std::pow(x/b2, 2.*c2-2.)+(c2-1.)*std::pow(x/b1, c2-2.));
//	T d2Gdro2 = ro0*ro0*d2Gdx2;
//	// p'(ro)
//	T A = .6726e9, beta = .3333, b = 11.55, K = 1.15e9, xi = .85;
//	T dp0dx = A*std::pow(x, -beta)*std::exp(b*(1.-std::pow(x, -beta))) - K*std::pow(x, xi) + A*beta*std::exp(b*(1.-std::pow(x, -beta)))*(-std::pow(x, -beta) + std::pow(x, -2.*beta)) - K*xi*std::pow(x, xi);
//	T dp0dro = ro*dp0dx;
//	T e0 = 1000./ro0*(A/beta/b * std::exp(b*(1.-std::pow(x, -beta))) - K/xi*std::pow(x, xi));
//	T de0dx = 1./ro0*(A*std::exp(b*1.-std::pow(x, -beta))*std::pow(x, -beta-1.) - K*std::pow(x, xi-1));
//	T de0dro = ro0*de0dx;
//	T dpdro = dp0dro + (G+ro*dGdro)*(e-e0) - ro*G*de0dro;
//	// p'(e)
//	T dpde = ro*G;
//	T d2p0dx2 = A*std::exp(b*(1.-std::pow(x, -beta)))*(-beta*(1.-beta)*std::pow(x, -beta-1.) + (b*beta-(2.+b)*beta*beta)*std::pow(x, -2.*beta-1.) + b*beta*beta*std::pow(x, -3.*beta-1.)) - K*xi*(xi+1.)*std::pow(x, xi-1.);
//	T d2p0dro2 = ro0*ro0*d2p0dx2;
//	T d2e0dx2 = 1./ro0*(A*std::exp(b*(1-std::pow(x, -beta)))*(b*beta*std::pow(x, 2.*beta-2.) - (beta+1.)*std::pow(x, -beta-2.)) - K*std::pow(x, xi-2.));
//	T d2e0dro2 = ro0*ro0*d2e0dx2;
//	T d2pdro2 = d2p0dro2 + (2.*G+ro*d2Gdro2)*(e-e0) - (G+ro*dGdro)*de0dro - ro*G*d2e0dro2;
//	T d2pdrode = G + ro*dGdro;
//	T d2pde2 = 0.;
//	T KSPrime = (ro*dpdro + ro*ro*d2pdro2 + dpdro*dpde + 2.*d2pdrode - p/ro*dpde - p/ro/ro*dpde*dpde + p*p/ro/ro*d2pde2)/KS;
//	return KSPrime;
//}




//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisenAl<T>::getp(T ro, T e) {
//	const T gamma = 3.9, G = 2.,
//				 B = 76.e9; // Bulk modulus of Al, [Pa];
//	T x = ro/ro0;
//	// p = (pc - G*ec) + G*E
//	// e c = (B/ (gamma (gamma -1 ) )) x^ gamma – (B /( gamma -1 )) x + B/ gamma
//	T p_c = B/gamma * (std::pow(x, gamma)-1.);
//	T e_c = (B/gamma/(gamma-1.)*std::pow(x, gamma) - B*x/(gamma-1.) + B/gamma)/ro;
//	T p  = (p_c-G*ro*e_c) + G*ro*e;
//	return p;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisenAl<T>::gete(T ro, T p) {
//	const T gamma = 3.9, G = 2.,
//				 B = 76.e9; // Bulk modulus of Al, [Pa];
//	T x = ro/ro0;
//	T p_c = B/gamma * (std::pow(x, gamma)-1.);
//	T e_c = (B/gamma/(gamma-1.)*std::pow(x, gamma) - B*x/(gamma-1.) + B/gamma)/ro;
//	T e = (p-p_c+G*ro*e_c)/G/ro;
//	return e;
//}


//template <ArithmeticWith<numeric_val> T>
//T FEOSMieGruneisenAl<T>::getc(T ro, T p) {
//	const T gamma = 3.9, G = 2.,
//				 B = 76.e9; // Bulk modulus of Al, [Pa];
//	T x = ro/ro0;
//	T p_c = B/gamma * (std::pow(x, gamma)-1.);
//	T c  = sqrt((G+1)*(p-p_c)/ro + B/ro0*std::pow(x, gamma-1));
//	return c;
//}
