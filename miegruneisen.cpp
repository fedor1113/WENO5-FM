//#include "miegruneisen.h"

//#include <cmath>


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
//	T c2 = 1./ro0*(p0Prime + (GPrime/G + 1./x + G/x)*(p-p0));
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
//	T c  = std::sqrt((G+1.)*(p-p_c)/ro + B/ro0*std::pow(x, gamma-1.));
//	return c;
//}


////template <>
////numeric_val FEOSMieGruneisenAl<numeric_val>::getc(numeric_val ro, numeric_val p) {
////	const numeric_val gamma = 3.9, G = 2.,
////				 B = 76.e9; // Bulk modulus of Al, [Pa];
////	numeric_val x = ro/ro0;
////	numeric_val p_c = B/gamma * (std::pow(x, gamma)-1.);
////	numeric_val c  = std::sqrt((G+1)*(p-p_c)/ro + B/ro0*std::pow(x, gamma-1));
////	return c;
////}


//// =======================================================================
//// =======================================================================
//// =======================================================================


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getpi(T ro, T ti) {
//	T t   = ti/27.2/11605.0; T x   = 0.64/(ro/1000.0);
//	T m1  = 0.00001; T m2  = 0.000004; T m3  = 0.00002; T m4  = 0.000002;
//	T tau = (m1/std::pow(x, 1.0/3.0) + m2/std::pow(x, 2.0/3.0) + m3/x + m4/std::pow(x, 4.0/3.0)) / t;
//	T tau1 =  -1.0/3.0 / t * (m1/std::pow(x, 4.0/3.0) + 2*m2/std::pow(x, 5.0/3.0) + 3*m3/(x*x) + 4*m4/std::pow(x, 7.0/3.0));
//	// p1
//	T alpha = 2.98; T beta  = 2.2; T A     = 4.483e-5;T B     = 1.378e-4;
//	T p1 = A/std::pow(x, alpha) - B/std::pow(x, beta);
//	// pi
//	T dzeta = 3.69; T d1    = 95.0163; T d2    = 440.63; T nc    = 0.237/111.8;
//	T pi = (-dzeta*nc*t*(1.0/tau + d1 - 2*d2*tau) * tau1 + p1) * 29400.0e10;
//	return pi / 10.0;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getpe(T ro, T ti, T te) {
//	T pe = __pe(ro, te) - __pe(ro, ti); return pe / 10.0;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::__pe(T ro, T teta) {
//	T NA    = 6.0e23; T kB    = 1.38e-16; T Matom = 27.0/NA; T nat   = ro/1000.0/Matom;
//	T eF    = 115000.0; T pF    = 2.0/5.0*3.*nat*kB*eF; T pCl   = 3.*nat*kB*teta;
//	T pe = std::sqrt(pF*pF + pCl*pCl);
//	return pe;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getei(T ro, T ti) {
//	// tau
//	T t   = ti/27.2/11605.0; T x   = 0.64/(ro/1000.0);
//	T m1  = 0.00001; T m2  = 0.000004; T m3  = 0.00002; T m4  = 0.000002;
//	T tau = (m1/std::pow(x, 1.0/3.0) + m2/std::pow(x, 2.0/3.0) + m3/x + m4/std::pow(x, 4.0/3.0)) / t;
//	// en0
//	T en0 = -1.05835;
//	// en1
//	T alpha = 2.98; T beta  = 2.2; T A     = 4.483e-5; T B     = 1.378e-4;
//	T en1 = A/(alpha-1)/std::pow(x, alpha-1) - B/(beta-1)/std::pow(x, beta-1);
//	// ei
//	T dzeta = 3.69; T d1    = 95.0163; T d2    = 440.63; T nc    = 0.237/111.8;
//	T ei = (en0 + t*dzeta*(1.0+d1*tau-2.0*d2*tau*tau) * 27.2 +	en1/nc*27.2) * 1.6e-12 * 6.0e23/27.0 * ro/1000.0;
//	return ei / (ro/1000.0) * 1.0e-4;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getee(T ro, T ti, T te) {
//	T ee = __ee(ro, te) - __ee(ro, ti);
//	return ee / (ro/1000.0) * 1.0e-4;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::__ee(T ro, T teta) {
//	T NA    = 6.0e23; T kB    = 1.38e-16; T Matom = 27.0/NA; T nat   = ro/1000.0/Matom;
////	Old gamma!
////	T gamma = 1450.0*std::pow(ro/2700.0, 1.0/3.0);
//	T   TF  = 115000.0; T   EF  = kB * TF; T gamma = 3.14159*3.14159*3.*nat*kB*kB/2.0/EF;
//	T ceCl  = 3.0/2.0*3.*nat*kB; T ee = ceCl/gamma * (std::sqrt(ceCl*ceCl + gamma*gamma*teta*teta) - ceCl);
//	return ee;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getti(T ro, T ei) {
//	T ti = solve_ti(ro, ei, 0.00001, 60000.0); return ti;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::gette(T ro, T ti, T ee) {
//	T NA    = 6.0e23; T kB    = 1.38e-16; T Matom = 27.0/NA; T nat   = ro/1000.0/Matom;
//// Old gamma
//// T gamma = 1450.0*std::pow(ro/2700.0, 1.0/3.0);
//	T   TF  = 115000.0; T   EF  = kB * TF;
//	T gamma = 3.14159*3.14159*3.*nat*kB*kB/2.0/EF; T ceCl  = 3.0/2.0*3.*nat*kB;
//	T A = std::sqrt(ceCl*ceCl+gamma*gamma*ti*ti) + gamma/ceCl*ee*(ro/1000.0)/1.0e-4;;
//	T te = 1.0/gamma * std::sqrt( A * A - ceCl * ceCl );
//	return te;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::solve_ti(T ro, T ei, T low_border, T high_border) {
//	T ti; T eps = 0.01; T ti_dei, low_dei, high_dei;
//	for(;;) {
//		ti = (low_border + high_border)/2.0;
//		if(std::abs(low_border-high_border) < eps)
//			return ti;
//		ti_dei   = getei(ro, ti)-ei;
//		low_dei  = getei(ro, low_border)-ei;
//		high_dei = getei(ro, high_border)-ei;
//		if(std::abs(ti_dei) < eps)
//			return ti;
//		if(low_dei * high_dei > 0) {
//			printf("\nsolve_ti: Error in initial data -- no equation root on the interval\n");
//			printf("ro=%e, ei=%e, F(a)=%e, F(b)=%e\n", ro, ei, low_dei, high_dei);
//			return -1.;
//		}
//		if(ti_dei * low_dei > 0)
//			low_border  += (high_border-low_border)/2.0;
//		else
//			high_border -= (high_border-low_border)/2.0;
//	}
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getC(T ro, T ti, T te) {
//	return std::sqrt( getdpdro(ro, ti, te) +
//				(getpi(ro, ti)+getpe(ro, ti, te)) * dpde(ro, ti, te) / ro / ro +
//				 getpi(ro, ti) * dpdei(ro, ti, te) / ro / ro );
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getc_p(T ro, T p) {
//	return 1000.0;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getci(T ro, T ti) {
///*	T NA    = 6.0e23;
//	T kB    = 1.38e-16;
//	T Matom = 27.0/NA;
//	T nat   = ro/1000.0/Matom;
//	T ci = 3.0*nat*kB;
//	return ci*1.0e-1; */
///*	if(ro!=2700.) {
//		cout << "Density is not equal to 2700 kg/m3. Possibly hydrodynamic model icluded. Try another approximation of ion heat capacity ci." << endl;
//		exit(1);
//	} */
//	T c1 = 0.002484e9, c2 = 0.00595e9, c3 = 0.0019e9;
//	T t1_minus = 1150., t1_plus = 1250., t2_minus = 1550., t2_plus = 1650.;
//	if(ti < t1_minus) return c1;
//	else if (ti<t1_plus) return c1 + (c2-c1)/(t1_plus-t1_minus)*(ti-t1_minus);
//	else if(ti<t2_minus) return c2;
//	else if (ti<t2_plus) return c2 + (c3-c2)/(t2_plus-t2_minus)*(ti-t2_minus);
//	else return c3;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getce(T ro, T te) {
//	//	ce(rho, T)=( ( gam T)^( - 2 ) + cCle^( - 2) )^( - 1/2)
//	T NA    = 6.0e23; T kB    = 1.38e-16; T Matom = 27.0/NA; T nat   = ro/1000.0/Matom;
//	//	T gamma = 1450.0*std::pow(ro/2700.0, 1.0/3.0);
//	T   TF  = 115000.0; T   EF  = kB * TF;
//	T gamma = 3.14159*3.14159*3.*nat*kB*kB/2.0/EF;
//	T ceCl  = 3.0/2.0*3.*nat*kB;
//	T ce = std::pow(1.0/gamma/gamma/te/te + 1.0/ceCl/ceCl, -0.5);
//	return ce * 1.0e-1;
//	// 4 testHeatStage()    : return 1.0;
//	// 4 testExchangeStage(): return 1.0;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getkappa(T ro, T ti, T te) {
//	T C      = 7.25; T beta   = 0.504; T TF     = 115000.0;
//	T theta  = te/TF; T thetai = ti/TF;
//	T sk54   = std::pow(theta*theta + 0.16, 1.25); T sk     = theta*theta + 0.44;
//	T sk12   = std::pow(theta*theta + 0.092, 0.5);
//	T kappa  = C*sk54*sk/sk12*theta / (theta*theta+beta*thetai) * 1.0e7;
//	return kappa * 1.0e-5;
//	// 4 heatTest(): return 1.0e-8;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getAlpha(T ro, T ti, T te) {
//	return -36.0e17 * 1.0e-1 * ro/ro0 ;
//	// 4 testExchangeStage() : return -1.0;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getdpdro(T ro, T ti, T te) {
//	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te)
//		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getdpdroe(T ro, T ti, T te) {
//	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te)
//		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getdpdroei(T ro, T ti, T te) {
//	return ( getpi(ro + EOS_EPS, ti) + getpe(ro + EOS_EPS, ti, te)
//		   - getpi(ro, ti)			 - getpe(ro, ti, te) ) / EOS_EPS;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getdpdro_rov_roE(T ro, T ti, T te, T v)
//{ return 0; }


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getdpdrov_ro_roE(T ro, T ti, T te, T v)
//{ return 0; }


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getdpdroE_ro_rov(T ro, T ti, T te, T v)
//{ return 0; }


//////////////// PRIVATE ////////////////

//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::dpdei(T ro, T ti, T te) {
//	return dpdti(ro, ti, te) / deidti(ro, ti, te) -
//		   dpdte(ro, ti, te) / dedte(ro, ti, te);
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::dpde(T ro, T ti, T te) {
//	return dpdte(ro, ti, te) / dedte(ro, ti, te);
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::dpdti(T ro, T ti, T te) {
//	return ( getpi(ro, ti + EOS_EPS) + getpe(ro, EOS_EPS, te)
//		   - getpi(ro, ti)           - getpe(ro, ti, te) ) / EOS_EPS;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::dedti(T ro, T ti, T te) {
//	return ( getei(ro, ti + EOS_EPS) + getee(ro, ti + EOS_EPS, te)
//		   - getei(ro, ti)           - getee(ro, ti, te) ) / EOS_EPS;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::deidti(T ro, T ti, T te) {
//	return ( getei(ro, ti + EOS_EPS) - getei(ro, ti) ) / EOS_EPS;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::dpdte(T ro, T ti, T te) {
//	return ( getpe(ro, ti, te + EOS_EPS) - getpe(ro, ti, te) ) / EOS_EPS;

//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::dedte(T ro, T ti, T te) {
//	return ( getee(ro, ti, te + EOS_EPS) - getee(ro, ti, te) ) / EOS_EPS;
//}


//template class EOSAnalytic<numeric_val>;


//template <>
//EOSAnalytic<numeric_val>::EOSAnalytic() {
//	MAX_T=100000.0;
//	MIN_T=300.0;
//	MAX_RO=5000.0;
//	MIN_RO=0.001;
//	ro0 = 2700.;  // [kg/m3]
//	M = 27.e-3;   // [kg/mol]
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getphase(T ro, T ti) {
//	return 1.0;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getmix(T ro, T ti) {
//	return 1.0;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getEntropy(T ro, T ti, T te) {
//	return 0.0;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getdpdt(T ro, T ti, T te) {
//	return (getpi(ro, ti+EOS_EPS)-getpi(ro, ti))/EOS_EPS;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSAnalytic<T>::getdedt(T ro, T ti, T te) {
//	return (getei(ro, ti+EOS_EPS)-getei(ro, ti))/EOS_EPS;
//}




//template <ArithmeticWith<numeric_val> T>
//T EOSOld<T>::getC_an(T ro, T ti, T te) {
//	T
//			   pi =        getpi(ro, ti),
//			pi_ro = getdpidro_ti(ro ,ti),
//			pi_ti = getdpidti_ro(ro ,ti),
//			ei_ro = getdeidro_ti(ro ,ti),
//			ei_ti = getdeidti_ro(ro ,ti);

//	//////////////DEBUG/////////////////
//	T c = std::sqrt(std::abs((pi_ro - pi_ti * (ei_ro - pi/ro/ro) / ei_ti)));


//	T c2=std::sqrt(5.3 * getpi(ro, ti)/ro);

//	////////////////////////////////////



//	return c;
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSOld<T>::getdpidro_ti(T ro, T ti) {
//	T h = 0.1;

//	return ( getpi(ro*(1.0+h), ti)-getpi(ro, ti) ) / (ro*h);
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSOld<T>::getdpidti_ro(T ro, T ti) {
//	T h = 0.1;
//	return ( getpi(ro, ti*(1.0+h))-getpi(ro, ti) ) / (ti*h);
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSOld<T>::getdeidro_ti(T ro, T ti) {
//	T h = 0.1;

//	return ( getei(ro*(1.0+h), ti)-getei(ro, ti) ) / (ro*h);
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSOld<T>::getdeidti_ro(T ro, T ti) {
//	T h = 0.1;
//	return ( getei(ro, ti*(1.0+h))-getei(ro, ti) ) / (ti*h);
//}


//template <ArithmeticWith<numeric_val> T>
//T EOSOld<T>::getdeedte_ro(T ro, T ti, T te) {
//	T h = 0.1;
//	return ( getee(ro, ti, te*(1.0+h))-getee(ro, ti, te) ) / (te*h);
//}
