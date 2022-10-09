#ifndef MIEGRUNEISEN_H
#define MIEGRUNEISEN_H

#include <string>

#include "arithmeticwith.h"


template <ArithmeticWith<numeric_val> T>
class FEOS {
public:
	// virtual T getp(T ro, T e) = 0;
	// virtual T gete(T ro, T p) = 0;
	// virtual T getc(T ro, T p) = 0;
	virtual std::string gettype(void) = 0;
	virtual T getdpdrho(T rho, T e) = 0;
	virtual T getdpde(T rho, T e) = 0;
};


template <ArithmeticWith<numeric_val> T>
class FEOSMieGruneisen : public FEOS<T> {
public:
	const T ro0;
	FEOSMieGruneisen(): ro0(998.2) {}     // ro0 = [kg/m3]
	FEOSMieGruneisen(T _ro0, T _e0) : ro0(_ro0) {}
	T getp(T ro, T e);
	T gete(T ro, T p);
	T getc(T ro, T p);
	T getG(T ro);
	// Derivatives of G, p0, KS
	T getGPrime(T ro);
	T getp0Prime(T ro);
	T getKSPrime(T ro, T e);
	// Cold components, Born-Meyer potential
	T getp0(T ro);
	T gete0(T ro);
	std::string gettype(void) {return std::string("mg"); }
	T getdpdrho(T rho, T e) {return 0.;}
	T getdpde(T rho, T e) {return 0.;}
};


template <ArithmeticWith<numeric_val> T>
class FEOSMieGruneisenAl : public FEOS<T> {
public:
	constexpr static const T ro0 = 2700.;
	// FEOSMieGruneisenAl(): /*ro0(2700.)*/ {}     // ro0 = [kg/m3]
	static T getp(T ro, T e);
	static T gete(T ro, T p);
	static T getc(T ro, T p);
	std::string gettype(void) {return std::string("mg"); }
};


template class FEOSMieGruneisenAl<numeric_val>;

#endif // MIEGRUNEISEN_H
