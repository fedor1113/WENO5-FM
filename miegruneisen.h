#ifndef MIEGRUNEISEN_H
#define MIEGRUNEISEN_H

#include <string>

#include "arithmeticwith.h"


class FEOS {
public:
	// virtual double getp(double ro, double e) = 0;
	// virtual double gete(double ro, double p) = 0;
	// virtual double getc(double ro, double p) = 0;
	virtual std::string gettype(void) = 0;
	virtual double getdpdrho(double rho, double e) = 0;
	virtual double getdpde(double rho, double e) = 0;
};


class FEOSMieGruneisen : public FEOS {
public:
	const double ro0;
	FEOSMieGruneisen(): ro0(998.2) {}     // ro0 = [kg/m3]
	FEOSMieGruneisen(double _ro0, double _e0) : ro0(_ro0) {}
	double getp(double ro, double e);
	double gete(double ro, double p);
	double getc(double ro, double p);
	double getG(double ro);
	// Derivatives of G, p0, KS
	double getGPrime(double ro);
	double getp0Prime(double ro);
	double getKSPrime(double ro, double e);
	// Cold components, Born-Meyer potential
	double getp0(double ro);
	double gete0(double ro);
	std::string gettype(void) {return std::string("mg"); }
	double getdpdrho(double rho, double e) {return 0.;}
	double getdpde(double rho, double e) {return 0.;}
};


template <ArithmeticWith<numeric_val> T>
class FEOSMieGruneisenAl : public FEOS {
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
