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
	constexpr static const T ro0 = 2700.;
//	FEOSMieGruneisen(): ro0(998.2) {}     // ro0 = [kg/m3]
	FEOSMieGruneisen() {}     // ro0 = [kg/m3]
//	FEOSMieGruneisen(T _ro0, T _e0) : ro0(_ro0) {}
	FEOSMieGruneisen(T _ro0, T _e0) {}
	static T getp(T ro, T e);
	static T gete(T ro, T p);
	static T getc(T ro, T p);
	static T getG(T ro);
	// Derivatives of G, p0, KS
	static T getGPrime(T ro);
	static T getp0Prime(T ro);
	static T getKSPrime(T ro, T e);
	// Cold components, Born-Meyer potential
	static T getp0(T ro);
	static T gete0(T ro);
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


// =======================================================================
// =======================================================================
// =======================================================================


#define EOS_EPS 1.0e-5
#define EOS_EPS_RO 10.0
#define EOS_EPS_TI 1.0

enum EOSType {table, analytic, ideal, test};


template <ArithmeticWith<numeric_val> T>
class EOSOld {
public:
	virtual EOSType getType()=0;

	virtual T getpi(T ro, T ti) = 0;
	virtual T getpe(T ro, T ti, T te) = 0;
			T getp(T ro, T ti, T te) {
				return getpi(ro, ti) + getpe(ro, ti, te);
			}
	virtual	T getp_e(T ro, T e) = 0;
	virtual	T gete_p(T ro, T p) = 0;

	virtual T getei(T ro, T ti) = 0;
	virtual T getee(T ro, T ti, T te) = 0;
	virtual T gete(T ro, T ti, T te=0.) {
				return getei(ro, ti) + getee(ro, ti, te);
			}

	virtual T getti(T ro, T ei) = 0;
	virtual T gette(T ro, T ti, T ee) = 0;

	virtual T getci(T ro, T ti) = 0;
	virtual T getce(T ro, T te) = 0;

	virtual T getC(T ro, T ti, T te) = 0;
	virtual T getc_p(T ro, T p) = 0;
	virtual T getAlpha(T ro, T ti, T te) = 0;
	virtual T getkappa(T ro, T ti, T te) = 0;

	virtual T getEntropy(T ro, T ti, T te) = 0;
	// Phase features

	virtual T getphase(T ro, T ti) = 0;
	virtual T getmix(T ro, T ti) = 0;

	virtual T getGamma() = 0;

	//////////////////////DEBUG////////////////////////////////////
	T	getC_an(T ro, T ti, T te);
	//////////////////////////////////////////////////////////////


	// Partial derivatives

	virtual T getdpdro(T ro, T ti, T te) = 0;
/*	virtual T getdpdroe(T ro, T ti, T te)=0;
	virtual T getdpdroei(T ro, T ti, T te)=0;*/

	T getdeedte_ro(T ro, T ti, T te);

	virtual T getdpdt(T ro, T ti, T te) = 0;
	virtual T getdedt(T ro, T ti, T te) = 0;

	// dp/dro with ro*v and ro*E held constant
	virtual T getdpdro_rov_roE(
			T ro, T ti, T te, T v) = 0;

	// dp/drov with ro and ro*E held constant
	virtual T getdpdrov_ro_roE(
			T ro, T ti, T te, T v) = 0;

	// dp/droE with ro and ro*v held constant
	virtual T getdpdroE_ro_rov(
			T ro, T ti, T te, T v) = 0;


	// Service and auxilary procedures

	T getMAX_T(void) { return MAX_T; }
	T getMIN_T(void) { return MIN_T; }
	T getMAX_RO(void) { return MAX_RO; }
	T getMIN_RO(void) { return MIN_RO; }
	T getro0(void) { return ro0; }

protected:
	T MAX_T, MIN_T;
	T MAX_RO, MIN_RO;

	T getdpidro_ti(T ro, T ti);
	T getdpidti_ro(T ro, T ti);
	T getdeidro_ti(T ro, T ti);
	T getdeidti_ro(T ro, T ti);

private:
	T ro0;
	T M;
	T gamma;
};


// Analytic EOS for gaseous Al
template <ArithmeticWith<numeric_val> T>
class EOSAnalytic : public EOSOld<T>
{
public:

	EOSAnalytic();

	EOSType getType() { return analytic; }

	T getpi(T ro, T ti);
	T getpe(T ro, T ti, T te);

	T getei(T ro, T ti);
	T getee(T ro, T ti, T te);

	T getti(T ro, T ei);
	T gette(T ro, T ti, T ee);

	T getC(T ro, T ti, T te);
	T getc_p(T ro, T p);

	T getci(T ro, T ti);
	T getce(T ro, T te);
	T getkappa(T ro, T ti, T te);
	T getAlpha(T ro, T ti, T te);

	T getphase(T ro, T ti);
	T getmix(T ro, T ti);

	T getEntropy(T ro, T ti, T te);

	T getdpdro  (T ro, T ti, T te);
	T getdpdroe (T ro, T ti, T te);
	T getdpdroei(T ro, T ti, T te);

	T getdpdt(T ro, T ti, T te);
	T getdedt(T ro, T ti, T te);

	T getdpdro_rov_roE(T ro, T ti, T te, T v);
	T getdpdrov_ro_roE(T ro, T ti, T te, T v);
	T getdpdroE_ro_rov(T ro, T ti, T te, T v);

protected:
	//T z;
	T __ee(T ro, T teta);
	T __pe(T ro, T teta);
	T solve_ti(
			T ro, T ei, T low_border, T high_border);
	// Auxilary partial derivatives
	T dpde(T ro, T ti, T te);
	T dpdei(T ro, T ti, T te);
	T dpdti(T ro, T ti, T te);
	T dedti(T ro, T ti, T te);
	T deidti(T ro, T ti, T te);
	T dpdte(T ro, T ti, T te);
	T dedte(T ro, T ti, T te);

private:
	T ro0;
	T M;
};

#endif // MIEGRUNEISEN_H
