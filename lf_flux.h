#ifndef LF_FLUX_H
#define LF_FLUX_H

#include <algorithm>
#include <execution>
#include <ranges>

#include "arithmeticwith.h"

template <ArithmeticWith<numeric_val> T>
void calcLaxFriedrichsNumericalFlux(
		const std::ranges::common_range auto& u_p,
		const std::ranges::common_range auto& u_m,
		std::ranges::common_range auto& res_f,
		auto&& calcPhysFlux,
		T alpha) {
	/* Global Lax-Friedrichs (LF) numerical flux. LF flux is
	 * attractive because it is extremely simple and easy to
	 * compute, but it is also probably the most dissipative
	 * of all the well-known fluxes, smearing the discontinuities
	 * quite noticeably for lower-order reconstruction methods.
	 *
	 * For high-order reconstructions this drawback becomes much
	 * less noticeable.
	 *
	 * P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations
	 * and Their Numerical Computation, Commun. Pure and Applied
	 * Mathematics, 7, 159-193, 1954.
	 */

	std::transform(std::ranges::begin(u_p), std::ranges::end(u_p),
		// std::ranges::begin(res_f),
		std::ranges::begin(res_f),
		[alpha, &calcPhysFlux](auto u_pt /* , T nf */) {
		return /* nf + */ 0.5 * (calcPhysFlux(u_pt) + alpha * u_pt);
	});

	// for (std::size_t k = 0; k < std::ranges::size(u_p); ++ k) {
	// 	res_f[k] = 0.5 * (calcPhysFlux(u_p[k], k) + alpha * u_p[k]);
	// }

	std::transform(std::ranges::begin(u_m), std::ranges::end(u_m),
		std::ranges::begin(res_f),
		std::ranges::begin(res_f),
		[alpha, &calcPhysFlux](auto u_pt, const auto nf) {
		return nf + 0.5 * (calcPhysFlux(u_pt) - alpha * u_pt);
	});

	// for (std::size_t k = 0; k < std::ranges::size(u_p); ++ k) {
	// 	res_f[k] += 0.5 * (calcPhysFlux(u_m[k], k) + alpha * u_m[k]);
	// }
}


template <ArithmeticWith<numeric_val> T>
struct SplitFluxes {
	std::valarray<T> f_plus;
	std::valarray<T> f_minus;
};


template <ArithmeticWith<numeric_val> T>
SplitFluxes<T> splitFluxAsLaxFriedrichs(
		const std::ranges::common_range auto& u,
		const std::ranges::common_range auto& f,
		T alpha) {
	/* Global Lax-Friedrichs (LF) flux splitting.
	 *
	 * For the purpose of linear stability (upwinding),
	 * a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	 * dfminus/du <= 0), needs to be performed in FD WENO5.
	 * LF flux splitting remains especially simple, while still
	 * retaining the necessary number of derivatives, so it's
	 * perfect for the job (even though LF is a very dissipative
	 * flux, for high-order reconstructions it is not particularly
	 * noticeable).
	 *
	 * P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations
	 * and Their Numerical Computation, Commun. Pure and Applied
	 * Mathematics, 7, 159-193, 1954.
	 */

	SplitFluxes monotone_lf_flux_components {
		std::valarray<T>(std::ranges::size(f)),
		std::valarray<T>(std::ranges::size(f))
	};

	std::transform(/*std::execution::par_unseq,*/
		std::ranges::begin(f), std::ranges::end(f),
		std::ranges::begin(u),
		std::ranges::begin(monotone_lf_flux_components.f_plus),
		[&alpha](T f_pt, T u_pt) {
		return 0.5 * (f_pt + alpha * u_pt);
	});

	std::transform(/*std::execution::par_unseq,*/
		std::ranges::begin(f), std::ranges::end(f),
		std::ranges::begin(u),
		std::ranges::begin(monotone_lf_flux_components.f_minus),
		[&alpha](T f_pt, T u_pt) {
		return 0.5 * (f_pt - alpha * u_pt);
	});

	return monotone_lf_flux_components;
}


//template <ArithmeticWith<numeric_val> T>
//SplitFluxes<T> splitFluxAsRusanovLocalLaxFriedrichs(
//		const std::ranges::common_range auto& u,
//		const std::ranges::common_range auto& f,
//		auto&& calcMaxLambda) {
//	/* Local Lax-Friedrichs (LLF) or Rusanov flux splitting.
//	 *
//	 * For the purpose of linear stability (upwinding),
//	 * a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
//	 * dfminus/du <= 0), needs to be performed in FD WENO5.
//	 * LF flux splitting remains especially simple, while still
//	 * retaining the necessary number of derivatives, so it's
//	 * perfect for the job (even though LF is a very dissipative
//	 * flux, for high-order reconstructions it is not particularly
//	 * noticeable).
//	 *
//	 * P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations
//	 * and Their Numerical Computation, Commun. Pure and Applied
//	 * Mathematics, 7, 159-193, 1954.
//	 */

//// TO DO: we need alpha = max lambda
//// = max{lambda_j, lambda_j+0.5, lambda_j+1}.

//	SplitFluxes monotone_lf_flux_components {
//		std::valarray<T>(std::ranges::size(f)),
//		std::valarray<T>(std::ranges::size(f))
//	};

//	std::transform(std::execution::par_unseq,
//		std::ranges::begin(f), std::ranges::end(f),
//		std::ranges::begin(u),
//		std::ranges::begin(monotone_lf_flux_components.f_plus),
//		[&calcMaxLambda](T f_pt, T u_pt) {
//		return 0.5 * (f_pt + calcMaxLambda(u_pt) * u_pt);
//	});

//	std::transform(std::execution::par_unseq,
//		std::ranges::begin(f), std::ranges::end(f),
//		std::ranges::begin(u),
//		std::ranges::begin(monotone_lf_flux_components.f_minus),
//		[&calcMaxLambda](T f_pt, T u_pt) {
//		return 0.5 * (f_pt - calcMaxLambda(u_pt) * u_pt);
//	});

//	return monotone_lf_flux_components;
//}

#endif // LF_FLUX_H
