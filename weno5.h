//#pragma GCC optimize("Ofast")
//#pragma GCC target("avx,avx2,fma")
//#pragma GCC optimize("unroll-loops")

#ifndef WENO5_H
#define WENO5_H

// #include "weno5coefs.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <ctime>
#include <execution>
// #include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <ranges>
// #include <span>
// #include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

#include "Eigen/Dense"

#include "_vector4.h"
#include "arithmeticwith.h"

// template class Vector4<double>;

// using Vec4 = Vector4<double>;


template <ArithmeticWith<numeric_val> T, typename T1, typename T2>
std::valarray<T> vecMatDot(const T1& vec, const T2& mat) {
	/* Multiply a vector by a matrix (a sum product over
	 * the last axis of mat and vec).
	 */

	std::valarray<T> result(std::ranges::size(mat));

	for (std::size_t k = 0; k < std::ranges::size(mat)
			&& k < std::ranges::size(vec); ++ k)
		result[k] = std::inner_product(std::ranges::begin(mat[k]),
									   std::ranges::end(mat[k]),
									   std::ranges::begin(vec), 0.);

	return result;
}


//template <ArithmeticWith<numeric_val> T>
////std::array<T, 3> smoothness_indicators(const T1& f_stencil) {
//std::valarray<T> betaSmoothnessIndicators(
//		const std::ranges::sized_range auto& f_stencil) {
//	/* Return the WENO5 smoothness indicators of Jiang and Shu (1996)
//	 * for each of the 3 substencils.
//	 * That is the sum of the normalized squares of the scaled
//	 * L2-norms of all the derivatives of 3 local interpolating
//	 * polynomials in the sub-stencils of 5-node `f_stencil`.
//	 *
//	 * This allows for (2*3-1)=5th order accuracy from the 3rd
//	 * order Eno schemes.
//	 */

//	// std::array<T, 3> res;
//	std::valarray<T> betas(3);

//	betas[0] = ((13./12.) * std::pow(
//					f_stencil[0] - 2.*f_stencil[1] + f_stencil[2], 2)
//				+ (1./4.) * std::pow(
//					f_stencil[0] - 4.*f_stencil[1] + 3.*f_stencil[2], 2
//					));

//	betas[1] = ((13./12.) * std::pow(
//					f_stencil[1] - 2.*f_stencil[2] + f_stencil[3], 2)
//				+ (1./4.) * std::pow((f_stencil[1] - f_stencil[3]), 2));

//	betas[2] = ((13./12.) * std::pow(
//					f_stencil[2] - 2.*f_stencil[3] + f_stencil[4], 2)
//				+ (1./4.) * std::pow(
//					3.*f_stencil[2] - 4.*f_stencil[3] + f_stencil[4], 2
//					));

//	return betas;
//}


template <ArithmeticWith<numeric_val> T>
//std::array<T, 3> smoothness_indicators(const T1& f_stencil) {
std::valarray<T> betaSmoothnessIndicators(
		const std::ranges::sized_range auto& f_stencil) {
	/* Return the WENO5 smoothness indicators of Jiang and Shu (1996)
	 * for each of the 3 substencils.
	 * That is the sum of the normalized squares of the scaled
	 * L2-norms of all the derivatives of 3 local interpolating
	 * polynomials in the sub-stencils of 5-node `f_stencil`.
	 *
	 * This allows for (2*3-1)=5th order accuracy from the 3rd
	 * order Eno schemes.
	 */

	// std::array<T, 3> res;
	std::valarray<T> betas(3);

	static const Eigen::DiagonalMatrix<T, 3> l0(
				2. * std::sqrt(3.) / 3.,
				std::sqrt(13.) / 4.,
				0.);
	static const Eigen::Matrix<T, 3, 3> c0 {
		{1., -19. / 8., 11. / 8.},
		{0.,        1.,      -1.},
		{0.,        0.,       0.}
	};
	static const Eigen::Matrix<T, 3, 3> b0 = l0 * c0;

	static const Eigen::DiagonalMatrix<T, 3> l1(
				2. * std::sqrt(3.) / 3.,
				std::sqrt(13.) / 4.,
				0.);
	static const Eigen::Matrix<T, 3, 3> c1 {
		{1., -13. / 8., 5. / 8.},
		{0.,        1.,     -1.},
		{0.,        0.,      0.}
	};
	static const Eigen::Matrix<T, 3, 3> b1 = l1 * c1;

	static const Eigen::DiagonalMatrix<T, 3> l2(
				std::sqrt(30.) / 3.,
				std::sqrt(130.) / 20.,
				0.);
	static const Eigen::Matrix<T, 3, 3> c2 {
		{1., -31. / 20., 11. / 20.},
		{0.,         1.,       -1.},
		{0.,         0.,        0.}
	};
	static const Eigen::Matrix<T, 3, 3> b2 = l2 * c2;

	Eigen::Matrix<T, 3, 1> substencil0 {
		f_stencil[0], f_stencil[1], f_stencil[2]
	};
	Eigen::Matrix<T, 3, 1> substencil1 {
		f_stencil[1], f_stencil[2], f_stencil[3]
	};
	Eigen::Matrix<T, 3, 1> substencil2 {
		f_stencil[2], f_stencil[3], f_stencil[4]
	};

	betas[0] = ((b0 * substencil0).transpose())
			* (b0 * substencil0);

	betas[1] = ((b1 * substencil1).transpose())
			* (b1 * substencil1);

	betas[2] = ((b2 * substencil2).transpose())
			* (b2 * substencil2);

	return betas;
}


//template <ArithmeticWith<numeric_val> T>
////std::array<T, 3> smoothness_indicators(const T1& f_stencil) {
//std::valarray<T> betaSmoothnessIndicatorsWENO7BS(
//		const std::ranges::sized_range auto& f_stencil) {
//	/* Return the WENO7 smoothness indicators of Balsara and Shu (2000)
//	 * for each of the 4 substencils.
//	 * That is the sum of the normalized squares of the scaled
//	 * L2-norms of all the derivatives of 4 local interpolating
//	 * polynomials in the sub-stencils of 7-node `f_stencil`.
//	 *
//	 * This allows for (2*4-1)=8th order accuracy from the 4th
//	 * order ENO schemes.
//	 */

//	// std::array<T, 3> res;
//	std::valarray<T> betas(4);

////	betas[0] = f_stencil[2] * (
////				134241. * f_stencil[2]
////				- 114894. * f_stencil[3])
////			+ f_stencil[0] * (
////				56694. * f_stencil[2]
////				- 47214. * f_stencil[1]
////				+ 6649. * f_stencil[0]
////				- 22778. * f_stencil[3])
////			+ 25729. * f_stencil[3] * f_stencil[3]
////			+ f_stencil[1] * (
////				-210282. * f_stencil[2]
////				+ 85641. * f_stencil[1]
////				+ 86214. * f_stencil[3]);

////	betas[1] = f_stencil[3] * (
////				41001. * f_stencil[3]
////				- 30414. * f_stencil[4])
////			+ f_stencil[1] * (
////				-19374. * f_stencil[2]
////				+ 3169. * f_stencil[1]
////				+ 19014. * f_stencil[3]
////				- 5978. * f_stencil[4])
////			+ 6649. * f_stencil[4] * f_stencil[4]
////			+ f_stencil[2] * (
////				33441. * f_stencil[2]
////				- 70602. * f_stencil[3]
////				+ 23094. * f_stencil[4]);

////	betas[2] = f_stencil[4] * (
////				33441. * f_stencil[4]
////				- 19374. * f_stencil[5])
////			+ f_stencil[2] * (
////				6649. * f_stencil[2]
////				- 30414. * f_stencil[3]
////				+ 23094. * f_stencil[4]
////				- 5978. * f_stencil[5])
////			+ 3169. * f_stencil[5] * f_stencil[5]
////			+ f_stencil[3] * (
////				41001. * f_stencil[3]
////				- 70602. * f_stencil[4]
////				+ 19014. * f_stencil[5]);

////	betas[3] = f_stencil[5] * (
////				85641. * f_stencil[5]
////				- 47214. * f_stencil[6])
////			+ f_stencil[3] * (
////				25729. * f_stencil[3]
////				- 114894. * f_stencil[4]
////				+ 86214. * f_stencil[5]
////				- 22778. * f_stencil[6])
////			+ 6649. * f_stencil[6] * f_stencil[6]
////			+ f_stencil[4] * (
////				134241. * f_stencil[4]
////				- 210282. * f_stencil[5]
////				+ 56694. * f_stencil[6]);

//	betas[0] = f_stencil[0] * (
//				547. * f_stencil[0]
//				- 3882. * f_stencil[1]
//				+ 4642. * f_stencil[2]
//				- 1854. * f_stencil[3])
//			+ f_stencil[1] * (
//				7043. * f_stencil[1]
//				- 17246. * f_stencil[2]
//				+ 7042. * f_stencil[3])
//			+ f_stencil[2] * (
//				11003. * f_stencil[2]
//				- 9402. * f_stencil[3])
//			+ 2107. * f_stencil[3] * f_stencil[3];

//	betas[1] = f_stencil[1] * (
//				267. * f_stencil[1]
//				- 1642. * f_stencil[2]
//				+ 1602. * f_stencil[3]
//				- 494. * f_stencil[4])
//			+ f_stencil[2] * (
//				2843. * f_stencil[2]
//				- 5966. * f_stencil[3]
//				+ 1922. * f_stencil[4])
//			+ f_stencil[3] * (
//				3443. * f_stencil[3]
//				- 2522. * f_stencil[4])
//			+ 547. * f_stencil[4] * f_stencil[4];

//	betas[2] = f_stencil[2] * (
//				547. * f_stencil[2]
//				- 2522. * f_stencil[3]
//				+ 1922. * f_stencil[4]
//				- 494. * f_stencil[5])
//			+ f_stencil[3] * (
//				3443. * f_stencil[3]
//				- 5966. * f_stencil[4]
//				+ 1602. * f_stencil[5])
//			+ f_stencil[4] * (
//				2843. * f_stencil[4]
//				- 1642. * f_stencil[5])
//			+ 267. * f_stencil[5] * f_stencil[5];

//	betas[3] = f_stencil[3] * (
//				2107. * f_stencil[3]
//				- 9402. * f_stencil[4]
//				+ 7042. * f_stencil[5]
//				- 1854. * f_stencil[6])
//			+ f_stencil[4] * (
//				11003. * f_stencil[4]
//				- 17246. * f_stencil[5]
//				+ 4642. * f_stencil[6])
//			+ f_stencil[5] * (
//				7043. * f_stencil[5]
//				- 3882. * f_stencil[6])
//			+ 547. * f_stencil[6] * f_stencil[6];

//	return betas;
//}


//template <ArithmeticWith<numeric_val> T>
//std::valarray<T> betaSmoothnessIndicatorsWENO7BS(
//		const std::ranges::sized_range auto& f_stencil) {
//	/* Return the WENO7 smoothness indicators of Balsara and Shu (2000)
//	 * for each of the 4 substencils.
//	 * That is the sum of the normalized squares of the scaled
//	 * L2-norms of all the derivatives of 4 local interpolating
//	 * polynomials in the sub-stencils of 7-node `f_stencil`.
//	 *
//	 * This allows for (2*4-1)=8th order accuracy from the 4th
//	 * order ENO schemes.
//	 */

//	std::valarray<T> betas(4);

//	static const Eigen::DiagonalMatrix<T, 4> l0(
//				std::sqrt(8205.) / 60.,
//				std::sqrt(1744383.) / 1641.,
//				std::sqrt(32377917.) / 6378.,
//				0.);
//	static const Eigen::Matrix<T, 4, 4> c0 {
//		{1., -1941./547.,    2321./547.,  -927./547.},
//		{0.,          1.,  -5293./2126., 3167./2126.},
//		{0.,          0.,            1.,         -1.},
//		{0.,          0.,            0.,          0.}
//	};
//	static const Eigen::Matrix<T, 4, 4> b0 = l0 * c0;

//	static const Eigen::DiagonalMatrix<T, 4> l1(
//				std::sqrt(445.) / 20.,
//				std::sqrt(94607.) / 267.,
//				std::sqrt(32377917.) / 6378.,
//				0.);
//	static const Eigen::Matrix<T, 4, 4> c1 {
//		{1.,  -821./267.,            3.,  -821./267.},
//		{0.,          1.,  -3471./2126., 1345./2126.},
//		{0.,          0.,            1.,         -1.},
//		{0.,          0.,            0.,          0.}
//	};
//	static const Eigen::Matrix<T, 4, 4> b1 = l1 * c1;

//	static const Eigen::DiagonalMatrix<T, 4> l2(
//				std::sqrt(8205.) / 60.,
//				std::sqrt(6014265.) / 1641.,
//				std::sqrt(111632235.) / 21990.,
//				0.);
//	static const Eigen::Matrix<T, 4, 4> c2 {
//		{1., -1261./547.,     961./547.,  -247./547.},
//		{0.,          1., -10497./7330., 3167./7330.},
//		{0.,          0.,            1.,         -1.},
//		{0.,          0.,            0.,          0.}
//	};
//	static const Eigen::Matrix<T, 4, 4> b2 = l2 * c2;

//	static const Eigen::DiagonalMatrix<T, 4> l3(
//				7. * std::sqrt(645.) / 60.,
//				std::sqrt(1747821.) / 903.,
//				std::sqrt(412688991.) / 81294.,
//				0.);
//	static const Eigen::Matrix<T, 4, 4> c3 {
//		{1., -4701./2107.,      503./301.,   -927./2107.},
//		{0.,           1., -40411./27098., 13313./27098.},
//		{0.,           0.,             1.,           -1.},
//		{0.,           0.,             0.,            0.}
//	};
//	static const Eigen::Matrix<T, 4, 4> b3 = l3 * c3;

//	Eigen::Matrix<T, 4, 1> substencil0 {
//		f_stencil[0], f_stencil[1], f_stencil[2], f_stencil[3]
//	};
//	Eigen::Matrix<T, 4, 1> substencil1 {
//		f_stencil[1], f_stencil[2], f_stencil[3], f_stencil[4]
//	};
//	Eigen::Matrix<T, 4, 1> substencil2 {
//		f_stencil[2], f_stencil[3], f_stencil[4], f_stencil[5]
//	};
//	Eigen::Matrix<T, 4, 1> substencil3 {
//		f_stencil[3], f_stencil[4], f_stencil[5], f_stencil[6]
//	};

//	betas[0] = ((b0 * substencil0).transpose())
//			* (b0 * substencil0);

//	betas[1] = ((b1 * substencil1).transpose())
//			* (b1 * substencil1);

//	betas[2] = ((b2 * substencil2).transpose())
//			* (b2 * substencil2);

//	betas[3] = ((b3 * substencil3).transpose())
//			* (b3 * substencil3);

//	return betas;
//}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> betaSmoothnessIndicatorsWENO7BS(
		const std::ranges::sized_range auto& f_stencil) {
	/* Return the WENO7 smoothness indicators of Balsara and Shu (2000)
	 * for each of the 4 substencils.
	 * That is the sum of the normalized squares of the scaled
	 * L2-norms of all the derivatives of 4 local interpolating
	 * polynomials in the sub-stencils of 7-node `f_stencil`.
	 *
	 * This allows for (2*4-1)=8th order accuracy from the 4th
	 * order ENO schemes.
	 */

	std::valarray<T> betas(4);

	static const Eigen::Matrix<T, 3, 4> g0 {
		{-1./3., 3./2., -3., 11./6.},
		{   -1.,    4., -5.,     2.},
		{   -1.,    3., -3.,     1.}
	};

	static const Eigen::Matrix<T, 3, 4> g1 {
		{1./6., -1., 1./2., 1./3.},
		{   0.,  1.,   -2.,    1.},
		{  -1.,  3.,   -3.,    1.}
	};

	static const Eigen::Matrix<T, 3, 4> g2 {
		{-1./3., -1./2.,  1., -1./6.},
		{    1.,    -2.,  1.,     0.},
		{   -1.,     3., -3.,     1.}
	};

	static const Eigen::Matrix<T, 3, 4> g3 {
		{-11./6.,  3.,  -3./2., 1./3.},
		{     2., -5.,      4.,   -1.},
		{    -1.,  3.,     -3.,    1.}
	};

	Eigen::Matrix<T, 4, 1> substencil0 {
		f_stencil[0], f_stencil[1], f_stencil[2], f_stencil[3]
	};
	Eigen::Matrix<T, 4, 1> substencil1 {
		f_stencil[1], f_stencil[2], f_stencil[3], f_stencil[4]
	};
	Eigen::Matrix<T, 4, 1> substencil2 {
		f_stencil[2], f_stencil[3], f_stencil[4], f_stencil[5]
	};
	Eigen::Matrix<T, 4, 1> substencil3 {
		f_stencil[3], f_stencil[4], f_stencil[5], f_stencil[6]
	};

	Eigen::Matrix<T, 3, 1> v0 = g0 * substencil0;
	Eigen::Matrix<T, 3, 1> v1 = g1 * substencil1;
	Eigen::Matrix<T, 3, 1> v2 = g2 * substencil2;
	Eigen::Matrix<T, 3, 1> v3 = g3 * substencil3;

	betas[0] = v0[0] * v0[0]
			+ 13./12. * v0[1] * v0[1]
			+ 781./720. * v0[2] * v0[2];

	betas[1] = v1[0] * v1[0]
			+ 13./12. * v1[1] * v1[1]
			+ 781./720. * v1[2] * v1[2];

	betas[2] = v2[0] * v2[0]
			+ 13./12. * v2[1] * v2[1]
			+ 781./720. * v2[2] * v2[2];

	betas[3] = v3[0] * v3[0]
			+ 13./12. * v3[1] * v3[1]
			+ 781./720. * v3[2] * v3[2];

	return betas;
}


template <ArithmeticWith<numeric_val> T>
//std::array<T, 3> smoothness_indicators(const T1& f_stencil) {
std::valarray<T> betaSmoothnessIndicatorsWENO9BS(
		const std::ranges::sized_range auto& f_stencil) {
	/* Return the WENO9 smoothness indicators of Balsara and Shu (2000)
	 * for each of the 5 substencils.
	 * That is the sum of the normalized squares of the scaled
	 * L2-norms of all the derivatives of 5 local interpolating
	 * polynomials in the sub-stencils of 9-node `f_stencil`.
	 *
	 * This allows for (2*5-1)=9th order accuracy from the 5th
	 * order ENO schemes.
	 */

	// std::array<T, 3> res;
	std::valarray<T> betas(5);

	betas[0] = f_stencil[0] * (
				   22658. * f_stencil[0]
				- 208501. * f_stencil[1]
				+ 364863. * f_stencil[2]
				- 288007. * f_stencil[3]
				+ 86329.  * f_stencil[4])
			+ f_stencil[1] * (
				   482963. * f_stencil[1]
				- 1704396. * f_stencil[2]
				+ 1358458. * f_stencil[3]
				- 411487.  * f_stencil[4])
			+ f_stencil[2] * (
				  1521393. * f_stencil[2]
				- 2462076. * f_stencil[3]
				+ 758823.  * f_stencil[4])
			+ f_stencil[3] * (
				 1020563. * f_stencil[3]
				- 649501. * f_stencil[4])
			+ 107918. * f_stencil[4] * f_stencil[4];

	betas[1] = f_stencil[1] * (
				  6908.  * f_stencil[1]
				- 60871. * f_stencil[2]
				+ 99213. * f_stencil[3]
				- 70237. * f_stencil[4]
				+ 18079. * f_stencil[5])
			+ f_stencil[2] * (
				  138563. * f_stencil[2]
				- 464976. * f_stencil[3]
				+ 337018. * f_stencil[4]
				- 88297.  * f_stencil[5])
			+ f_stencil[3] * (
				406293. * f_stencil[3]
				- 611976. * f_stencil[4]
				+ 165153. * f_stencil[5])
			+ f_stencil[4] * (
				  242723. * f_stencil[4]
				- 140251. * f_stencil[5])
			+ 22658. * f_stencil[5] * f_stencil[5];

	betas[2] = f_stencil[2] * (
				6908. * f_stencil[2]
				- 51001. * f_stencil[3]
				+ 67923. * f_stencil[4]
				- 38947. * f_stencil[5]
				+ 8209.  * f_stencil[6])
			+ f_stencil[3] * (
				104963. * f_stencil[3]
				- 299076. * f_stencil[4]
				+ 179098. * f_stencil[5]
				- 38947.  * f_stencil[6])
			+ f_stencil[4] * (
				231153. * f_stencil[4]
				- 299076. * f_stencil[5]
				+ 67923.  * f_stencil[6])
			+ f_stencil[5] * (
				 104963. * f_stencil[5]
				- 51001. * f_stencil[6])
			+ 6908. * f_stencil[6] * f_stencil[6];

	betas[3] = f_stencil[3] * (
				  22658.  * f_stencil[3]
				- 140251. * f_stencil[4]
				+ 165153. * f_stencil[5]
				- 88297.  * f_stencil[6]
				+ 18079.  * f_stencil[7])
			+ f_stencil[4] * (
				  242723. * f_stencil[4]
				- 611976. * f_stencil[5]
				+ 337018. * f_stencil[6]
				- 70237.  * f_stencil[7])
			+ f_stencil[5] * (
				  406293. * f_stencil[5]
				- 464976. * f_stencil[6]
				+ 99213.  * f_stencil[7])
			+ f_stencil[6] * (
				 138563. * f_stencil[6]
				- 60871. * f_stencil[7])
			+ 6908. * f_stencil[7] * f_stencil[7];

	betas[4] = f_stencil[4] * (
				  107918. * f_stencil[4]
				- 649501. * f_stencil[5]
				+ 758823. * f_stencil[6]
				- 411487. * f_stencil[7]
				+ 86329.  * f_stencil[8])
			+ f_stencil[5] * (
				  1020563. * f_stencil[5]
				- 2462076. * f_stencil[6]
				+ 1358458. * f_stencil[7]
				- 288007.  * f_stencil[8])
			+ f_stencil[6] * (
				  1521393. * f_stencil[6]
				- 1704396. * f_stencil[7]
				+ 364863.  * f_stencil[8])
			+ f_stencil[7] * (
				  482963. * f_stencil[7]
				- 208501. * f_stencil[8])
			+ 22658. * f_stencil[8] * f_stencil[8];

	return betas;
}


//template <ArithmeticWith<numeric_val> T>
////std::array<T, 3> smoothness_indicators(const T1& f_stencil) {
//std::valarray<T> betaSmoothnessIndicatorsWENO9BS(
//		const std::ranges::sized_range auto& f_stencil) {
//	/* Return the WENO9 smoothness indicators of Balsara and Shu (2000)
//	 * for each of the 5 substencils.
//	 * That is the sum of the normalized squares of the scaled
//	 * L2-norms of all the derivatives of 5 local interpolating
//	 * polynomials in the sub-stencils of 9-node `f_stencil`.
//	 *
//	 * This allows for (2*5-1)=9th order accuracy from the 5th
//	 * order ENO schemes.
//	 */

//	// std::array<T, 3> res;
//	std::valarray<T> betas(5);
//	static const Eigen::DiagonalMatrix<T, 5> l0(
//				std::sqrt(793030.)/420.,
//				std::sqrt(193716587562.)/543792.,
//				2.*std::sqrt(248759756350023.)/42747945.,
//				std::sqrt(96904093308502161.)/349153284.,
//				0.);
//	static const Eigen::Matrix<T, 5, 5> c0 {
//		{1., -208501./45316., 364863./45316., -288007./45316., 86329./45316.},
//		{0., 1., -55338513./14249315., 71911201./14249315., -1813059./838195.},
//		{0., 0., 1., -583582777./232768856., 350813921./232768856.},
//		{0., 0., 0., 1., -1.},
//		{0., 0., 0., 0., 0.}
//	};
//	static const Eigen::Matrix<T, 5, 5> b0 = l0 * c0;

//	static const Eigen::DiagonalMatrix<T, 5> l1(
//				std::sqrt(60445.)/210.,
//				std::sqrt(677061715.)/27632.,
//				2.*std::sqrt(11406983268815.)/5880675.,
//				std::sqrt(96904093308502161.)/349153284.,
//				0.);
//	static const Eigen::Matrix<T, 5, 5> c1 {
//		{1., -60871./13816., 99213./13816., -70237./13816., 18079./13816.},
//		{0., 1., -18329233./5880675., 67923./22025., -5686883./5880675.},
//		{0., 0., 1., -379530087./232768856., 146761231./232768856.},
//		{0., 0., 0., 1., -1.},
//		{0., 0., 0., 0., 0.}
//	};
//	static const Eigen::Matrix<T, 5, 5> b1 = l1 * c1;

//	/*constexpr*/ static Eigen::DiagonalMatrix<T, 5> l2(
//				std::sqrt(60445.)/210.,
//				std::sqrt(14765140203.)/82896.,
//				std::sqrt(12846941339118.)/2514585.,
//				std::sqrt(1446303987150483714.)/2605575108.,
//				0.);
//	static const Eigen::Matrix<T, 5, 5> c2 {
//		{1., -51001./13816., 67923./13816., -38947./13816., 8209./13816.},
//		{0., 1., -1870849./838195., 23242001./14249315., -5686883./14249315.},
//		{0., 0., 1., -71725821./51089708., 20636113./51089708.},
//		{0., 0., 0., 1., -1.},
//		{0., 0., 0., 0., 0.}
//	};
//	static const Eigen::Matrix<T, 5, 5> b2 = l2 * c2;

//	static const Eigen::DiagonalMatrix<T, 5> l3(
//				std::sqrt(793030.)/420.,
//				7.*std::sqrt(30758438922.)/543792.,
//				std::sqrt(2862368513038002.)/47512815.,
//				std::sqrt(7022472797218781694.)/12651268668.,
//				0.);
//	static const Eigen::Matrix<T, 5, 5> c3 {
//		{1., -140251./45316., 165153./45316., -88297./45316., 18079./45316.},
//		{0., 1., -217591953./110863235., 19650103./15837605., -30822003./110863235.},
//		{0., 0., 1., -5962732027./4217089556., 1745642471./4217089556.},
//		{0., 0., 0., 1., -1.},
//		{0., 0., 0., 0., 0.}
//	};
//	static const Eigen::Matrix<T, 5, 5> b3 = l3 * c3;

//	static const Eigen::DiagonalMatrix<T, 5> l4(
//				std::sqrt(3777130.)/420.,
//				std::sqrt(6405506236662.)/863344.,
//				2.*std::sqrt(132391829280434178.)/890329635.,
//				std::sqrt(7428632088185797566.)/26765962104.,
//				0.);
//	static const Eigen::Matrix<T, 5, 5> c4 {
//		{1., -649501./215836., 758823./215836., -411487./215836., 86329./215836.},
//		{0., 1., -1835635153./890329635., 411792427./296776545., -290071763./890329635.},
//		{0., 0., 1., -26228937057./17843974736., 8384962321./17843974736.},
//		{0., 0., 0., 1., -1.},
//		{0., 0., 0., 0., 0.}
//	};
//	static const Eigen::Matrix<T, 5, 5> b4 = l4 * c4;

//	Eigen::Matrix<T, 5, 1> substencil0 {
//		f_stencil[0],
//		f_stencil[1],
//		f_stencil[2],
//		f_stencil[3],
//		f_stencil[4]
//	};
//	Eigen::Matrix<T, 5, 1> substencil1 {
//		f_stencil[1],
//		f_stencil[2],
//		f_stencil[3],
//		f_stencil[4],
//		f_stencil[5]
//	};
//	Eigen::Matrix<T, 5, 1> substencil2 {
//		f_stencil[2],
//		f_stencil[3],
//		f_stencil[4],
//		f_stencil[5],
//		f_stencil[6]
//	};
//	Eigen::Matrix<T, 5, 1> substencil3 {
//		f_stencil[3],
//		f_stencil[4],
//		f_stencil[5],
//		f_stencil[6],
//		f_stencil[7]
//	};
//	Eigen::Matrix<T, 5, 1> substencil4 {
//		f_stencil[4],
//		f_stencil[5],
//		f_stencil[6],
//		f_stencil[7],
//		f_stencil[8]
//	};

//	betas[0] = ((b0 * substencil0).transpose())
//				* (b0 * substencil0);
//	betas[1] = ((b1 * substencil1).transpose())
//				* (b1 * substencil1);
//	betas[2] = ((b2 * substencil2).transpose())
//				* (b2 * substencil2);
//	betas[3] = ((b3 * substencil3).transpose())
//				* (b3 * substencil3);
//	betas[4] = ((b4 * substencil4).transpose())
//				* (b4 * substencil4);

//	return betas;
//}


//template <ArithmeticWith<numeric_val> T>
//std::valarray<T> betaSmoothnessIndicatorsMat(
//		const std::ranges::common_range auto& f_stencil,
//		std::array<std::array<const T, 6>, 6> _coefs[] = plus_coefs) {
//	/* Return the smoothness indicators beta_k, k=0,1,2
//	 * for each of the 3 substencils of `f_stencil`.
//	 */

//	std::valarray<T> res(3);

//	for (std::size_t k = 0; k < 3; ++ k) {
//		// for (std::size_t k = 0; k < half_size + 1; ++ k)
//		res[k] = std::inner_product(
//			std::ranges::begin(f_stencil),
//			std::ranges::end(f_stencil),
//			std::ranges::begin(vecMatDot(f_stencil, _coefs[k])),
//			0.
//		);
//	}


//	return res;
//}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> f3OrdReconstructionFromStencil(
		const std::ranges::sized_range auto& f_stencil) {
	/* 3rd order reconstructions of f(j) from all the 3 3-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *               (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')).
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 */
	//		{{2./6, -7./6, 11./6, 0., 0., 0.}},
	//		{{0., -1./6, 5./6, 2./6, 0., 0.}},
	//		{{0., 0., 2./6, 5./6, -1./6, 0.}}

	std::valarray<T> q_res(3);

	q_res[0] = (2. * f_stencil[0]
			  - 7. * f_stencil[1]
			 + 11. * f_stencil[2]) / 6.;

	q_res[1] = (-1. * f_stencil[1]
			   + 5. * f_stencil[2]
			   + 2. * f_stencil[3]) / 6.;

	q_res[2] = (2. * f_stencil[2]
			  + 5. * f_stencil[3]
			  - 1. * f_stencil[4]) / 6.;

	return q_res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> f4OrdReconstructionFromStencil(
		const std::ranges::sized_range auto& f_stencil) {
	/* 4th order reconstructions of f(j) from all the 4 4-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 7 values
	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')).
	 *       ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    |
	 */

	std::valarray<T> q_res(4);

	q_res[0] = ( - 3. * f_stencil[0]
				+ 13. * f_stencil[1]
				- 23. * f_stencil[2]
				+ 25. * f_stencil[3]) / 12.;

	q_res[1] = ( + 1. * f_stencil[1]
				 - 5. * f_stencil[2]
				+ 13. * f_stencil[3]
				 + 3. * f_stencil[4]) / 12.;

	q_res[2] = ( - 1. * f_stencil[2]
				 + 7. * f_stencil[3]
				 + 7. * f_stencil[4]
				 - 1. * f_stencil[5]) / 12.;

	q_res[3] = ( + 3. * f_stencil[3]
				+ 13. * f_stencil[4]
				 - 5. * f_stencil[5]
				 + 1. * f_stencil[6]) / 12.;

	return q_res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> f5OrdReconstructionFromStencil(
		const std::ranges::sized_range auto& f_stencil) {
	/* 5th order reconstructions of f(j) from all the 4 4-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 9 values
	 *     [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4, ...] for '+'
	 * (or [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...] for '-')).
	 *       ^    ^    ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    7    8    |
	 */

	std::valarray<T> q_res(5);

	q_res[0] = (+  12. * f_stencil[0]
				-  63. * f_stencil[1]
				+ 137. * f_stencil[2]
				- 163. * f_stencil[3]
				+ 137. * f_stencil[4]) / 60.;

	q_res[1] = (-  3. * f_stencil[1]
				+ 17. * f_stencil[2]
				- 43. * f_stencil[3]
				+ 77. * f_stencil[4]
				+ 12. * f_stencil[5]) / 60.;

	q_res[2] = (+  2. * f_stencil[2]
				- 13. * f_stencil[3]
				+ 47. * f_stencil[4]
				+ 27. * f_stencil[5]
				-  3. * f_stencil[6]) / 60.;

	q_res[3] = (-  3. * f_stencil[3]
				+ 27. * f_stencil[4]
				+ 47. * f_stencil[5]
				- 13. * f_stencil[6]
				+  2. * f_stencil[7]) / 60.;

	q_res[4] = (+ 12. * f_stencil[4]
				+ 77. * f_stencil[5]
				- 43. * f_stencil[6]
				+ 17. * f_stencil[7]
				-  3. * f_stencil[8]) / 60.;

	return q_res;
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> f6OrdReconstructionFromStencil(
		const std::ranges::sized_range auto& f_stencil) {
	/* 5th order reconstructions of f(j) from all the 4 4-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 11 values...).
	 */

	std::valarray<T> q_res(6);

	q_res[0] = (-  10. * f_stencil[0]
				+  62. * f_stencil[1]
				- 163. * f_stencil[2]
				+ 237. * f_stencil[3]
				- 213. * f_stencil[4]
				+ 147. * f_stencil[5]) / 60.;

	q_res[1] = (+  2. * f_stencil[1]
				- 13. * f_stencil[2]
				+ 37. * f_stencil[3]
				- 63. * f_stencil[4]
				+ 87. * f_stencil[5]
				+ 10. * f_stencil[6]) / 60.;

	q_res[2] = (-  1. * f_stencil[2]
				+  7. * f_stencil[3]
				- 23. * f_stencil[4]
				+ 57. * f_stencil[5]
				+ 22. * f_stencil[6]
				-  2. * f_stencil[7]) / 60.;

	q_res[3] = (+  1. * f_stencil[3]
				-  8. * f_stencil[4]
				+ 37. * f_stencil[5]
				+ 37. * f_stencil[6]
				-  8. * f_stencil[7]
				+  1. * f_stencil[8]) / 60.;

	q_res[4] = (-  2. * f_stencil[4]
				+ 22. * f_stencil[5]
				+ 57. * f_stencil[6]
				- 23. * f_stencil[7]
				+  7. * f_stencil[8]
				-  1. * f_stencil[9]) / 60.;

	q_res[5] = (+ 10. * f_stencil[5]
				+ 87. * f_stencil[6]
				- 63. * f_stencil[7]
				+ 37. * f_stencil[8]
				- 13. * f_stencil[9]
				+  2. * f_stencil[10]) / 60.;

	return q_res;
}


template <ArithmeticWith<numeric_val> T>
T henrickGMappingForLambda(T lambda_weno_weight,
						   T lambda_ideal = 1./3.) {
	/* The mapping function g by Henrick modified for symmetric
	 * lambda-weights by Zheng Hong, Zhengyin Ye and Kun Ye.
	 */

	T square_ideal = lambda_ideal * lambda_ideal;

	return lambda_weno_weight
			* (lambda_ideal
			   + square_ideal
			   - 3. * lambda_ideal * lambda_weno_weight
			   + lambda_weno_weight * lambda_weno_weight)
			/ (square_ideal
			   + lambda_weno_weight
			   * (1. - 2. * lambda_ideal));
}


//template <ArithmeticWith<numeric_val> T>
//T henrickGMapping(T omega_weno_weight,
//				  T d_ideal = 1./3.) {
//	/* The original mapping function g by Henrick et al. */

//	T d_square = d_ideal * d_ideal;
//	return omega_weno_weight
//			* (d_ideal
//			   + d_square
//			   - 3. * d_ideal * omega_weno_weight
//			   + omega_weno_weight * omega_weno_weight)
//			/ (d_square
//			   + omega_weno_weight
//			   * (1. - 2.*d_ideal));
//}


// const std::ranges::common_range auto&&
template <ArithmeticWith<numeric_val> T>
T alphaWENO5FMWeight(
		T beta_IS_coefficient,
		T epsilon = 1e-40,
		T p = 2.) {
	/* Compute appropriate alpha(α)-weights for the WENO5-FM scheme,
	 * by which the inverses of smoothness indicators are meant,
	 * so inverse beta(β) with the caveat of aritificially finite
	 * answers using the added epsilon-parameter to the beta-weights.
	 *
	 * `p` controls (increases) the amount of numerical dissipation
	 * (it's recommended to take it = r-1 for 2r-1 order schemes,
	 * so 2 for WENO5).
	 */

	return 1. / std::pow(epsilon + beta_IS_coefficient, p);
}


template <ArithmeticWith<numeric_val> T>
T alphaWENO5ZMWeight(
		T beta_IS_coefficient,
		T tau_5,
		T epsilon = 1e-40,
		T p = 2.) {
	/* Compute appropriate alpha(α)-weights for the WENO5-ZM scheme,
	 * by which the inverses of smoothness indicators are meant,
	 * so inverse beta(β) with the caveat of aritificially finite
	 * answers using the added epsilon-parameter to the beta-weights.
	 *
	 * `p` controls (increases) the amount of numerical dissipation
	 * (it's recommended to take it = r-1 for 2r-1 order schemes,
	 * so 2 for WENO5).
	 */

	return 1. + std::pow(tau_5 / (beta_IS_coefficient + epsilon), p);
}


template <ArithmeticWith<numeric_val> T>
void lambdaWENO5FMWeights(
		const std::ranges::common_range auto&& alpha_weights,
		std::ranges::common_range auto&& res) {
	/* FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
	 * for WENO5-FM or WENO5-ZM
	 * due to Zheng Hong, Zhengyin Ye and Kun Ye:
	 * lambda_weights = alpha_weights / alpha_weights.sum();
	 */

	// return alpha_weights / alpha_weights.sum();
	T sum = std::reduce(
				/*std::execution::par_unseq,*/
				std::ranges::begin(alpha_weights),
				std::ranges::end(alpha_weights),
				0.);

	std::transform(
				/*std::execution::par_unseq,*/
				std::ranges::begin(alpha_weights),
				std::ranges::end(alpha_weights),
				std::ranges::begin(res),
				[sum](const auto alpha) {
		return alpha / sum;
	});
}


template <ArithmeticWith<numeric_val> T>
std::valarray<T> prediscretizeWENO5LambdaMapping(
		std::size_t N,
		T lambda_ideal = 1./3.) {
	/* Pre-discrete mapping method for WENO5-FM
	 * (idea due to Hong et al.) to increase performace.
	 * Construct the lambda mapping `valarray` of size `N + 1`
	 * for lambdas in the range [0..1].
	 */

	std::valarray<T> res_lookup_table(N + 1);

//	auto ns = std::ranges::common_view(
//			std::ranges::views::iota(std::size_t(0))
//				| std::views::take(N + 1)
//	);

//	std::transform(
//				std::execution::par_unseq,
//				std::ranges::begin(ns),
//				std::ranges::end(ns),
//				std::ranges::begin(res_lookup_table),
//				[N](std::size_t n) {
//		return henrickGMappingForLambda<T>(
//					static_cast<T>(n) / static_cast<T>(N));
//	});

	for (std::size_t n = 0; n <= N; ++ n)
		res_lookup_table[n] = henrickGMappingForLambda<T>(
					static_cast<T>(n) / static_cast<T>(N),
					lambda_ideal);

	return res_lookup_table;
}


std::valarray<numeric_val> DISCRETE_LAMBDA5
			= prediscretizeWENO5LambdaMapping<numeric_val>(1000000/*00*/,
														   1./3.);

std::valarray<numeric_val> DISCRETE_LAMBDA7
			= prediscretizeWENO5LambdaMapping<numeric_val>(100000000,
														   1./4.);

std::valarray<numeric_val> DISCRETE_LAMBDA9
			= prediscretizeWENO5LambdaMapping<numeric_val>(100000000,
														   1./5.);

std::valarray<numeric_val> DISCRETE_LAMBDA11
			= prediscretizeWENO5LambdaMapping<numeric_val>(100000000,
														   1./6.);


template <ArithmeticWith<numeric_val> T>
std::ranges::common_range auto omegaWENOFMWeights(
		const std::ranges::common_range auto&& lambda_weights,
		const std::valarray<numeric_val>& d_ideal_lin_weights,
		const std::valarray<numeric_val>& discrete_lambda
					/*= DISCRETE_LAMBDA5*/) {
	/* From Henrick et al.'s mappings of g(λ_k) for the improved
	 * symmetric normalized lambda-weights of Hong, Ye & Ye
	 * and linear weights d_k we get the new corrected resultant
	 * normalized WENO5-FM (WENO5-ZM) omega (ω_k-)weights for WENO5-FM
	 * (again due to Zheng Hong, Zhengyin Ye and Kun Ye).
	 */

	// From them WENO5-Z and WENO-M will calculate the non-linear
	// alpha and omega weights.

	// In WENO5-FM, further, we have one ideal value for λ
	// \overbar{lambda} = 1/3
	// T lambda_ideal = 1/3;
	// In the smooth region the smoothness indicators β_k ought to
	// be equal for all sub-stencils, and thus the weight's ideal
	// value must be unique.

	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
	// omega_weights = d_lin_weights * alpha_weights;
	// lambda_weights

	// And only to the λ-weights a mapping in the spirit of
	// Henrick et al. is applied:
//	auto gMap = [](T x) -> T {
//		return henrickGMappingForLambda(x);
//	};

//	std::valarray<T> alpha_weights(3);
//	std::ranges::transform(
//				lambda_weights,
//				std::ranges::begin(alpha_weights),
//				gMap);
	const std::size_t n = std::ranges::size(discrete_lambda) - 1;
	std::valarray<T> alpha_weights(std::ranges::size(d_ideal_lin_weights));
	std::ranges::transform(
				lambda_weights,
				std::ranges::begin(alpha_weights),
				[n, &discrete_lambda](T lambda) -> T {
		return discrete_lambda[
				static_cast<std::size_t>(
					static_cast<T>(n) * lambda)];
	});
	// α*-weights

	// From α*=g(λ_k) and d_k we get the new corrected resultant
	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
	// omega_weights = d_lin_weights * alpha_weights;
	std::valarray<T> omega_weights = d_ideal_lin_weights
			* alpha_weights;
	omega_weights /= omega_weights.sum();

	return omega_weights;
}


template <ArithmeticWith<numeric_val> T>
std::ranges::common_range auto omegaWENO5FMWeights(
		const std::ranges::common_range auto&& lambda_weights) {
	// The ideal weights (they generate the central upstream fifth-order
	// scheme for the 5-point stencil), which are in WENO usu. called
	// (optimal) linear weights:
	std::valarray<T> d_lin_weights = {0.1, 0.6, 0.3};

	return omegaWENOFMWeights<T>(
				std::move(lambda_weights), d_lin_weights,
				DISCRETE_LAMBDA5);
}


template <ArithmeticWith<numeric_val> T>
std::ranges::common_range auto omegaWENO7FMWeights(
		const std::ranges::common_range auto&& lambda_weights) {
	// The ideal weights (they generate the central upstream seventh-order
	// scheme for the 7-point stencil), which are in WENO usu. called
	// (optimal) linear weights:
	// std::valarray<T> d_lin_weights = {4./35., 18./35., 12./35., 1./35.};
	std::valarray<T> d_lin_weights = {1./35., 12./35., 18./35., 4./35.};

	return omegaWENOFMWeights<T>(
				std::move(lambda_weights), d_lin_weights,
				DISCRETE_LAMBDA7);
}


template <ArithmeticWith<numeric_val> T>
std::ranges::common_range auto omegaWENO9FMWeights(
		const std::ranges::common_range auto&& lambda_weights) {
	// The ideal weights (they generate the central upstream ninth-order
	// scheme for the 9-point stencil), which are in WENO usu. called
	// (optimal) linear weights:
	std::valarray<T> d_lin_weights = {
		1./126., 10./63., 10./21., 20./63., 5./126.
	};

	return omegaWENOFMWeights<T>(
				std::move(lambda_weights), d_lin_weights,
				DISCRETE_LAMBDA9);
}


template <ArithmeticWith<numeric_val> T>
std::ranges::common_range auto omegaWENO11FMWeights(
		const std::ranges::common_range auto&& lambda_weights) {
	// The ideal weights (they generate the central upstream ninth-order
	// scheme for the 9-point stencil), which are in WENO usu. called
	// (optimal) linear weights:
	std::valarray<T> d_lin_weights = {
		1./462., 30./462., 150./462., 200./462., 75./462., 6./462.
	};

	return omegaWENOFMWeights<T>(
				std::move(lambda_weights), d_lin_weights,
				DISCRETE_LAMBDA11);
}


//template <ArithmeticWith<numeric_val> T>
//T postINDEX(T omega_weight1, T omega_weight2, auto&& gMap) {
//	return (omega_weight1 - omega_weight2)
//			* (gMap(omega_weight1) - gMap(omega_weight2));
//}


//template <ArithmeticWith<numeric_val> T>
//bool isOPPoint(T omega_weight1, T omega_weight2, auto&& gMap) {
//	const T ind = postINDEX<T>(omega_weight1, omega_weight2, gMap);
//	const T zero = static_cast<T>(0.);

//	return ind > zero
//			|| ((omega_weight1 - omega_weight2) == zero && ind == zero);
//}


//template <ArithmeticWith<numeric_val> T>
//void applyLOPWeightMap(
//		const std::ranges::common_range auto&& lambda_weights,
//		auto&& gMap) {
//	const std::size_t len = std::ranges::size(lambda_weights);
//	bool is_op = true;
//	for (std::size_t j = 0; j < len - 1 && !is_op; ++ j)
//		for (std::size_t k = j; k < len && !is_op; ++ k)
//			is_op = isOPPoint<T>(
//						lambda_weights[j], lambda_weights[k], gMap);

//	if (is_op)
//		for (std::size_t j = 0; j < len && !is_op; ++ j)
//			lambda_weights[j] = gMap(lambda_weights[j]);
//}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5JSReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);
//	beta_IS_coefs[0] = (13.0/12.0)*std::pow(
//				f_stencil[0]-2.0*f_stencil[1]+f_stencil[2], 2)
//			+ 0.25*std::pow(
//				f_stencil[0]-4.0*f_stencil[1]+3.0*f_stencil[2], 2);
//	beta_IS_coefs[1] = (13.0/12.0)*std::pow(
//				f_stencil[1]-2.0*f_stencil[2]+f_stencil[3], 2)
//			+ 0.25*std::pow(f_stencil[1]-f_stencil[3], 2);
//	beta_IS_coefs[2] = (13.0/12.0)*std::pow(
//				f_stencil[2]-2.0*f_stencil[3]+f_stencil[4], 2)
//			+ 0.25*std::pow(
//				3.0*f_stencil[2]-4.0*f_stencil[3]+f_stencil[4], 2);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {0.1, 0.6, 0.3};
	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
//	std::valarray<T> alpha_weights(3);
//	alpha_weights[0] = 1.0e-1/std::pow((eps+beta_IS_coefs[0]), 2);
//	alpha_weights[1] = 6.0e-1/std::pow((eps+beta_IS_coefs[1]), 2);
//	alpha_weights[2] = 3.0e-1/std::pow((eps+beta_IS_coefs[2]), 2);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();
//	std::valarray<T> omega_weights(3);
//	omega_weights[0] = alpha_weights[0]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);
//	omega_weights[1] = alpha_weights[1]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);
//	omega_weights[2] = alpha_weights[2]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<T> eno_reconstructed_f
			= f3OrdReconstructionFromStencil<T>(f_stencil);

//	 eno_reconstructed_f[0] = f_stencil[0]/3.0
//			 - 7.0/6.0*f_stencil[1] + 11.0/6.0*f_stencil[2];
//	 eno_reconstructed_f[1] =-f_stencil[1]/6.0
//			 + 5.0/6.0*f_stencil[2] + f_stencil[3]/3.0;
//	 eno_reconstructed_f[2] = f_stencil[2]/3.0
//			 + 5.0/6.0*f_stencil[3] - f_stencil[4]/6.0;

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5JSReconstructionKernelRev(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {0.3, 0.6, 0.1};
	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights / alpha_weights.sum();

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of f_{i+1/2}
	// via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<T> eno_reconstructed_f(3);
	eno_reconstructed_f[0] = (-1. * f_stencil[0]
						   + 5. * f_stencil[1]
						   + 2. * f_stencil[2]) / 6.;
	eno_reconstructed_f[1] = (2. * f_stencil[1]
						   + 5. * f_stencil[2]
						   - 1. * f_stencil[3]) / 6.;
	eno_reconstructed_f[2] = (11. * f_stencil[2]
						   - 7. * f_stencil[3]
						   + 2. * f_stencil[4]) / 6.;

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];


	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5MReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience). Use the WENO-HAP (WENO5-M)
	 * procedure.
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
//	std::valarray<numeric_val> DISCRETE_LAMBDA5
//				= prediscretizeWENO5LambdaMapping<numeric_val>(1000000/*00*/,
//															   1./3.);
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {0.1, 0.6, 0.3};
	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::transform(
				/*std::execution::par_unseq,*/
				std::ranges::begin(omega_weights),
				std::ranges::end(omega_weights),
				std::ranges::begin(d_lin_weights),
				std::ranges::begin(omega_weights),
				[](auto w, auto d) {
		return henrickGMappingForLambda(w, d);
	});  // we obtain a new non-normalized weight alpha*

	omega_weights = omega_weights
			/ omega_weights.sum();  // normalize it

	std::valarray<T> eno_reconstructed_f
			= f3OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


//template <ArithmeticWith<numeric_val> T>
//T computeFHatWENO5MReconstructionKernel(
//		const std::ranges::sized_range auto&& f_stencil,
//		T eps = 1e-40, T p = 2.) {
//	/* Calculate (reconstruct) one of the two split monotone numerical
//	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
//	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
//	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
//	 *                      ^    ^    ^    ^    ^    ^
//	 *                      0    1    2    3    4    |
//	 * in either case for convenience).
//	 *
//	 * I.e. this function implements the upwind reconstruction which
//	 * should be used for positive fluxes (with information propagated
//	 * from left to right) if the nodes are passed in order. However,
//	 * the downwind reconstruction should obviously look the same
//	 * modulo flipping the values with respect to j+0, so that it
//	 * becomes downwind biased and takes one extra point to the right
//	 * instead of taking one extra to the left. In other words, to get
//	 * this to behave as a downwind reconstrution we need to pass
//	 * the points symmetric to those of upwind reconstruction with
//	 * respect to j+0:
//	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
//	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
//	 *                  |
//	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
//	 */

//	T beta0 = (13./12.)
//			* (f_stencil[0] - 2. * f_stencil[1] + f_stencil[2])
//			* (f_stencil[0] - 2. * f_stencil[1] + f_stencil[2])
//			+ (1./4.)
//			* (f_stencil[0] - 4. * f_stencil[1] + 3. * f_stencil[2])
//			* (f_stencil[0] - 4. * f_stencil[1] + 3. * f_stencil[2]);

//	T beta1 = (13./12.)
//			* (f_stencil[1] - 2. * f_stencil[2] + f_stencil[3])
//			* (f_stencil[1] - 2. * f_stencil[2] + f_stencil[3])
//			+ (1./4.)
//			* (f_stencil[3] - f_stencil[1])
//			* (f_stencil[3] - f_stencil[1]);

//	T beta2 = (13./12.)
//			* (f_stencil[2] - 2. * f_stencil[3] + f_stencil[4])
//			* (f_stencil[2] - 2. * f_stencil[3] + f_stencil[4])
//			+ (1./4.)
//			* (3. * f_stencil[2] - 4. * f_stencil[3] + f_stencil[4])
//			* (3. * f_stencil[2] - 4. * f_stencil[3] + f_stencil[4]);

//	T f_hat = 0.;

////	T alpha0 = 	(1./10.) / std::pow(eps + beta0, p);
////	T alpha1 = 	(6./10.) / std::pow(eps + beta1, p);
////	T alpha2 = 	(3./10.) / std::pow(eps + beta2, p);
//	T alpha0 = 	(1./10.) / (1e-40 + beta0) / (1e-40 + beta0);
//	T alpha1 = 	(6./10.) / (1e-40 + beta1) / (1e-40 + beta1);
//	T alpha2 = 	(3./10.) / (1e-40 + beta2) / (1e-40 + beta2);

//	T alpha_sum = alpha0 + alpha1 + alpha2;

//	T omega0_ast = alpha0 / alpha_sum;
//	T omega1_ast = alpha1 / alpha_sum;
//	T omega2_ast = alpha2 / alpha_sum;

//	T alpha0_g_ast = omega0_ast
//			* ((1./10.)
//			   + (1./10.) * (1./10.)
//			   - 3. * (1./10.) * omega0_ast
//			   + omega0_ast * omega0_ast)
//			/ ((1./10.) * (1./10.)
//			   + (1. - 2. * (1./10.)) * omega0_ast);

//	T alpha1_g_ast = omega1_ast
//			* ((6./10.)
//			   + (6./10.) * (6./10.)
//			   - 3. * (6./10.) * omega1_ast
//			   + omega1_ast * omega1_ast)
//			/ ((6./10.) * (6./10.)
//			   + (1. - 2. * (6./10.)) * omega1_ast);

//	T alpha2_g_ast = omega2_ast
//			* ((3./10.)
//			   + (3./10.) * (3./10.)
//			   - 3. * (3./10.) * omega2_ast
//			   + omega2_ast * omega2_ast)
//			/ ((3./10.) * (3./10.)
//			   + (1. - 2. * (3./10.)) * omega2_ast);
////	T alpha0_g_ast = henrickGMappingForLambda<T>(omega0_ast, 0.1);
////	T alpha1_g_ast = henrickGMappingForLambda<T>(omega1_ast, 0.6);
////	T alpha2_g_ast = henrickGMappingForLambda<T>(omega2_ast, 0.3);

//	T alpha_ast_sum = alpha0_g_ast + alpha1_g_ast + alpha2_g_ast;
//	T omega0 = alpha0_g_ast / alpha_ast_sum;
//	T omega1 = alpha1_g_ast / alpha_ast_sum;
//	T omega2 = alpha2_g_ast / alpha_ast_sum;

//	T q0 = (1./6.) * (2. * f_stencil[0]
//			- 7. * f_stencil[1]
//			+ 11. * f_stencil[2]);
//	T q1 = (1./6.) * (-1. * f_stencil[1]
//			+ 5. * f_stencil[2]
//			+ 2. * f_stencil[3]);
//	T q2 = (1./6.) * (2. * f_stencil[2]
//			+ 5. * f_stencil[3]
//			- 1. * f_stencil[4]);

//	f_hat = omega0 * q0
//			+ omega1 * q1
//			+ omega2 * q2;

//	return f_hat;
//}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5FMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);

	std::array<T, 3> alpha_weights;
	std::ranges::transform(
				beta_IS_coefs,
				std::ranges::begin(alpha_weights),
				[eps, p](auto beta) {
		return alphaWENO5FMWeight(beta, eps, p);
	});

	std::array<T, 3> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	std::valarray<T> omega_weights = omegaWENO5FMWeights<T>(
				std::move(lambda_weights));

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<T> eno_reconstructed_f
			= f3OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO7BSReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 7 values
	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...]. (We reverse the points in
	 *      [j-3, j-2, j-1, j+0, j+1, j+2, j+3] j+4 and get
	 *                       |
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2] j-3.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	std::valarray<T> beta_IS_coefs(4);

	T f_hat = 0.;

	beta_IS_coefs = betaSmoothnessIndicatorsWENO7BS<T>(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {1./35., 12./35., 18./35., 4./35.};

	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	std::valarray<T> eno_reconstructed_f
			= f4OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO7MReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 7 values
	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...]. (We reverse the points in
	 *      [j-3, j-2, j-1, j+0, j+1, j+2, j+3] j+4 and get
	 *                       |
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2] j-3.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(4);

	T f_hat = 0.;
	beta_IS_coefs = betaSmoothnessIndicatorsWENO7BS<T>(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {1./35., 12./35., 18./35., 4./35.};
	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	std::transform(
				/*std::execution::par_unseq,*/
				std::ranges::begin(omega_weights),
				std::ranges::end(omega_weights),
				std::ranges::begin(d_lin_weights),
				std::ranges::begin(omega_weights),
				[](auto w, auto d) {
		return henrickGMappingForLambda(w, d);
	});  // we obtain a new non-normalized weight alpha*

	omega_weights = omega_weights
			/ omega_weights.sum();  // normalize it

	std::valarray<T> eno_reconstructed_f
			= f4OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO7FMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 7 values
	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...]. (We reverse the points in
	 *      [j-3, j-2, j-1, j+0, j+1, j+2, j+3] j+4 and get
	 *                       |
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2] j-3.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(4);

	T f_hat = 0.;
	beta_IS_coefs = betaSmoothnessIndicatorsWENO7BS<T>(f_stencil);

	std::array<T, 4> alpha_weights;
	std::ranges::transform(
				beta_IS_coefs,
				std::ranges::begin(alpha_weights),
				[eps, p](auto beta) {
		return alphaWENO5FMWeight(beta, eps, p);
	});

	std::array<T, 4> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	std::valarray<T> omega_weights = omegaWENO7FMWeights<T>(
				std::move(lambda_weights));

	std::valarray<T> eno_reconstructed_f
			= f4OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO7ZMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 7 values
	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...]. (We reverse the points in
	 *      [j-3, j-2, j-1, j+0, j+1, j+2, j+3] j+4 and get
	 *                       |
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2] j-3.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(4);

	T f_hat = 0.;
	beta_IS_coefs = betaSmoothnessIndicatorsWENO7BS<T>(f_stencil);

	T tau_7 = std::abs(-beta_IS_coefs[0]
				- 3. * beta_IS_coefs[1]
				+ 3. * beta_IS_coefs[2]
				+ beta_IS_coefs[3]);

	std::array<T, 4> alpha_weights;
	std::ranges::transform(
				beta_IS_coefs,
				std::ranges::begin(alpha_weights),
				[tau_7, eps, p](auto beta) {
		return alphaWENO5ZMWeight(beta, tau_7, eps, p);
	});

	std::array<T, 4> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	std::valarray<T> omega_weights = omegaWENO7FMWeights<T>(
				std::move(lambda_weights));

	std::valarray<T> eno_reconstructed_f
			= f4OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO7ZReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 7 values
	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...]. (We reverse the points in
	 *      [j-3, j-2, j-1, j+0, j+1, j+2, j+3] j+4 and get
	 *                       |
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2] j-3.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	std::valarray<T> beta_IS_coefs(4);

	T f_hat = 0.;

	beta_IS_coefs = betaSmoothnessIndicatorsWENO7BS<T>(f_stencil);

	T tau_7 = std::abs(-beta_IS_coefs[0]
				- 3. * beta_IS_coefs[1]
				+ 3. * beta_IS_coefs[2]
				+ beta_IS_coefs[3]);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {1./35., 12./35., 18./35., 4./35.};

	std::valarray<T> alpha_weights = d_lin_weights
			* (1. + std::pow(tau_7 / (eps + beta_IS_coefs), p));

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	std::valarray<T> eno_reconstructed_f
			= f4OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO7SReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-80, T p = 1.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 7 values
	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...]. (We reverse the points in
	 *      [j-3, j-2, j-1, j+0, j+1, j+2, j+3] j+4 and get
	 *                       |
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2] j-3.)
	 */

	std::valarray<T> beta_IS_coefs(4);
	std::valarray<T> c_k_coefs(4);

	T f_hat = 0.;
	auto f_stencil_view = std::ranges::views::all(f_stencil);

	for (std::size_t k = 0; k <= std::ranges::size(f_stencil) - 4; ++ k) {
		auto substencil = f_stencil_view
				| std::ranges::views::drop(k)
				| std::ranges::views::take(4);

		c_k_coefs[k] = -substencil[0]
					+ 3. * substencil[1]
					- 3. * substencil[2]
					+ substencil[3];

		beta_IS_coefs[k] = std::pow(
						substencil[0]
						- substencil[1]
						- substencil[2]
						+ substencil[3], 2)
					+ std::abs((
						-substencil[0]
							- substencil[1]
							+ substencil[2]
							+ substencil[3])
						* c_k_coefs[k]
					);
	}

	T tau_7s = std::pow(
					c_k_coefs[0]
					- c_k_coefs[1]
					- c_k_coefs[2]
					+ c_k_coefs[3], 2)
				+ std::abs(
					(-c_k_coefs[0]
						- c_k_coefs[1]
						+ c_k_coefs[2]
						+ c_k_coefs[3])
					* (-c_k_coefs[0]
						+ 3. * c_k_coefs[1]
						- 3. * c_k_coefs[2]
						+ c_k_coefs[3]));

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {1./35., 12./35., 18./35., 4./35.};

	std::valarray<T> alpha_weights = d_lin_weights
			* (1. + std::pow(tau_7s / (eps + beta_IS_coefs), p));

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	std::valarray<T> eno_reconstructed_f
			= f4OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO7SMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-80, T p = 1.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 7 values
	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...]. (We reverse the points in
	 *      [j-3, j-2, j-1, j+0, j+1, j+2, j+3] j+4 and get
	 *                       |
	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2] j-3.)
	 */

	std::valarray<T> beta_IS_coefs(4);
	std::valarray<T> c_k_coefs(4);

	T f_hat = 0.;
	auto f_stencil_view = std::ranges::views::all(f_stencil);

	for (std::size_t k = 0; k <= std::ranges::size(f_stencil) - 4; ++ k) {
		auto substencil = f_stencil_view
				| std::ranges::views::drop(k)
				| std::ranges::views::take(4);

		c_k_coefs[k] = -substencil[0]
					+ 3. * substencil[1]
					- 3. * substencil[2]
					+ substencil[3];

		beta_IS_coefs[k] = std::pow(
						substencil[0]
						- substencil[1]
						- substencil[2]
						+ substencil[3], 2)
					+ std::abs((
						-substencil[0]
							- substencil[1]
							+ substencil[2]
							+ substencil[3])
						* c_k_coefs[k]
					);
	}

	T tau_7s = std::pow(
					c_k_coefs[0]
					- c_k_coefs[1]
					- c_k_coefs[2]
					+ c_k_coefs[3], 2)
				+ std::abs(
					(-c_k_coefs[0]
						- c_k_coefs[1]
						+ c_k_coefs[2]
						+ c_k_coefs[3])
					* (-c_k_coefs[0]
						+ 3. * c_k_coefs[1]
						- 3. * c_k_coefs[2]
						+ c_k_coefs[3]));

	// non-linear non-scaled (α-)weights
	std::valarray<T> alpha_weights = 1.
				+ std::pow(tau_7s / (eps + beta_IS_coefs), p);

	std::array<T, 4> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = omegaWENO7FMWeights<T>(
				std::move(lambda_weights));

	std::valarray<T> eno_reconstructed_f
			= f4OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO9BSReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 9 values
	 *     [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4, ...] for '+'
	 * (or [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    7    8    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...].
	 * (We reverse the points in
	 *      [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4] j+5 and get
	 *                            |
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3] j-4.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(5);

	T f_hat = 0.;

	beta_IS_coefs = betaSmoothnessIndicatorsWENO9BS<T>(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {
		1./126., 10./63., 10./21., 20./63., 5./126.
	};

	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	std::valarray<T> eno_reconstructed_f
			= f5OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3]
			+ omega_weights[4] * eno_reconstructed_f[4];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO9MReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 9 values
	 *     [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4, ...] for '+'
	 * (or [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    7    8    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...].
	 * (We reverse the points in
	 *      [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4] j+5 and get
	 *                            |
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3] j-4.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(5);

	T f_hat = 0.;

	beta_IS_coefs = betaSmoothnessIndicatorsWENO9BS<T>(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {
		1./126., 10./63., 10./21., 20./63., 5./126.
	};
	std::valarray<T> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	std::transform(
				/*std::execution::par_unseq,*/
				std::ranges::begin(omega_weights),
				std::ranges::end(omega_weights),
				std::ranges::begin(d_lin_weights),
				std::ranges::begin(omega_weights),
				[](auto w, auto d) {
		return henrickGMappingForLambda(w, d);
	});  // we obtain a new non-normalized weight alpha*

	omega_weights = omega_weights
			/ omega_weights.sum();  // normalize it

	std::valarray<T> eno_reconstructed_f
			= f5OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3]
			+ omega_weights[4] * eno_reconstructed_f[4];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO9FMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 9 values
	 *     [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4, ...] for '+'
	 * (or [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    7    8    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...].
	 * (We reverse the points in
	 *      [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4] j+5 and get
	 *                            |
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3] j-4.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)
	constexpr const std::size_t substencil_order = 5;

	std::valarray<T> beta_IS_coefs(substencil_order);

	T f_hat = 0.;
	beta_IS_coefs = betaSmoothnessIndicatorsWENO9BS<T>(f_stencil);

	std::array<T, substencil_order> alpha_weights;
	std::ranges::transform(
				beta_IS_coefs,
				std::ranges::begin(alpha_weights),
				[eps, p](auto beta) {
		return alphaWENO5FMWeight(beta, eps, p);
	});

	std::array<T, substencil_order> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	std::valarray<T> omega_weights = omegaWENO9FMWeights<T>(
				std::move(lambda_weights));

	std::valarray<T> eno_reconstructed_f
			= f5OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3]
			+ omega_weights[4] * eno_reconstructed_f[4];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO9SReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-80, T p = 1.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 9 values
	 *     [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4, ...] for '+'
	 * (or [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    7    8    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...].
	 * (We reverse the points in
	 *      [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4] j+5 and get
	 *                            |
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3] j-4.)
	 */

	constexpr const std::size_t substencil_order = 5;

	std::valarray<T> beta_IS_coefs(substencil_order);
	std::valarray<T> c_k_coefs(substencil_order);

	T f_hat = 0.;
	auto f_stencil_view = std::ranges::views::all(f_stencil);

	auto a_func = [](auto substencilv) -> T {
		return substencilv[0]
					- 2. * substencilv[2]
					+ substencilv[4];
	};

	auto b_func = [](auto substencilv) -> T {
		return -substencilv[0]
					+ 2. * substencilv[1]
					- 2. * substencilv[3]
					+ substencilv[4];
	};

	auto c_func = [](auto substencilv) -> T {
		return substencilv[0]
					- 4. * substencilv[1]
					+ 6. * substencilv[2]
					- 4. * substencilv[3]
					+ substencilv[4];
	};

	T a;
	T b;
	for (std::size_t k = 0;
			k <= std::ranges::size(f_stencil) - substencil_order; ++ k) {
		auto substencil = f_stencil_view
				| std::ranges::views::drop(k)
				| std::ranges::views::take(substencil_order);

		a = a_func(substencil);
		b = b_func(substencil);
		c_k_coefs[k] = c_func(substencil);

		beta_IS_coefs[k] = std::pow(b, 2) + std::abs(a * c_k_coefs[k]);
	}

	T tau_s = std::pow(b_func(c_k_coefs), 2)
				+ std::abs(a_func(c_k_coefs) * c_func(c_k_coefs));

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {
		1./126., 10./63., 10./21., 20./63., 5./126.
	};

	std::valarray<T> alpha_weights = d_lin_weights
			* (1. + std::pow(tau_s / (eps + beta_IS_coefs), p));

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	std::valarray<T> eno_reconstructed_f
			= f5OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3]
			+ omega_weights[4] * eno_reconstructed_f[4];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO9SMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-80, T p = 1.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 9 values
	 *     [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4, ...] for '+'
	 * (or [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    7    8    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...].
	 * (We reverse the points in
	 *      [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4] j+5 and get
	 *                            |
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3] j-4.)
	 */

	constexpr const std::size_t substencil_order = 5;

	std::valarray<T> beta_IS_coefs(substencil_order);
	std::valarray<T> c_k_coefs(substencil_order);

	T f_hat = 0.;
	auto f_stencil_view = std::ranges::views::all(f_stencil);

	auto a_func = [](auto substencilv) -> T {
		return substencilv[0]
					- 2. * substencilv[2]
					+ substencilv[4];
	};

	auto b_func = [](auto substencilv) -> T {
		return -substencilv[0]
					+ 2. * substencilv[1]
					- 2. * substencilv[3]
					+ substencilv[4];
	};

	auto c_func = [](auto substencilv) -> T {
		return substencilv[0]
					- 4. * substencilv[1]
					+ 6. * substencilv[2]
					- 4. * substencilv[3]
					+ substencilv[4];
	};

	T a;
	T b;

	for (std::size_t k = 0;
			k <= std::ranges::size(f_stencil) - substencil_order; ++ k) {
		auto substencil = f_stencil_view
				| std::ranges::views::drop(k)
				| std::ranges::views::take(substencil_order);

		a = a_func(substencil);
		b = b_func(substencil);
		c_k_coefs[k] = c_func(substencil);

		beta_IS_coefs[k] = std::pow(b, 2) + std::abs(a * c_k_coefs[k]);
	}

	T tau_s = std::pow(b_func(c_k_coefs), 2)
				+ std::abs(a_func(c_k_coefs) * c_func(c_k_coefs));

	// non-linear non-scaled (α-)weights
	std::valarray<T> alpha_weights = 1.
				+ std::pow(tau_s / (eps + beta_IS_coefs), p);

	std::array<T, substencil_order> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = omegaWENO9FMWeights<T>(
				std::move(lambda_weights));

	std::valarray<T> eno_reconstructed_f
			= f5OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3]
			+ omega_weights[4] * eno_reconstructed_f[4];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO11SReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-80, T p = 1.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 11 values
	 *     [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4, ...] for '+'
	 * (or [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    7    8    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...].
	 * (We reverse the points in
	 *      [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4] j+5 and get
	 *                            |
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3] j-4.)
	 */

	constexpr const std::size_t substencil_order = 6;

	std::valarray<T> beta_IS_coefs(substencil_order);
	std::valarray<T> c_k_coefs(substencil_order);

	T f_hat = 0.;
	auto f_stencil_view = std::ranges::views::all(f_stencil);

	auto a_func = [](auto substencilv) -> T {
		return -substencilv[0]
					+ substencilv[1]
					+ 2. * substencilv[2]
					- 2. * substencilv[3]
					- substencilv[4]
					+ substencilv[5];
	};

	auto b_func = [](auto substencilv) -> T {
		return substencilv[0]
					- 3. * substencilv[1]
					+ 2. * substencilv[2]
					+ 2. * substencilv[3]
					- 3. * substencilv[4]
					+ substencilv[5];
	};

	auto c_func = [](auto substencilv) -> T {
		return -substencilv[0]
					+ 5. * substencilv[1]
					- 10. * substencilv[2]
					+ 10. * substencilv[3]
					- 5. * substencilv[4]
					+ substencilv[5];
	};

	T a;
	T b;
	for (std::size_t k = 0;
			k <= std::ranges::size(f_stencil) - substencil_order; ++ k) {
		auto substencil = f_stencil_view
				| std::ranges::views::drop(k)
				| std::ranges::views::take(substencil_order);

		a = a_func(substencil);
		b = b_func(substencil);
		c_k_coefs[k] = c_func(substencil);

		beta_IS_coefs[k] = std::pow(b, 2) + std::abs(a * c_k_coefs[k]);
	}

	T tau_s = std::pow(b_func(c_k_coefs), 2)
				+ std::abs(a_func(c_k_coefs) * c_func(c_k_coefs));

	// non-linear non-scaled (α-)weights
	std::valarray<T> d_lin_weights = {
		1./462., 30./462., 150./462., 200./462., 75./462., 6./462.
	};
//	std::valarray<T> d_lin_weights = {
//		1./462., 5./77., 25./77., 100./231., 25./154., 1./77.
//	};

	std::valarray<T> alpha_weights = d_lin_weights
			* (1. + std::pow(tau_s / (eps + beta_IS_coefs), p));

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = alpha_weights
			/ alpha_weights.sum();

	std::valarray<T> eno_reconstructed_f
			= f6OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3]
			+ omega_weights[4] * eno_reconstructed_f[4]
			+ omega_weights[5] * eno_reconstructed_f[5];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO11SMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-80, T p = 1.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives the following 11 values
	 *     [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4, ...] for '+'
	 * (or [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...] for '-')
	 *       ^    ^    ^    ^    ^    ^    ^    ^    ^    ^
	 *       0    1    2    3    4    5    6    7    8    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3, ...].
	 * (We reverse the points in
	 *      [j-4, j-3, j-2, j-1, j+0, j+1, j+2, j+3, j+4] j+5 and get
	 *                            |
	 * [j+5, j+4, j+3, j+2, j+1, j+0, j-1, j-2, j-3] j-4.)
	 */

	constexpr const std::size_t substencil_order = 6;

	std::valarray<T> beta_IS_coefs(substencil_order);
	std::valarray<T> c_k_coefs(substencil_order);

	T f_hat = 0.;
	auto f_stencil_view = std::ranges::views::all(f_stencil);

	auto a_func = [](auto substencilv) -> T {
		return -substencilv[0]
					+ substencilv[1]
					+ 2. * substencilv[2]
					- 2. * substencilv[3]
					- substencilv[4]
					+ substencilv[5];
	};

	auto b_func = [](auto substencilv) -> T {
		return substencilv[0]
					- 3. * substencilv[1]
					+ 2. * substencilv[2]
					+ 2. * substencilv[3]
					- 3. * substencilv[4]
					+ substencilv[5];
	};

	auto c_func = [](auto substencilv) -> T {
		return -substencilv[0]
					+ 5. * substencilv[1]
					- 10. * substencilv[2]
					+ 10. * substencilv[3]
					- 5. * substencilv[4]
					+ substencilv[5];
	};

	T a;
	T b;
	for (std::size_t k = 0;
			k <= std::ranges::size(f_stencil) - substencil_order; ++ k) {
		auto substencil = f_stencil_view
				| std::ranges::views::drop(k)
				| std::ranges::views::take(substencil_order);

		a = a_func(substencil);
		b = b_func(substencil);
		c_k_coefs[k] = c_func(substencil);

		beta_IS_coefs[k] = std::pow(b, 2) + std::abs(a * c_k_coefs[k]);
	}

	T tau_s = std::pow(b_func(c_k_coefs), 2)
				+ std::abs(a_func(c_k_coefs) * c_func(c_k_coefs));

	// non-linear non-scaled (α-)weights
	std::valarray<T> alpha_weights = 1.
				+ std::pow(tau_s / (eps + beta_IS_coefs), p);

	std::array<T, substencil_order> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<T> omega_weights = omegaWENO11FMWeights<T>(
				std::move(lambda_weights));

	std::valarray<T> eno_reconstructed_f
			= f6OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2]
			+ omega_weights[3] * eno_reconstructed_f[3]
			+ omega_weights[4] * eno_reconstructed_f[4]
			+ omega_weights[5] * eno_reconstructed_f[5];

	return f_hat;
}


template <ArithmeticWith<numeric_val> T>
T computeFHatWENO5ZMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		T eps = 1e-40, T p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 *
	 * Borges et al.'s WENO5-Z. Then in improves this by the symmetric
	 * mapping of Hong et al.
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// and in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)). For WENO7-Z(M)
	// p = 4.
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<T> beta_IS_coefs(3);

	T f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators<T>(f_stencil);
	T tau_5 = std::abs(beta_IS_coefs[2] - beta_IS_coefs[0]);

	std::array<T, 3> alpha_weights;
	std::ranges::transform(
				beta_IS_coefs,
				std::ranges::begin(alpha_weights),
				[tau_5, eps, p](auto beta) {
		return alphaWENO5ZMWeight(beta, tau_5, eps, p);
	});

	std::array<T, 3> lambda_weights;
	lambdaWENO5FMWeights<T>(std::move(alpha_weights), lambda_weights);

	std::valarray<T> omega_weights = omegaWENO5FMWeights<T>(
				std::move(lambda_weights));

	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<T> eno_reconstructed_f
			= f3OrdReconstructionFromStencil<T>(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


//template <ArithmeticWith<numeric_val> T>
//T computeFHatWENO5FMReconstructionKernelRev(std::span<T, 5> f_stencil,
//											T eps = 1e-40, T p = 2.) {
//	/* Calculate (reconstruct) one of the two split monotone numerical
//	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
//	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
//	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
//	 *                      ^    ^    ^    ^    ^    ^
//	 *                      0    1    2    3    4    |
//	 * in either case for convenience).
//	 *
//	 * I.e. this function implements the downwind reconstruction which
//	 * should be used for negative fluxes (with information propagated
//	 * from right to left) if the nodes are passed in order.
//	 */

//	std::valarray<T> beta_IS_coefs(3);

//	T f_hat = 0.;

//	// smoothness indicators of the stencil
//	// (measure how smooth u is in the stencil)
////	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

//	// The non-matrix variant seems to be faster(?)
//	// beta_IS_coefs = betaSmoothnessIndicators(f_stencil);

//	T f_prev2 = f_stencil[0];
//	T f_prev1 = f_stencil[1];
//	T f_curr0 = f_stencil[2];
//	T f_next1 = f_stencil[3];
//	T f_next2 = f_stencil[4];
////	T f_prev2 = f_stencil[4];
////	T f_prev1 = f_stencil[3];
////	T f_curr0 = f_stencil[2];
////	T f_next1 = f_stencil[1];
////	T f_next2 = f_stencil[0];

//	beta_IS_coefs[0] = ((13./12.) * std::pow(
//							f_prev2 - 2.*f_prev1 + f_curr0, 2)
//				+ (1./4.) * std::pow(
//							f_prev2 - 4.*f_prev1 + 3.*f_curr0, 2));

//	beta_IS_coefs[1] = ((13./12.) * std::pow(
//							f_prev1 - 2.*f_curr0 + f_next1, 2)
//				+ (1./4.) * std::pow((f_/usr/bin/prev1 - f_next1), 2));

//	beta_IS_coefs[2] = ((13./12.) * std::pow(
//							f_curr0 - 2.*f_next1 + f_next2, 2)
//				+ (1./4.) * std::pow(
//							3.*f_curr0 - 4.*f_next1 + f_next2, 2));

//	// non-linear non-scaled (α-)weights
//	std::valarray<T> d_lin_weights = {0.3, 0.6, 0.1};
//	std::valarray<T> alpha_weights = d_lin_weights
//			/ std::pow(eps + beta_IS_coefs, p);
//	// — we have no need of them in WENO-FM!

//	// Instead we use:
////	std::valarray<T> alpha_weights = alphaWENO5FMWeights(
////		std::move(beta_IS_coefs), eps, p
////	);

//	// scaled (normalized) non-linear (ω-)weights (ENO weights)
//	 std::valarray<T> omega_weights = alpha_weights
//			 / alpha_weights.sum();

//	// FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
//	// due to Zheng Hong, Zhengyin Ye and Kun Ye:
////	 lambda_weights = alpha_weights / alpha_weights.sum();
////	std::valarray<T> lambda_weights = lambdaWENO5FMWeights(
////		std::move(alpha_weights)
////	);

////	std::valarray<T> omega_weights = omegaWENO5FMWeights(
////		std::move(lambda_weights)
////	);

//	// vecMatDot<T>(u_..., WmN...) stores a 3-rd order estimate of
//	// f_{i+1/2} via linear combinations with WmNplus coefficients
//	// for each substencil which is then used to calculate
//	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
//	// using the nonlinear weights [ω]
////	f_hat = std::inner_product(
////		std::ranges::begin(omega_weights),
////		std::ranges::end(omega_weights),
////		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)),
////		0.
////	);

//	 std::valarray<
//		 T
//	 > eno_reconstructed_f(3);
//	 eno_reconstructed_f[0] = (-1. * f_stencil[0]
//							   + 5. * f_stencil[1]
//							   + 2. * f_stencil[2]) / 6.;

//	 eno_reconstructed_f[1] = (2. * f_stencil[1]
//							   + 5. * f_stencil[2]
//							   - 1. * f_stencil[3]) / 6.;

//	 eno_reconstructed_f[2] = (11. * f_stencil[2]
//							   - 7. * f_stencil[3]
//							   + 2. * f_stencil[4]) / 6.;

//	f_hat = omega_weights[0] * eno_reconstructed_f[0]
//			+ omega_weights[1] * eno_reconstructed_f[1]
//			+ omega_weights[2] * eno_reconstructed_f[2];

//	return f_hat;
//}


// FD WENO5FM (WENO5-FM) - method
// Reconstruction based on LF flux splitting + improved mapped WENO
// of 5th order
// (see Mapped weighted essentially non-oscillatory schemes:
// achieving optimal order near critical points, 2005 by Henrick et al.)
// and 'An improved WENO-Z scheme with symmetry-preserving mapping'
// by Zheng Hong, Zhengyin Ye and Kun Ye, 2020
template <ArithmeticWith<numeric_val> T, std::size_t N>
void calcHydroStageFDWENO(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		auto&& computeWENOReconstructionKernel,
		std::size_t n_ghost_cells = (N + 1) / 2,
		T eps = 1e-40,
		T p = 2.) {
	/* Component-wise finite-difference WENO5FM (FD WENO5-FM) - space
	 * reconstruction method with the global Lax-Friedrichs (LF) flux
	 * splitting.
	 *
	 * Usually, componentwise reconstruction produces satisfactory
	 * results for schemes up to third-order accuracy, while characteristic
	 * reconstruction produces better nonoscillatory results for
	 * higher-order accuracy, albeit with an increased computational cost.
	 */

	const unsigned order = N;
	assert(order % 2 != 0);
	const std::size_t stencil_size = order;
	const std::size_t _actual_stencil_size = stencil_size + 1;
	const std::size_t half_size = order / 2;

	const std::size_t r = _actual_stencil_size / 2;
	assert(n_ghost_cells >= r);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(
				numerical_flux) - n_ghost_cells - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	auto shifted_index_range = std::ranges::common_view(
				std::views::iota(mini - 1)
					| std::views::take(maxi + 1 - (mini - 1)/* + 1*/));

	T fhatminus = 0.;
	T fhatplus = 0.;

	auto j_it_p = std::ranges::begin(f_plus);  // f_plus
	auto j_it_m = std::ranges::begin(f_minus);  // f_minus

	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size);
	std::advance(j_it_m, mini - 1 + half_size + 1 - stencil_size);
	auto u_plus = std::ranges::views::counted(
				j_it_p, _actual_stencil_size);
	auto u_minus = std::ranges::views::counted(
				j_it_m, _actual_stencil_size)
					| std::ranges::views::reverse;

	std::for_each(std::execution::par_unseq,
				std::ranges::begin(shifted_index_range),
				std::ranges::end(shifted_index_range),
				  [&](std::size_t j) {
		j_it_p = std::ranges::begin(f_plus);  // f_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size);
		u_plus = std::ranges::views::counted(j_it_p, _actual_stencil_size);

		j_it_m = std::ranges::begin(f_minus);  // f_minus
		std::advance(j_it_m, j + half_size + 1 - stencil_size);
		u_minus = std::ranges::views::counted(j_it_m, _actual_stencil_size)
					| std::ranges::views::reverse;

		fhatplus = computeWENOReconstructionKernel(
			std::ranges::views::counted(
						std::ranges::begin(u_plus), stencil_size), eps, p
		);

		fhatminus = computeWENOReconstructionKernel(
			std::ranges::subrange(
				std::ranges::begin(u_minus),
						std::ranges::end(u_minus) - 1), eps, p
		);

		numerical_flux[j] = fhatplus + fhatminus;
	});
}


// FD WENO5FM (WENO5-FM) - method
// Reconstruction based on LF flux splitting + improved mapped WENO
// of 5th order
// (see Mapped weighted essentially non-oscillatory schemes:
// achieving optimal order near critical points, 2005 by Henrick et al.)
// and 'An improved WENO-Z scheme with symmetry-preserving mapping'
// by Zheng Hong, Zhengyin Ye and Kun Ye, 2020
template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO5(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		auto&& computeWENOReconstructionKernel,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	/* Component-wise finite-difference WENO5FM (FD WENO5-FM) - space
	 * reconstruction method with the global Lax-Friedrichs (LF) flux
	 * splitting.
	 *
	 * Usually, componentwise reconstruction produces satisfactory
	 * results for schemes up to third-order accuracy, while characteristic
	 * reconstruction produces better nonoscillatory results for
	 * higher-order accuracy, albeit with an increased computational cost.
	 */

	const unsigned order = 5;
	const std::size_t stencil_size = order;
	const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = order / 2;  // 2

	// r = (order + 1) / 2 = 3
	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(
				numerical_flux) - n_ghost_cells - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	auto shifted_index_range = std::ranges::common_view(
				std::views::iota(mini - 1)
					| std::views::take(maxi + 1 - (mini - 1) + 1));
	// [g      g      g      i      i      i      i      i      i      ...]
	// {0      1      2      3      4      5}     6      7      8      ...
	//  |             |      |      |
	// itr            j nGhostCells end()

	// WENO5 stencils

	// Coefficients WmN(+/-) before fluxes at the stencil nodes
	// to find component stencils
	// [q_k] = WmN(+/-) * [ f[j-2] ... f[j+2] ]
	// of the WENO interpolator.
	// So 3rd order approximation coefficients associated with
	// the substencils.

	// Calculation of f_hat, the numerical flux of u (whichever
	// is chosen), requires the approximation of u that uses at
	// the j-th cell of u-discretization a group of cell average
	// values on the left (`u_minus`) and on the right (`u_plus`).
	// So left- and right-biased approximations respectively.
	// `u_plus` represents the cells [j-2, j-1, j, j+1, j+2],
	// and `u_minus` represents the cells [j-1, j, j+1, j+2, j+3];
	// for convenience and uniformity we represent both using the
	// same combined structure of    [j-2, j-1, j, j+1, j+2, j+3].
	// std::valarray<T> u_plus(_actual_stencil_size);   // f_plus
	// std::valarray<T> u_minus(_actual_stencil_size);  // f_minus

	// For the purpose of linear stability (upwinding),
	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	// dfminus/du <= 0), is performed.
	// Lax-Friedrichs (LF) flux splitting (because it is the simplest
	// and smooth - we need the positive and negative fluxes to have
	// as many derivatives as the order of our finite-difference WENO)
	// is chosen here. `f_plus` uses a biased stencil with 1 point to
	// the left, while `f_minus` uses a biased stencil with 1 point to
	// the right.

	T fhatminus = 0.;
	T fhatplus = 0.;

	// So an LF flux	`numerical_flux`, f_hat(u_minus, u_plus),
	// a monotone numerical flux consistent with the physical one
	// (f_hat(u, u) = f(u)), will be construced at the end of
	// the loop below from `fhatminus` and `fhatplus` using
	// `f_minus` and `f_plus`.
	// N.B.! This can be replaced by an exact or approximate
	// Riemann solver (see Toro, 2009). Somehow...
	// Not every monotone flux can be writtenin the flux split form.
	// For example, the Godunov flux cannot.
	// std::valarray<T> numerical_flux(0., u.size());
	// f = std::valarray<T>(u.size());

	auto j_it_p = std::ranges::begin(f_plus);  // f_plus
	auto j_it_m = std::ranges::begin(f_minus);  // f_minus

//	std::ranges::transform(
//			monotone_flux_components[0],
//			monotone_flux_components[1],
//			std::ranges::begin(f), [](const auto fp, const auto fm) {

//	})
	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size);
	std::advance(j_it_m, mini - 1 + half_size + 1 - stencil_size);
	auto u_plus = std::ranges::views::counted(j_it_p, 6);
	auto u_minus = std::ranges::views::counted(j_it_m, 6)
					| std::ranges::views::reverse;

	std::for_each(std::execution::par_unseq,
				std::ranges::begin(shifted_index_range),
				std::ranges::end(shifted_index_range),
				  [&](std::size_t j) {
		j_it_p = std::ranges::begin(f_plus);  // f_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size);
		u_plus = std::ranges::views::counted(j_it_p, 6);

		j_it_m = std::ranges::begin(f_minus);  // f_minus
		std::advance(j_it_m, j + half_size + 1 - stencil_size);
		u_minus = std::ranges::views::counted(j_it_m, 6)
					| std::ranges::views::reverse;

		fhatplus = computeWENOReconstructionKernel(
			std::ranges::views::counted(
						std::ranges::begin(u_plus), 5), eps, p
		);

		fhatminus = computeWENOReconstructionKernel(
			std::ranges::subrange(
				std::ranges::begin(u_minus),
						std::ranges::end(u_minus) - 1), eps, p
		);
//		fhatminus = computeFHatWENO5JSReconstructionKernelRev<T>(
//			std::ranges::views::counted(
//						std::ranges::begin(u_minus)+1, 5), eps, p
//		);

		numerical_flux[j] = fhatplus + fhatminus;
	});

//	std::valarray<T> u_plus(6);
//	std::valarray<T> u_minus(6);
//	for (std::size_t j = mini - 1; j < maxi + 1; ++ j) {
//		u_plus[0] = f_plus[j - 2];
//		u_plus[1] = f_plus[j - 1];
//		u_plus[2] = f_plus[j + 0];
//		u_plus[3] = f_plus[j + 1];
//		u_plus[4] = f_plus[j + 2];

//		u_minus[0] = f_minus[j + 3];
//		u_minus[1] = f_minus[j + 2];
//		u_minus[2] = f_minus[j + 1];
//		u_minus[3] = f_minus[j + 0];
//		u_minus[4] = f_minus[j - 1];

//		numerical_flux[j] = computeFHatWENO5MReconstructionKernel<T>(
//					std::ranges::views::all(u_plus), eps, p)
//				+ computeFHatWENO5MReconstructionKernel<T>(
//					std::ranges::views::all(u_minus), eps, p);
//	}

	// return numerical_flux;
	// std::cout << " done!" << "\n";std::ranges::views::all|
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO5JS(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO5<T>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO5JSReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO7BS(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 4,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO<T, 7>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO7BSReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO5M(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO5<T>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO5MReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO9M(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 5,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO<T, 9>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO9MReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


//template <ArithmeticWith<numeric_val> T>
//void calcHydroStageFDWENO5FM(
//		const std::ranges::common_range auto&& f_plus,
//		const std::ranges::common_range auto&& f_minus,
//		T t,
//		std::ranges::common_range auto&& numerical_flux,
//		std::size_t n_ghost_cells = 3,
//		T eps = 1e-40,
//		T p = 2.) {
//	calcHydroStageFDWENO5<T>(
//				std::ranges::views::all(f_plus),
//				std::ranges::views::all(f_minus), t,
//				std::ranges::views::all(numerical_flux),
//				[](const std::ranges::sized_range auto&& stencil,
//							T eps, T p) -> T {
//					return computeFHatWENO5FMReconstructionKernel<T>(
//								std::ranges::views::all(stencil), eps, p);
//				}, n_ghost_cells, eps, p);
//}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO5FM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO<T, 5>(
				std::move(f_plus),
				std::move(f_minus), t,
				std::move(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO5FMReconstructionKernel<T>(
								std::move(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO7FM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 4,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO<T, 7>(
				std::move(f_plus),
				std::move(f_minus), t,
				std::move(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO7FMReconstructionKernel<T>(
								std::move(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO5ZM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO5<T>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO5ZMReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO7Z(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 4,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO<T, 7>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO7ZReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO7ZM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 4,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO<T, 7>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO7ZMReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO7S(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 4,
		T eps = 1e-80,
		T p = 1.) {
	calcHydroStageFDWENO<T, 7>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO7SReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO7SM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 4,
		T eps = 1e-80,
		T p = 1.) {
	calcHydroStageFDWENO<T, 7>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO7SMReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO9S(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 5,
		T eps = 1e-80,
		T p = 1.) {
	calcHydroStageFDWENO<T, 9>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO9SReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO9SM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 5,
		T eps = 1e-80,
		T p = 1.) {
	calcHydroStageFDWENO<T, 9>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO9SMReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO11S(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 6,
		T eps = 1e-80,
		T p = 1.) {
	calcHydroStageFDWENO<T, 11>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO11SReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO11SM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 6,
		T eps = 1e-100,
		T p = 2.) {
	calcHydroStageFDWENO<T, 11>(
				std::ranges::views::all(f_plus),
				std::ranges::views::all(f_minus), t,
				std::ranges::views::all(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO11SMReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFVWENO5(
		const std::ranges::common_range auto&& u,
		T t,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		auto&& computeWENOReconstructionKernel,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	/* Component-wise finite-volume WENO5FM (FV WENO5-FM) - space
	 * reconstruction method with the global Lax-Friedrichs (LF) flux
	 * splitting.
	 *
	 * Usually, componentwise reconstruction produces satisfactory
	 * results for schemes up to third-order accuracy, while characteristic
	 * reconstruction produces better nonoscillatory results for
	 * higher-order accuracy, albeit with an increased computational cost.
	 */

	const unsigned order = 5;
	const std::size_t stencil_size = order;
	// const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = order / 2;  // 2

	// r = (order + 1) / 2 = 3
	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;  // 3
	const std::size_t mini = n_ghost_cells;  // at least 3
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(u) - n_ghost_cells - 1;
	auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	// [g      g      g      i      i      i      i      i      i      ...]
	// {0      1      2      3      4      5}     6      7      8      ...
	//  |             |      |      |
	// itr            j nGhostCells end()

	// WENO5 stencils

	// Coefficients WmN(+/-) before fluxes at the stencil nodes
	// to find component stencils
	// [q_k] = WmN(+/-) * [ f[j-2] ... f[j+2] ]
	// of the WENO interpolator.
	// So 3rd order approximation coefficients associated with
	// the substencils.

	// Calculation of f_hat, the numerical flux of u (whichever
	// is chosen), requires the approximation of u that uses at
	// the j-th cell of u-discretization a group of cell average
	// values on the left (`u_minus`) and on the right (`u_plus`).
	// So left- and right-biased approximations respectively.
	// `u_plus` represents the cells [j-2, j-1, j, j+1, j+2],
	// and `u_minus` represents the cells [j-1, j, j+1, j+2, j+3];
	// for convenience and uniformity we represent both using the
	// same combined structure of    [j-2, j-1, j, j+1, j+2, j+3].
	// std::valarray<double> u_plus(_actual_stencil_size);   // f_plus
	// std::valarray<double> u_minus(_actual_stencil_size);  // f_minus

	// For the purpose of linear stability (upwinding),
	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	// dfminus/du <= 0), is performed.
	// Lax-Friedrichs (LF) flux splitting (because it is the simplest
	// and smooth - we need the positive and negative fluxes to have
	// as many derivatives as the order of our finite-difference WENO)
	// is chosen here. `f_plus` uses a biased stencil with 1 point to
	// the left, while `f_minus` uses a biased stencil with 1 point to
	// the right.

	// So an LF flux	`numerical_flux`, f_hat(u_minus, u_plus),
	// a monotone numerical flux consistent with the physical one
	// (f_hat(u, u) = f(u)), will be construced at the end of
	// the loop below from `fhatminus` and `fhatplus` using
	// `f_minus` and `f_plus`.
	// N.B.! This can be replaced by an exact or approximate
	// Riemann solver (see Toro, 2009). Somehow...
	// Not every monotone flux can be writtenin the flux split form.
	// For example, the Godunov flux cannot.
	// std::valarray<double> numerical_flux(0., u.size());
	// f = std::valarray<double>(u.size());

	auto j_it_p = std::ranges::begin(u);  // f_plus

//	std::ranges::transform(
//			monotone_flux_components[0],
//			monotone_flux_components[1],
//			std::ranges::begin(f), [](const auto fp, const auto fm) {

//	})
	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size);
	auto u_plus = std::ranges::views::counted(j_it_p, 6);
	auto u_minus = u_plus | std::ranges::views::reverse;

	for (std::size_t j : shifted_index_range) {
		j_it_p = std::ranges::begin(u);  // f_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size);
		u_plus = std::ranges::views::counted(j_it_p, 6);

		u_plus_rec[j] = computeWENOReconstructionKernel(
			std::ranges::views::counted(
						std::ranges::begin(u_plus), 5), eps, p
		);


		u_minus = u_plus | std::ranges::views::reverse;
		u_minus_rec[j] = computeWENOReconstructionKernel(
			std::ranges::subrange(
				std::ranges::begin(u_minus),
						std::ranges::end(u_minus) - 1), eps, p
		);
//		u_minus_rec[j] = computeFHatWENO5JSReconstructionKernelRev(
//			std::ranges::views::counted(
//						std::ranges::begin(u_minus)+1, 5), eps, p
//		);
	}

	// std::cout << " done!" << "\n";
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFVWENO5JS(
		const std::ranges::common_range auto&& u,
		T t,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFVWENO5<T>(
				std::ranges::views::all(u), t,
				std::ranges::views::all(u_plus_rec),
				std::ranges::views::all(u_minus_rec),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO5JSReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				},
				n_ghost_cells,
				eps,
				p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFVWENO5M(
		const std::ranges::common_range auto&& u,
		T t,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFVWENO5<T>(
				std::ranges::views::all(u), t,
				std::ranges::views::all(u_plus_rec),
				std::ranges::views::all(u_minus_rec),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO5MReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				},
				n_ghost_cells,
				eps,
				p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFVWENO5FM(
		const std::ranges::common_range auto&& u,
		T t,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFVWENO5<T>(
				std::ranges::views::all(u), t,
				std::ranges::views::all(u_plus_rec),
				std::ranges::views::all(u_minus_rec),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO5FMReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				},
				n_ghost_cells,
				eps,
				p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFVWENO5ZM(
		const std::ranges::common_range auto&& u,
		T t,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFVWENO5<T>(
				std::ranges::views::all(u), t,
				std::ranges::views::all(u_plus_rec),
				std::ranges::views::all(u_minus_rec),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO5ZMReconstructionKernel<T>(
								std::ranges::views::all(stencil), eps, p);
				},
				n_ghost_cells,
				eps,
				p);
}


template <ArithmeticWith<numeric_val> T, ArithmeticWith<numeric_val> VT>
void calcHydroStageCharWiseFDWENO5FM(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		const std::ranges::common_range auto&& flux,
		std::ranges::common_range auto&& numerical_flux,
		T t,
		auto&& project,
//		auto&& deproject,
		T alpha,
		std::size_t n_ghost_cells = 3,
		T eps = 1e-40,
		T p = 2.) {
	const unsigned order = 5;
	const std::size_t stencil_size = order;
	const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = order / 2;  // 2

	// r = (order + 1) / 2 = 3
	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(
				numerical_flux) - n_ghost_cells - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	auto shifted_index_range = std::ranges::common_view(
				std::views::iota(mini - 1)
					| std::views::take(maxi + 1 - (mini - 1) + 1));
	VT fhatminus = static_cast<VT>(0.);
	VT fhatplus = static_cast<VT>(0.);

	auto j_it_q = std::ranges::begin(u);
	auto j_it_f = std::ranges::begin(flux);

	std::advance(j_it_q, mini - 1 + half_size + 1 - stencil_size);
	std::advance(j_it_f, mini - 1 + half_size + 1 - stencil_size);
	auto f_stencil = std::ranges::views::counted(j_it_f, 6);
	auto q_stencil = std::ranges::views::counted(j_it_q, 6);

	std::valarray<VT> u_plus(_actual_stencil_size);
	std::valarray<VT> u_minus(_actual_stencil_size);

	auto components = {
		std::make_pair(0, &Vector4<T>::x),
		std::make_pair(1, &Vector4<T>::y),
		std::make_pair(2, &Vector4<T>::z)
//		&Vector4<T>::w
	};

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(shifted_index_range),
				std::ranges::end(shifted_index_range),
				[&](std::size_t j) {
		j_it_q = std::ranges::begin(u);
		j_it_f = std::ranges::begin(flux);
		std::advance(j_it_q, j + half_size + 1 - stencil_size);
		q_stencil = std::ranges::views::counted(j_it_q, 6);

		std::advance(j_it_f, j + half_size + 1 - stencil_size);
		f_stencil = std::ranges::views::counted(j_it_f, 6);

		auto proj_u_j = [j, &q_avg, &project](auto u) -> decltype(u) {
			const auto q = q_avg[j];
			return project(q, u);
		};

//		auto deproj_u_j = [j, &q_avg, &deproject](auto u) -> decltype(u) {
//			const auto q = q_avg[j];
//			return deproject(q, u);
//		};

		// u_plus[0] = proj_u_j(f_stencil[0]);
		// std::cout << u_plus[0] << "\n";

		for (std::size_t k = 0; k < 6; ++ k)
			u_plus[k] = (proj_u_j(f_stencil[k])
					+ /*1.1 * */alpha * proj_u_j(q_stencil[k])) * 0.5;

		for (std::size_t k = 0; k < 6; ++ k)
			u_minus[k] = (proj_u_j(f_stencil[6 - k - 1])
					- /*1.1 * */alpha * proj_u_j(q_stencil[6 - k - 1])) * 0.5;

		for (auto& comp : components) {
			fhatplus[comp.first] = computeFHatWENO5FMReconstructionKernel<T>(
				std::ranges::views::counted(
							std::ranges::begin(u_plus), 5)
						| std::ranges::views::transform(
							comp.second), eps, p
			);
//			if (j == 3)
//			std::cout << t << " " << j << " ["
//					<< q_avg[j] << "\n"
//					<< proj_u_j(f_stencil[0]) << "\n"
//					<< proj_u_j(q_stencil[0]) << "\n"
//					<< alpha << "\n"
//					<< alpha * proj_u_j(q_stencil[0]) << "\n"
//					<< u_plus[0] << "\n"
//					<< fhatplus << "\n"
//					<< deproj_u_j(fhatplus)
//					<< "]" << "\n\n";

			fhatminus[comp.first] = computeFHatWENO5FMReconstructionKernel<T>(
				std::ranges::views::counted(
							std::ranges::begin(u_minus), 5)
						| std::ranges::views::transform(
							comp.second), eps, p
			);
		}


		numerical_flux[j] = fhatplus + fhatminus;
	});
	// std::cout << numerical_flux[3] << "\n";
}


template <ArithmeticWith<numeric_val> T, ArithmeticWith<numeric_val> VT,
			std::size_t N>
void calcHydroStageCharWiseFDWENO(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		const std::ranges::common_range auto&& flux,
		std::ranges::common_range auto&& numerical_flux,
		T t,
		auto&& project,
//		auto&& deproject,
		T alpha,
		auto&& computeWENOReconstructionKernel,
		std::size_t n_ghost_cells = (N + 1) / 2,
		T eps = 1e-40,
		T p = 2.) {
	const unsigned order = N;
	const std::size_t stencil_size = order;
	const std::size_t _actual_stencil_size = stencil_size + 1;
	const std::size_t half_size = order / 2;

	const std::size_t r = (order + 1) / 2;
	assert(n_ghost_cells >= r);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;
	const std::size_t mini = n_ghost_cells;
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(
				numerical_flux) - n_ghost_cells - 1;
	// auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	auto shifted_index_range = std::ranges::common_view(
				std::views::iota(mini - 1)
					| std::views::take(maxi + 1 - (mini - 1)/* + 1*/));
	VT fhatminus = static_cast<VT>(0.);
	VT fhatplus = static_cast<VT>(0.);

	auto j_it_q = std::ranges::begin(u);
	auto j_it_f = std::ranges::begin(flux);

	std::advance(j_it_q, mini - 1 + half_size + 1 - stencil_size);
	std::advance(j_it_f, mini - 1 + half_size + 1 - stencil_size);
	auto f_stencil = std::ranges::views::counted(j_it_f,
												 _actual_stencil_size);
	auto q_stencil = std::ranges::views::counted(j_it_q,
												 _actual_stencil_size);

	std::valarray<VT> u_plus(_actual_stencil_size);
	std::valarray<VT> u_minus(_actual_stencil_size);

	auto components = {
		std::make_pair(0, &Vector4<T>::x),
		std::make_pair(1, &Vector4<T>::y),
		std::make_pair(2, &Vector4<T>::z)
//		&Vector4<T>::w
	};

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(shifted_index_range),
				std::ranges::end(shifted_index_range),
				[&](std::size_t j) {
		j_it_q = std::ranges::begin(u);
		j_it_f = std::ranges::begin(flux);
		std::advance(j_it_q, j + half_size + 1 - stencil_size);
		q_stencil = std::ranges::views::counted(j_it_q,
												_actual_stencil_size);

		std::advance(j_it_f, j + half_size + 1 - stencil_size);
		f_stencil = std::ranges::views::counted(j_it_f,
												_actual_stencil_size);

		auto proj_u_j = [j, &q_avg, &project](auto u) -> decltype(u) {
			const auto q = q_avg[j];
			return project(q, u);
		};

		for (std::size_t k = 0; k < _actual_stencil_size; ++ k)
			u_plus[k] = (proj_u_j(f_stencil[k])
					+ /*1.1 * */alpha * proj_u_j(q_stencil[k])) * 0.5;

		for (std::size_t k = 0; k < _actual_stencil_size; ++ k)
			u_minus[k] = (proj_u_j(f_stencil[_actual_stencil_size - k - 1])
					- /*1.1 * */alpha * proj_u_j(
						q_stencil[_actual_stencil_size - k - 1])) * 0.5;

		for (auto& comp : components) {
			fhatplus[comp.first] = computeWENOReconstructionKernel(
				std::ranges::views::counted(
							std::ranges::begin(u_plus), order)
						| std::ranges::views::transform(
							comp.second), eps, p
			);

			fhatminus[comp.first] = computeWENOReconstructionKernel(
				std::ranges::views::counted(
							std::ranges::begin(u_minus), order)
						| std::ranges::views::transform(
							comp.second), eps, p
			);
		}

		numerical_flux[j] = fhatplus + fhatminus;
	});
}


template <ArithmeticWith<numeric_val> T, ArithmeticWith<numeric_val> VT>
void calcHydroStageCharWiseFDWENO7FM(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		const std::ranges::common_range auto&& flux,
		std::ranges::common_range auto&& numerical_flux,
		T t,
		auto&& project,
//		auto&& deproject,
		T alpha,
		std::size_t n_ghost_cells = 4,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageCharWiseFDWENO<T, VT, 7>(
				std::move(u),
				std::move(q_avg),
				std::move(flux),
				std::move(numerical_flux),
				t, project,
				// auto&& deproject,
				alpha,
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO7FMReconstructionKernel<T>(
								std::move(stencil), eps, p);
				},
				n_ghost_cells,
				eps,
				p);
}


template <ArithmeticWith<numeric_val> T, ArithmeticWith<numeric_val> VT>
void calcHydroStageCharWiseFDWENO9BS(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		const std::ranges::common_range auto&& flux,
		std::ranges::common_range auto&& numerical_flux,
		T t,
		auto&& project,
//		auto&& deproject,
		T alpha,
		std::size_t n_ghost_cells = 5,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageCharWiseFDWENO<T, VT, 9>(
				std::move(u),
				std::move(q_avg),
				std::move(flux),
				std::move(numerical_flux),
				t, project,
				// auto&& deproject,
				alpha,
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO9BSReconstructionKernel<T>(
								std::move(stencil), eps, p);
				},
				n_ghost_cells,
				eps,
				p);
}


template <ArithmeticWith<numeric_val> T, ArithmeticWith<numeric_val> VT>
void calcHydroStageCharWiseFDWENO9FM(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		const std::ranges::common_range auto&& flux,
		std::ranges::common_range auto&& numerical_flux,
		T t,
		auto&& project,
//		auto&& deproject,
		T alpha,
		std::size_t n_ghost_cells = 5,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageCharWiseFDWENO<T, VT, 9>(
				std::move(u),
				std::move(q_avg),
				std::move(flux),
				std::move(numerical_flux),
				t, project,
				// auto&& deproject,
				alpha,
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO9FMReconstructionKernel<T>(
								std::move(stencil), eps, p);
				},
				n_ghost_cells,
				eps,
				p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO9FM(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 5,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO<T, 9>(
				std::move(f_plus),
				std::move(f_minus), t,
				std::move(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO9FMReconstructionKernel<T>(
								std::move(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


template <ArithmeticWith<numeric_val> T>
void calcHydroStageFDWENO9BS(
		const std::ranges::common_range auto&& f_plus,
		const std::ranges::common_range auto&& f_minus,
		T t,
		std::ranges::common_range auto&& numerical_flux,
		std::size_t n_ghost_cells = 5,
		T eps = 1e-40,
		T p = 2.) {
	calcHydroStageFDWENO<T, 9>(
				std::move(f_plus),
				std::move(f_minus), t,
				std::move(numerical_flux),
				[](const std::ranges::sized_range auto&& stencil,
							T eps, T p) -> T {
					return computeFHatWENO9BSReconstructionKernel<T>(
								std::move(stencil), eps, p);
				}, n_ghost_cells, eps, p);
}


// template <ArithmeticWith<numeric_val> T>
void updateGhostPointsTransmissive(
		std::ranges::common_range auto&& U,
		std::size_t left_bound_size = 3,
		std::optional<std::size_t> right_bound_size = std::nullopt) {
	/* Update ghost points in U with transmissive (Neumann) b.c.s. */

	if (!right_bound_size)
		right_bound_size.emplace(left_bound_size);

	const std::size_t n_full_size = std::ranges::size(U);
	const std::size_t right_start_index = (n_full_size
										   - right_bound_size.value());
	auto left_boundary_start = std::ranges::begin(U);
	auto right_boundary_start = std::ranges::begin(U);
	std::advance(right_boundary_start, right_start_index);

	// Transmissive b.c.s
	std::for_each_n(/*std::execution::par_unseq,*/
					left_boundary_start,
					left_bound_size,
					[&U, left_bound_size](auto& n) {
		n = U[left_bound_size];
	});
//	const std::size_t mini = left_bound_size;
//	U[2] = U[mini]; U[1] = U[mini]; U[0] = U[mini];

	std::for_each_n(/*std::execution::par_unseq,*/
					right_boundary_start,
					right_bound_size.value(),
					[&U, right_start_index](auto& n) {
		n = U[right_start_index-1];
	});
//	const std::size_t maxi = n_full_size - right_bound_size.value() - 1;
//	U[maxi+1] = U[maxi]; U[maxi+2] = U[maxi]; U[maxi+3] = U[maxi];
}


// template <ArithmeticWith<numeric_val> T>
void updateGhostPointsPeriodic(
		std::ranges::common_range auto&& U,
		std::size_t left_bound_size = 5,
		std::optional<std::size_t> right_bound_size = std::nullopt) {
	/* Update ghost points in U with periodic b.c.s. */

	if (!right_bound_size)
		right_bound_size.emplace(left_bound_size);

	const std::size_t n_full_size = std::ranges::size(U);
//	const std::size_t right_start_index = (n_full_size
//										   - right_bound_size.value());
//	const std::size_t left_start_index = left_bound_size;

	// Periodic b.c.s
//	std::ranges::transform(
//				U | std::ranges::views::take(left_bound_size),
//				std::ranges::begin(U | std::ranges::views::drop(
//						right_start_index - right_bound_size.value())),
//				[](const auto& u) {
//		return u;
//	});
//	std::cout << "Left\n";
//	std::cout << U[0] << " " << U[1] << " " << U[2] << "\n";
//	std::cout << U[n_full_size-2-3-1]
//			<< " " << U[n_full_size-1-3-1]
//			<< " " << U[n_full_size-3-1] << "\n";
//	U[0] = U[n_full_size - 1 - 3 - 2];
//	U[1] = U[n_full_size - 1 - 3 - 1];
//	U[2] = U[n_full_size - 1 - 3 - 0];
	for (std::size_t k = 0; k < left_bound_size; ++ k)
		U[k] = U[n_full_size - 1 - left_bound_size
					- (left_bound_size - k)];
//	U[0] = U[n_full_size - 1 - 5 - 5];
//	U[1] = U[n_full_size - 1 - 5 - 4];
//	U[2] = U[n_full_size - 1 - 5 - 3];
//	U[3] = U[n_full_size - 1 - 5 - 2];
//	U[4] = U[n_full_size - 1 - 5 - 1];

//	std::ranges::transform(
//				U | std::ranges::views::drop(right_start_index),
//				std::ranges::begin(U | std::ranges::views::take(
//						left_bound_size + left_bound_size)),
//				[](const auto& u) {
//		return u;
//	});
//	std::cout << "Right\n";
//	std::cout << U[n_full_size-2-1]
//			<< " " << U[n_full_size-1-1]
//			<< " " << U[n_full_size-1] << "\n";
//	std::cout << U[3]
//			<< " " << U[4]
//			<< " " << U[5] << "\n";
//	U[n_full_size - 1 - 2] = U[3 + 0];
//	U[n_full_size - 1 - 1] = U[3 + 1];
//	U[n_full_size - 1 - 0] = U[3 + 2];
	for (std::size_t k = 0; k < right_bound_size.value(); ++ k)
		U[n_full_size - 1 - (right_bound_size.value() - k - 1)]
				= U[left_bound_size + k + 1];
//	U[n_full_size - 1 - 4] = U[5 + 1];
//	U[n_full_size - 1 - 3] = U[5 + 2];
//	U[n_full_size - 1 - 2] = U[5 + 3];
//	U[n_full_size - 1 - 1] = U[5 + 4];
//	U[n_full_size - 1 - 0] = U[5 + 5];
}


template <ArithmeticWith<numeric_val> T, typename... Args>
void timeOperator(
	std::ranges::common_range auto& U,
	std::ranges::common_range auto& flux,
	std::ranges::common_range auto& intermediate_fluxes,
	T t0, T dx, std::size_t n_ghost_points, T fin_t,
	auto&& timeStepFunction,
	auto&& calcMaxWaveSpeed,
	T cfl = 0.4,
	Args... opts_args
) {
	/* Time Operator of some problem: perform the time loop
	 * and solve it, storing the result in `U` and the numerical flux
	 * needed for the calculation in `flux`.
	 */

	T t = t0;
	T dt = 0.;

	T cpu = calcMaxWaveSpeed(U, dt);
	std::valarray<T> lam = std::valarray(cpu, 4);
	long long counter = 0;
	clock_t t_start = 0;
	clock_t t_end = 0;
	double t_calc = 0.;
	clock_t t_start_global = 0;
	clock_t t_end_global = 0;

	t_start_global = clock();
	while (t < fin_t) {
		t_start = clock();

		dt = cfl * dx / lam[0];

		if (counter < 5) {
			// cfl = cfl * .2;
			dt *= .2;
		}
		// dt = /*8. **/ cfl * std::pow(dx, 5./3.);
		// dt = 8. * std::pow(dx, 3.);;

		if (t + dt > fin_t)
			dt = fin_t - t;

		timeStepFunction(U, flux, intermediate_fluxes,
				t, dt, dx, lam, n_ghost_points, opts_args...);

		cpu = calcMaxWaveSpeed(U, dt);
		lam = std::valarray(cpu, 4);

		t_end = clock();

		t_calc = static_cast<double>(t_end - t_start) / CLOCKS_PER_SEC;
		std::cout << "iteration=" << counter << "\t"
				  << "t=" << t << "\t"
				  << "dt=" << dt << "\t"
				  << "CFL=" << cfl << "\t"
				  << std::setprecision(2)
				  << "exec. t=" << t_calc << "s" << "\n";
		t += dt;

		++ counter;
	}
	t_end_global = clock();
	t_calc = static_cast<double>(
				t_end_global - t_start_global) / CLOCKS_PER_SEC;
	std::cout << "total execution time=" << t_calc << "s" << "\n";
	// std::cout << "t = " << t << "\n";

	// return U;
}

#endif // WENO5_H
