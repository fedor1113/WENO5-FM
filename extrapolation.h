#include <algorithm>
#include <cassert>
#include <iostream>
#include <span>
#include <valarray>

// #include <Eigen/Dense>
// #include <Eigen/QR>

#include "Eigen/Dense"
#include "Eigen/QR"

#include "arithmeticwith.h"


template <ArithmeticWith<numeric_val> T, std::size_t N = 4>
void polyfit(std::span<T> const argument_data,
		std::span<T> const function_value_data,
		std::span<T, N + 1> fitted_polynomial_coefficients,
		Eigen::Matrix<T, Eigen::Dynamic, N + 1>& vandermonde_mat) {
	/* Efficient linear least square polynomial fit using
	 * Eigen for fast matrix operations.
	 *
	 * `order` = `N`;
	 * Fit a polynomial of degree `order` (4 by default,
	 * as this is usually sufficient, and a polynomial
	 * of too high a degree would pose a problem with
	 * the Vandermonde-matrix-computation approach taken
	 * here) to interpolate/extrapolate the function
	 * given by a discrete collection of its values
	 * `function_value_data` at nodes `argument_data`.
	 *
	 * Store the resulting fitted coeeficients
	 * at the first (`order` + 1) cells
	 * of the span `fitted_polynomial_coefficients`
	 * (which naturally should represent a big enough
	 * object to be able to store them).
	 *
	 * The coefficients are stored starting from the
	 * constant and ending with the leading coefficient,
	 * so that the fitted polynomial can be restored by
	 * `T f_x = fitted_polynomial_coefficients[0]
	 *        + fitted_polynomial_coefficients[1] * x
	 *        + fitted_polynomial_coefficients[2] * std::pow(x, 2)
	 *        + fitted_polynomial_coefficients[3] * std::pow(x, 3)
	 *        + ........................................................
	 *        + fitted_polynomial_coefficients[order] * std::pow(x,
	 *                                                          order);`
	 * at a point (with the argument) x.
	 *
	 * The Vandermonde matrix used for computation is stored in
	 * `vandermonde_mat`, so it should be a matrix of size
	 * 'number of nodes' x 'order of polynomial + 1'.
	 */

	const std::size_t order = N;

	// Create Matrix Placeholder of size
	// 'number of nodes' x 'order of polynomial + 1' ('+1' to store
	// the constant coefficient, i.e. c * (x^0)):
	// Eigen::MatrixXd vandermonde_mat(argument_data.size(), order + 1);
	Eigen::Matrix<
		T, Eigen::Dynamic, 1
	> values_vec = Eigen::Matrix<T, Eigen::Dynamic, 1>::Map(
			&function_value_data.front(), function_value_data.size());
	Eigen::Matrix<T, N + 1, 1> result_vec;

	assert(static_cast<std::size_t>(vandermonde_mat.rows())
		   == argument_data.size());
	assert(argument_data.size() == function_value_data.size());
	assert(argument_data.size() >= order + 1);

	for (std::size_t j = 0 ; j < argument_data.size(); ++ j)
		for (std::size_t k = 0; k < order + 1; ++ k)
			vandermonde_mat(j, k) = std::pow(argument_data[j], k);
			// vandermonde_mat(j, k) = powq(argument_data[j], k);

	// std::cout << vandermonde_mat << "\n";

	// Solve for linear least square fit
	result_vec  = vandermonde_mat.householderQr().solve(values_vec);
	// fitted_polynomial_coefficients.resize(order+1);

	for (std::size_t k = 0; k < order + 1; ++ k) {
		fitted_polynomial_coefficients[k] = result_vec[k];
	}

}
