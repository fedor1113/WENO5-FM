#ifndef ARITHMETICWITH_H
#define ARITHMETICWITH_H

// #include <cmath>
#include <concepts>

// #include <quadmath.h>


//template<class N>
//concept SignedNumber = std::is_signed_v<N>
//	&& std::is_arithmetic<N>
//	&& std::is_floating_point<N>;


template<typename T, typename N>
concept AddableWith = requires (T x, N y) {
	x + x; x + y; y + x; x += y;
};


template<typename T, typename N>
concept SubtractableWith = requires (T x, N y) {
	x - x; x - y; y - x; x -= y;
};


template<typename T, typename N>
concept MultipliableWith = requires (T x, N y) {
	x * x; x * y; y * x; x *= y;
};


template<typename T, typename N>
concept DivisibleWith = requires (T x, N y) {
	x / x; x / y; y / x; x *= y;
};


template<typename T, typename N>
concept ArithmeticWith = AddableWith<T, N>
		&& SubtractableWith<T, N>
		&& MultipliableWith<T, N>
		&& DivisibleWith<T, N>;


using numeric_val = long double;
// using numeric_val = double_t;
// using numeric_val = __float128;


//namespace std {
////numeric_val cos(numeric_val qn) {
////	return cosq(qn);
////}


////numeric_val sin(numeric_val qn) {
////	return sinq(qn);
////}


//numeric_val sqrt(numeric_val qn) {
//	return sqrtq(qn);
//}


//numeric_val pow(numeric_val qn1, numeric_val qn2) {
//	return powq(qn1, qn2);
//}


////numeric_val exp(numeric_val qn) {
////	return expq(qn);
////}
//}


#endif // ARITHMETICWITH_H
