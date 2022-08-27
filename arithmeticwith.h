#ifndef ARITHMETICWITH_H
#define ARITHMETICWITH_H

#include <concepts>


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


using numeric_val = double;


#endif // ARITHMETICWITH_H
