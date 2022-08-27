#include <valarray>

#include <functional>


template <typename T>
std::valarray<T> operator + (
		const std::valarray<T>& arr,
		const std::ranges::common_range auto& some_range) {
	std::valarray<T> res(arr.size());
	std::transform(
				std::ranges::begin(arr), std::ranges::end(arr),
				std::ranges::begin(some_range),
				std::ranges::begin(res), std::plus<>{});

	return res;
}


template <typename T>
std::valarray<T> operator - (
		const std::valarray<T>& arr,
		const std::ranges::common_range auto& some_range) {
	std::valarray<T> res(arr.size());
	std::transform(
				std::ranges::begin(arr), std::ranges::end(arr),
				std::ranges::begin(some_range),
				std::ranges::begin(res), std::minus<>{});

	return res;
}


template <typename T>
std::valarray<T>& operator += (
		std::valarray<T>& arr,
		const std::ranges::common_range auto& some_range) {
	std::transform(
				std::ranges::begin(arr), std::ranges::end(arr),
				std::ranges::begin(some_range),
				std::ranges::begin(arr), std::plus<>{});

	return arr;
}


template <typename T>
std::valarray<T>& operator -= (
		std::valarray<T>& arr,
		const std::ranges::common_range auto& some_range) {
	std::transform(
				std::ranges::begin(arr), std::ranges::end(arr),
				std::ranges::begin(some_range),
				std::ranges::begin(arr), std::minus<>{});

	return arr;
}
