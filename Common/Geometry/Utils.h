#pragma once
#include <functional>
#include <iterator>

namespace FemLibrary {

	class Utils
	{
	public:
		template<class T> static T NextInd(T i, T n) { return i == n - 1 ? 0 : i + 1; }
		template<class T> static T PrevInd(T i, T n) { return i == 0 ? n - 1 : i - 1; }

		static size_t hash_combine(size_t seed, size_t v) {
			// http://www.boost.org/doc/libs/1_53_0/doc/html/hash/reference.html#boost.hash_combine
			return seed ^ (v + 0x9e3779b9 + (seed << 6) + (seed >> 2));
		}

		//template <typename T, typename R>
		//static std::vector<T> Convert(const std::vector<R>& src, std::function func)
		//{
		//	std::vector<T> result;
		//	result.reserve(src.size());
		//	std::transform(src.begin(), src.end(), std::back_inserter(result), func);
		//	return result;
		//}
	};
}
