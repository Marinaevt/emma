#include "StdAfx.h"
#include "DoubleDispatcher.h"

//     Без указания инстанцировать этот шаблон специализированные инициализации будут выдавать
// ошибку доступа к памяти из-за того, что умный указатель на бинарную функцию не успевает
// инициализироваться.
template class xtd::smart_ptr<xtd::binary_function<xtd::smart_ptr<MM::CUnit> > >;

namespace MM {

// Создание статической переменной.
template <template <typename FirstArgumentType, typename SecondArgumentType, typename PredicateResultType> class BinaryPredicate>
CUnitDoubleDispatcher<BinaryPredicate> CUnitDoubleDispatcher<BinaryPredicate>::instance;

// ---- Специализации инициализаций ---------------------------------------------------------------

//!		Инициализация диспетчера сложения.
template <>
inline void CUnitDoubleDispatcher<xtd::plus>::Init () {
	Put<CNumber>();                         // Сумма двух чисел, результат — число.
	Put<CParameter, CParameter, CNumber>(); // Сумма двух параметров, результат — число.
	Put<CPoint>();                          // Сумма двух точек, результат — точка.
	Put<CParameter, CNumber, CNumber>();    // Сумма параметра и числа, результат — число.
	Put<CNumber, CParameter, CNumber>();    // Сумма числа и параметра, результат — число.
}

//!		Инициализация диспетчера вычитания.
template <>
inline void CUnitDoubleDispatcher<xtd::minus>::Init () {
	Put<CNumber>();                         // Разность двух чисел, результат — число.
	Put<CParameter, CParameter, CNumber>(); // Разность двух параметров, результат — число.
	Put<CPoint>();                          // Разность двух точек, результат — точка.
	Put<CParameter, CNumber, CNumber>();    // Разность параметра и числа, результат — число.
	Put<CNumber, CParameter, CNumber>();    // Разность числа и параметром, результат — число.
}

//!		Инициализация диспетчера умножения.
template <>
inline void CUnitDoubleDispatcher<xtd::multiplies>::Init () {
	Put<CNumber>();                         // Произведение двух чисел, результат — число.
	Put<CParameter, CParameter, CNumber>(); // Произведение двух параметров, результат — число.
	Put<CParameter, CNumber, CNumber>();    // Произведение параметра и числа, результат — число.
	Put<CNumber, CParameter, CNumber>();    // Произведение числа и параметра, результат — число.
	Put<CPoint, CNumber, CPoint>();         // Произведение точки и числа, результат — точка.
	Put<CNumber, CPoint, CPoint>();         // Произведение числа и точки, результат — точка.
}

//!		Инициализация диспетчера деления.
template <>
inline void CUnitDoubleDispatcher<xtd::divides>::Init () {
	Put<CNumber>();                         // Частное двух чисел, результат — число.
	Put<CParameter, CParameter, CNumber>(); // Частное двух параметров, результат — число.
	Put<CParameter, CNumber, CNumber>();    // Частное параметра и числа, результат — число.
	Put<CNumber, CParameter, CNumber>();    // Частное числа и параметра, результат — число.
	Put<CPoint, CNumber, CPoint>();         // Частное точки и числа, результат — точка.
}

//!		Инициализация диспетчера возведения в степень.
template <>
inline void CUnitDoubleDispatcher<xtd::raise>::Init() {
	Put<CNumber>();							// Возведение числа в степень другого числа.
	Put<CParameter, CParameter, CNumber>(); // Возведение параметра в степень другого параметра, результат — число.
	Put<CParameter, CNumber, CNumber>();    // Возведение параметра в степень числа, результат — число.
	Put<CNumber, CParameter, CNumber>();    // Возведение числа в степень параметра, результат — число.
}

} // namespace MM