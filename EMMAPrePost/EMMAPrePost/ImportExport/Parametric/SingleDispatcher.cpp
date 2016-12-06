#include "StdAfx.h"
#include "SingleDispatcher.h"

//     Без указания инстанцировать этот шаблон специализированные инициализации будут выдавать
// ошибку доступа к памяти из-за того, что умный указатель на унарную функцию не успевает
// инициализироваться.
template class xtd::smart_ptr<xtd::unary_function<xtd::smart_ptr<MM::CUnit> > >;

namespace MM {

// Создание статической переменной.
template <template <typename ArgumentType, typename ResultType> class UnaryPredicate>
CUnitSingleDispatcher<UnaryPredicate> CUnitSingleDispatcher<UnaryPredicate>::instance;

// ---- Специализация инициализаций ---------------------------------------------------------------

//!		Инициализация диспетчера отрицания.
template <>
inline void CUnitSingleDispatcher<xtd::negate>::Init () {
	Put<CNumber>();             // Отрицание числа.
	Put<CParameter, CNumber>(); // Отрицание параметра, результат — число.
//	Put<CPoint>();
}

template <>
inline void CUnitSingleDispatcher<xtd::sine>::Init () {
	Put<CNumber>();
	Put<CParameter, CNumber>(); // Взятие синуса параметра, результат — число.
}

template <>
inline void CUnitSingleDispatcher<xtd::cosine>::Init () {
	Put<CNumber>();
	Put<CParameter, CNumber>(); // Взятие косинуса параметра, результат — число.
}

template <>
inline void CUnitSingleDispatcher<xtd::logarithm>::Init () {
	Put<CNumber>();
	Put<CParameter, CNumber>(); // Взятие логарифма параметра, результат — число.
}

} // namespace MM