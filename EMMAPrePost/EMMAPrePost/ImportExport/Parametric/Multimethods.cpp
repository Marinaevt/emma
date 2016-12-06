#include "StdAfx.h"
#include "xtl/xtl_log.h"
#include "xtl/xtl_functional.h"

#include "Multimethods.h"
#include "DoubleDispatcher.h"
#include "SingleDispatcher.h"

namespace MM {

CUnitSptr operator + (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) {
	CUnitSptr e = CUnitDoubleDispatcher<xtd::plus>::Dispatch(leftOperand, rightOperand);

	if (e == nullptr) {
		xtd::error::report("Сложение элементов " + leftOperand->Definition() + " и " + rightOperand->Definition() + " не зарегистрировано в диспетчере.");
	}

	return e;
}

CUnitSptr operator - (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) {
	CUnitSptr e = CUnitDoubleDispatcher<xtd::minus>::Dispatch(leftOperand, rightOperand);

	if (e == nullptr) {
		xtd::error::report("Вычитание элемента " + leftOperand->Definition() + " из " + rightOperand->Definition() + " не зарегистрировано в диспетчере.");
	}

	return e;
}

CUnitSptr operator * (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) {
	CUnitSptr e = CUnitDoubleDispatcher<xtd::multiplies>::Dispatch(leftOperand, rightOperand);

	if (e == nullptr) {
		xtd::error::report("Умножение элементов " + leftOperand->Definition() + " и " + rightOperand->Definition() + " не зарегистрировано в диспетчере.");
	}

	return e;
}

CUnitSptr operator / (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) {
	CUnitSptr e = CUnitDoubleDispatcher<xtd::divides>::Dispatch(leftOperand, rightOperand);

	if (e == nullptr) {
		xtd::error::report("Деление элемента " + leftOperand->Definition() + " на элемент " + rightOperand->Definition() + " не зарегистрировано в диспетчере.");
	}

	return e;
}

CUnitSptr Power (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) {
	CUnitSptr e = CUnitDoubleDispatcher<xtd::raise>::Dispatch(leftOperand, rightOperand);

	if (e == nullptr) {
		xtd::error::report("Возведение элемента " + leftOperand->Definition() + " в степень элемента " + rightOperand->Definition() + " не зарегистрировано в диспетчере.");
	}

	return e;
}

CUnitSptr operator - (const CUnitSptr & operand) {
	CUnitSptr e = CUnitSingleDispatcher<xtd::negate>::Dispatch(operand);

	if (e == nullptr) {
		xtd::error::report("Отрицание элемента " + operand->Definition() + " не зарегистрировано в диспетчере.");
	}
	
	return e;
}

CUnitSptr Sine (const CUnitSptr & operand) {
	CUnitSptr e = CUnitSingleDispatcher<xtd::sine>::Dispatch(operand);

	if (e == nullptr) {
		xtd::error::report("Взятие синуса элемента " + operand->Definition() + " не зарегистрировано в диспетчере.");
	}
	
	return e;
}

CUnitSptr Cosine (const CUnitSptr & operand) {
	CUnitSptr e = CUnitSingleDispatcher<xtd::cosine>::Dispatch(operand);

	if (e == nullptr) {
		xtd::error::report("Взятие косинуса элемента " + operand->Definition() + " не зарегистрировано в диспетчере.");
	}
	
	return e;
}

CUnitSptr Logarithm (const CUnitSptr & operand) {
	CUnitSptr e = CUnitSingleDispatcher<xtd::logarithm>::Dispatch(operand);

	if (e == nullptr) {
		xtd::error::report("Взятие логарифма элемента " + operand->Definition() + " не зарегистрировано в диспетчере.");
	}
	
	return e;
}

} // namespace MM