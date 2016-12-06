#ifndef _XTD_FUNCTIONAL_H
#define _XTD_FUNCTIONAL_H

namespace xtd {

// ----	Базовые классы ----------------------------------------------------------------------------
/*
	Отличие от аналогичных классов из STL состоит в существовании значений по-умолчанию для второго и третьего параметров шаблона.
*/

template <typename first_argument_type, typename second_argument_type = first_argument_type, typename result_type = first_argument_type>
struct binary_function {
	virtual ~binary_function () {}

	virtual result_type operator () (const first_argument_type & x, const second_argument_type & y) const = 0;
};

template <typename argument_type, typename result_type = argument_type>
struct unary_function {
	virtual ~unary_function () {}

	virtual result_type operator () (const argument_type & x) const = 0;
};

// ---- Арифметические операции -------------------------------------------------------------------
/*
	Отличие от аналогичных классов из STL состоит в возможности оперировать с произвольными типами и получать произвольный результат.
*/

template <typename first_argument_type, typename second_argument_type = first_argument_type, typename result_type = first_argument_type>
struct plus : binary_function <first_argument_type, second_argument_type, result_type> {
	result_type operator () (const first_argument_type & x, const second_argument_type & y) const { return x + y; }
};

template <typename first_argument_type, typename second_argument_type = first_argument_type, typename result_type = first_argument_type>
struct minus : binary_function <first_argument_type, second_argument_type, result_type> {
	result_type operator () (const first_argument_type & x, const second_argument_type & y) const { return x - y; }
};

template <typename first_argument_type, typename second_argument_type = first_argument_type, typename result_type = first_argument_type>
struct multiplies : binary_function <first_argument_type, second_argument_type, result_type> {
	result_type operator () (const first_argument_type & x, const second_argument_type & y) const { return x * y; }
};

template <typename first_argument_type, typename second_argument_type = first_argument_type, typename result_type = first_argument_type>
struct divides : binary_function <first_argument_type, second_argument_type, result_type> {
	result_type operator () (const first_argument_type & x, const second_argument_type & y) const { return x / y; }
};

template <typename first_argument_type, typename second_argument_type = first_argument_type, typename result_type = first_argument_type>
struct raise : binary_function <first_argument_type, second_argument_type, result_type> {
	result_type operator () (const first_argument_type & x, const second_argument_type & y) const { return pow(x, y); }
};

template <typename argument_type, typename result_type = argument_type>
struct negate : unary_function <argument_type, result_type> {
	result_type operator () (const argument_type & x) const { return -x; }
};

template <typename argument_type, typename result_type = argument_type>
struct sine : unary_function <argument_type, result_type> {
	result_type operator () (const argument_type & x) const { return sin(x); }
};

template <typename argument_type, typename result_type = argument_type>
struct cosine : unary_function <argument_type, result_type> {
	result_type operator () (const argument_type & x) const { return cos(x); }
};

template <typename argument_type, typename result_type = argument_type>
struct logarithm : unary_function <argument_type, result_type> {
	result_type operator () (const argument_type & x) const { return log(x); }
};

// ---- Операции сравнения ------------------------------------------------------------------------
/*
	Отличие от аналогичных классов из STL состоит в возможности сравнивать разные типы.
*/

template <typename first_argument_type, typename second_argument_type = first_argument_type>
struct less : binary_function <first_argument_type, second_argument_type, bool> {
	bool operator () (const first_argument_type & x, const second_argument_type & y) const { return x < y; }
};

}

#endif