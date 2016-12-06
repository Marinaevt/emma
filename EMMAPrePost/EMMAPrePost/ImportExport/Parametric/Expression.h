#pragma once

#include "xtl/xtl_memory.h"

#include "Unit.h"

namespace MM {

//!		Арифметическое выражение.
/*!
		Базовый класс для всех арифметических выражений, операторов и простых операнд: чисел и
	символьных констант.
*/
class CExpression : public CUnit {
public:
	// ---- Определения ---------------------------------------------------------------------------
	
	//!		Приоритет.
	/*!
			Класс, предназначенный для удобной работы с приоритетами.
	*/
	class CPriority {
	public:
		//!		Перечисление приоритетов.
		/*!
				Здесь указаны все приоритеты в порядке возрастания. Для того, чтобы изменить обработку приоритетов,
			достаточно поменять местами значения в этом списке.
		*/
		typedef enum _EPriority {
			Lowest,						//! Нижняя граница приоритетов.
			
			AdditiveOperator,			//! Приоритет аддитивных операторов (умножения, сложения, вычитания).
			MultiplicativeOperator,		//! Приоритет мультипликативных операторов (умножения, деления, деления с остатком и т.д.).
			PrefixOperator,				//! Приоритет префиксных унарных операторов (изменение знака и т.п.).
			Power,						//! Приоритет оператора возведения в степень.
			Function,					//! Приоритет функций (синуса, косинуса и т.п.).
			PostfixOperator,			//! Приоритет постфиксных унарных операторов (факториала и т.п.).
			Brackets,					//! Приоритет скобок.
			Number,						//! Приоритет всех чисел, переменных и символьных констант (e, pi и т.д.).
			
			Highest						//! Верхняя граница приоритетов.
		} EPriority;
		
		CPriority () {
			m_value = Lowest;	//по-умолчанию приоритет чисел
		}
		CPriority (const EPriority & value) : m_value(value) {}
		
		bool operator == (const CPriority &     p) const { return m_value == p.m_value; }
		bool operator == (const EPriority & value) const { return m_value ==     value; }
		
		bool operator  < (const CPriority & p) const { return    m_value < p.m_value; }
		bool operator <= (const CPriority & p) const { return !(p.m_value < m_value); }
		
	private:
		EPriority m_value;
	};
	
	// ---- Модификация ---------------------------------------------------------------------------
	
	void SetPriority (const CPriority & priority) { m_priority = priority; }
	
	// ---- Доступ к переменным -------------------------------------------------------------------
	
	CPriority Priority () const { return m_priority; }
	
	// ---- Ввод-вывод ----------------------------------------------------------------------------
	
	virtual void SaveAsExpression (std::ostream & stream) const { Save(stream); }
	
protected:
	CPriority m_priority;
};

} // namespace MM