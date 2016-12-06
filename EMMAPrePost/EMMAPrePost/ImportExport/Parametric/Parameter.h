#pragma once

#include "Element.h"
#include "Operands.h"

namespace MM {

//!		Параметр метамодели.
/*!
		Тот самый элемент, который является характеристикой конкретной модели, и варьируя который можно
	автоматически получать необходимую геометрию.
*/
class CParameter : public CElement {
public:
	CParameter () : m_value(new CNumber) { Init(); }
	CParameter (const CNumber & value) : m_value(new CNumber(value)) { Init(); }
	CParameter (const CUnitSptr & unit) : m_value(unit) { Init(); }
	
	void SetValue (const CNumber & value) { m_value = new CNumber(value); }
	
	// ---- Доступ к переменным -------------------------------------------------------------------
	
	CNumber Value () const { return dynamic_cast<const CNumber &>(*m_value->Result()); }
	double DoubleValue () const { return Value().Value(); }
	
	CUnitSptr Result () const { return new CNumber(Value()); }
	
	// ---- Ввод-вывод ----------------------------------------------------------------------------
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const;

private:
	void Init () { m_definition = "parameter"; }
	
	CUnitSptr m_value; //! Указатель на выражение, результат которого обязан являться указателем на экземпляр класса "CNumber".
};

} // namespace MM