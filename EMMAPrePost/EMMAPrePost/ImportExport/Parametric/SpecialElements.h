#pragma once

#include "Element.h"
#include "Operands.h"
#include "Reader.h"

namespace MM {

class CVariable : public CElement {
public:
	CVariable () : m_value(new CNumber) { Init(); }
	CVariable (const CUnitSptr & unit) : m_value(unit) { Init(); }

	CUnitSptr Result () const { return m_value->Result(); }

	int ExportToOutline (C2DOutline * outline) const {
		return dynamic_cast<const CElement &>(*m_value).ExportToOutline(outline);
	}
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const;

private:
	void Init ();

	CUnitSptr m_value;
};

class CVariableLoader : public CElement {
public:
	CVariableLoader () { Init(); }

	CUnitSptr Result () const { return nullptr; }
	
	bool DefinitionIs (const std::string & string) const { return xtd::contains_all_of(m_definition, string); }
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const { xtd::warning::report("Загрузчику переменной нечего сохранять. Это служебный символ"); }
	
private:
	void Init ();
};

class CElementEnd : public CElement {
public:
	CElementEnd () { Init(); }

	CUnitSptr Result () const { return nullptr; }
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const { xtd::warning::report("Символу окончания элемента нечего сохранять. Это служебный символ"); }

private:
	void Init () { m_definition = ";"; }
	//void SetDefaultValues () { xtd::warning::report("Функция пока не работает."); } //не используется
};

} // namespace MM