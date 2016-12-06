#pragma once

#include "Unit.h"
#include "Element.h"

namespace MM {

// ---- Унарные функторы сравнения ----------------------------------------------------------------

class FUnitDefinitionIs {
public:
	explicit FUnitDefinitionIs (const std::string & characters) : m_characters(characters) {}

	bool operator () (const CUnit & unit) const { return unit.DefinitionIs(m_characters); }
	
private:
	std::string m_characters;
};

class FElementNameIs {
public:
	explicit FElementNameIs (const std::string & characters) : m_characters(characters) {}

	bool operator () (const CElement & element) const { return element.NameIs(m_characters); }
	
private:
	std::string m_characters;
};

// ---- Бинарные функторы сравнения ---------------------------------------------------------------

struct FUnitDefinitionIsLess {
	bool operator () (const CUnit & unit1, const CUnit & unit2) const { return unit1.Definition() < unit2.Definition(); }
};

struct FUnitSptrDefinitionIsLess {
	bool operator () (const CUnitSptr & unit1, const CUnitSptr & unit2) const { return unit1->Definition() < unit2->Definition(); }
};

//struct FElementNameIsLess : public FElementBinaryPredicate {
//	bool operator () (const CElement & element1, const CElement & element2) const { return element1.Name() < element2.Name(); }
//};

// ---- Сравнения наоборот ------------------------------------------------------------------------

class FStringIsDefinitionOf {
public:
	explicit FStringIsDefinitionOf (const CUnit & unit) : m_unit(unit) {}
	
	bool operator () (const std::string & string) const { return m_unit.DefinitionIs(string); }

private:
	const CUnit & m_unit;
};

} // namespace MM