#pragma once

#include <iostream>

#include "Element.h"
#include "Operands.h"
#include "../../../../Common/2DOutline.h"

namespace MM {

	class CPoint : public CElement {
	public:
		CPoint() : m_x(new CNumber), m_y(new CNumber) { Init(); }
		CPoint(double x, double y) : m_x(new CNumber(x)), m_y(new CNumber(y)) { Init(); }
		CPoint(const CNumber & x, const CNumber & y) : m_x(new CNumber(x)), m_y(new CNumber(y)) { Init(); }
		CPoint(const CUnitSptr & x, const CUnitSptr & y) : m_x(x), m_y(y) { Init(); }
	
		// ---- Доступ к переменным -------------------------------------------------------------------
	
		CNumber X() const;
		CNumber Y() const;
	
		CUnitSptr Result() const;

		int ExportToOutline(C2DOutline * outline) const;
	
		// ---- Ввод-вывод ----------------------------------------------------------------------------
	
		bool Load(CReader & reader);
		void Save(std::ostream & stream) const;

	private:
		void Init() { m_definition = "point"; }
	
		CUnitSptr m_x; //! Указатель на выражение, результат которого обязан являться указателем на экземпляр класса "CNumber".
		CUnitSptr m_y; //! Указатель на выражение, результат которого обязан являться указателем на экземпляр класса "CNumber".
	};

} // namespace MM