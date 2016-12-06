#pragma once

#include <vector>
#include "Element.h"
#include "Curves.h"
#include "../../../../Common/2DOutline.h"

namespace MM {

	class CContour : public CElement {
	public:
		CContour() { Init(); }
		CContour(const std::vector<CUnitSptr> & curves) : m_curves(curves) { Init(); }
	
		void AddCurve(const CUnitSptr & curve) { m_curves.push_back(curve); }
	
		// ---- Доступ к переменным -------------------------------------------------------------------
	
		CUnitSptr Result() const;
	
		int ExportToOutline(C2DOutline * outline) const;
	
		// ---- Ввод-вывод ----------------------------------------------------------------------------
	
		bool Load(CReader & reader);
		void Save(std::ostream & stream) const;

	private:
		void Init();
	
		parametersContainer m_curves;
	};

} // namespace MM