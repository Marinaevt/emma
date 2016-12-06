#include "StdAfx.h"
#include "Contour.h"

namespace MM {
	
CUnitSptr CContour::Result () const {
	CContour contour;
	
	for (std::vector<CUnitSptr>::const_iterator i = m_curves.begin(); i != m_curves.end(); ++i) {
		CUnitSptr pUnitSptr = *i;
		contour.AddCurve(pUnitSptr->Result());
	}
	
	return new CContour(contour);
	
}

int CContour::ExportToOutline(C2DOutline * outline) const {
	C2DContour * contour = dynamic_cast<C2DContour *>(IO::CreateObject(C2DCONTOUR));
	contour->SetParent(outline);

	bool bAdd = false;
	for (parametersContainer::const_iterator i = m_curves.begin(); i != m_curves.end(); ++i) {
		int curve_index = dynamic_cast<const CElement  &>(**i).ExportToOutline(outline);
		if (curve_index < 0) {
			return -1;
		}
		bAdd = contour->AddCurve(curve_index);
	}

	int nAddContour = outline->AddContour(contour);
	return nAddContour;
}

bool CContour::Load (CReader & reader) {
	parametersContainer parameters;

	if (CElement::Load(reader, parameters) == xtd::good) {
		if (parameters.size() != 0) {
			m_curves = parameters;
		}
	}
	else {
		return xtd::error::report("Элемент \"" + m_definition + "\" описан неправильно.");
	}
	
	return FinishLoading(reader);
}

void CContour::Save (std::ostream & stream) const {
	SaveBeginning(stream);
	
	if (!m_curves.empty()) {
		parametersContainer::const_iterator i = m_curves.begin();
		
		(*i++)->SaveAsExpression(stream);
	
		while (i != m_curves.end()) {
			SaveDelimeter(stream);
			(*i++)->SaveAsExpression(stream);
		}
	}
	
	SaveEnding(stream);
}

void CContour::Init () {
	m_definition = "contour";
	m_name.clear();
}

} // namespace MM