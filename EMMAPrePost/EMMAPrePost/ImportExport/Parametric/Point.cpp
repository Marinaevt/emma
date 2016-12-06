#include "StdAfx.h"
#include "Point.h"

namespace MM {

bool CPoint::Load (CReader & reader) {
	parametersContainer parameters;

	if (CElement::Load(reader, parameters) == xtd::good) {
		if (parameters.size() == 2) {
			m_x = parameters[0];
			m_y = parameters[1];
		}
		else if (parameters.size() != 0) {
			bool bUndoLoading = UndoLoading(reader, parameters);
			return bUndoLoading;
		}
	}
	else {
		return xtd::error::report("Элемент \"" + m_definition + "\" описан неправильно.");
	}
	
	bool bFinishLoading = FinishLoading(reader);
	return bFinishLoading;
}

int CPoint::ExportToOutline(C2DOutline * outline) const {
	C2DNode * node = dynamic_cast<C2DNode *>(IO::CreateObject(C2DNODE));

	node->SetPoint(
		Math::C2DPoint(
			X().Value(),
			Y().Value()
			)
		);

	return outline->AddNode(node);
}

void CPoint::Save (std::ostream & stream) const {
	SaveBeginning(stream);
	m_x->SaveAsExpression(stream);
	SaveDelimeter(stream);
	m_y->SaveAsExpression(stream);
	SaveEnding(stream);
}

CUnitSptr CPoint::Result() const {
	return new CPoint(m_x->Result(), m_y->Result());
}

CNumber CPoint::X() const {
	CNumber cNumX = dynamic_cast<CNumber &>(*m_x->Result());
	return cNumX;
}

CNumber CPoint::Y() const {
	CNumber cNumY = dynamic_cast<CNumber &>(*m_y->Result());
	return cNumY;
}

} // namespace MM