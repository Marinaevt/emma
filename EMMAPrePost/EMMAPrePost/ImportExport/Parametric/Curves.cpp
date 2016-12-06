#include "StdAfx.h"
#include "Curves.h"

namespace MM {

int CLine::ExportToOutline(C2DOutline * outline) const {
		C2DCurve * line = dynamic_cast<C2DCurve *>(IO::CreateObject(C2DLINESEG));

		int nStartNode = StartPoint().ExportToOutline(outline);
		int nEndNode = EndPoint().ExportToOutline(outline);

		line->SetNodes(nStartNode, nEndNode);

		return outline->AddCurve(line);
	}


bool CLine::Load (CReader & reader) {
	parametersContainer parameters;

	if (CElement::Load(reader, parameters) == xtd::good) {
		if (parameters.size() == 2) {
			m_startPoint = parameters[0];
			m_endPoint = parameters[1];
		}
		else if (parameters.size() != 0) {
			return UndoLoading(reader, parameters);
		}
	}
	else {
		return xtd::error::report("Элемент \"" + m_definition + "\" описан неправильно.");
	}
	
	return FinishLoading(reader);
}

void CLine::Save (std::ostream & stream) const {
	SaveBeginning(stream);
	m_startPoint->SaveAsExpression(stream);
	SaveDelimeter(stream);
	m_endPoint->SaveAsExpression(stream);
	SaveEnding(stream);
}

/*	--------------------
	CCircle - окружности
*/

bool CCircle::Load (CReader & reader) {
	parametersContainer parameters;

	if (CElement::Load(reader, parameters) == xtd::good) {
		if (parameters.size() == 2) {
			m_center = parameters[0];
			m_radius = parameters[1];
		}
		else if (parameters.size() != 0) {
			return UndoLoading(reader, parameters);
		}
	}
	else {
		return xtd::error::report("Элемент \"" + m_definition + "\" описан неправильно.");
	}
	
	return FinishLoading(reader);
}

void CCircle::Save (std::ostream & stream) const {
	SaveBeginning(stream);
	m_center->SaveAsExpression(stream);
	SaveDelimeter(stream);
	m_radius->SaveAsExpression(stream);
	SaveEnding(stream);
}

/*
int CCircle::ExportToOutline (C2DOutline * outline) const {
	//
	//TODO: дописать экспорт (когда-нибудь)
	//
	return 0;
}
//*/

/*	-----------
	CArc - дуги
*/

bool CArc::Load (CReader & reader) {
	parametersContainer parameters;
	
	if (CElement::Load(reader, parameters) == xtd::good) {
		if (parameters.size() == 3) {
			m_startPoint = parameters[0];
			m_endPoint = parameters[1];
			m_radius = parameters[2];
		} else if (parameters.size() != 0) return UndoLoading(reader, parameters);
	}
	else {
		return xtd::error::report("Элемент \"" + m_definition + "\" описан неправильно.");
	}
	
	return FinishLoading(reader);
}

void CArc::Save (std::ostream & stream) const {
	SaveBeginning(stream);
	m_startPoint->SaveAsExpression(stream);
	SaveDelimeter(stream);
	m_endPoint->SaveAsExpression(stream);
	SaveDelimeter(stream);
	m_radius->SaveAsExpression(stream);
	SaveEnding(stream);
}

int CArc::ExportToOutline(C2DOutline * outline) const {
	/*
	C2DCircleArc *arc = dynamic_cast<C2DCircleArc *>(IO::CreateObject(C2DCIRCARC));

	arc->SetNodes(
		StartPoint().ExportToOutline(outline),
		EndPoint().ExportToOutline(outline)
		);

	arc->SetRadius(Radius().Value());
	arc->SetCW(false);
	*/
	int n1 = StartPoint().ExportToOutline(outline);
	int n2 = EndPoint().ExportToOutline(outline);
	DBL dRad = Radius().Value() + EPS;
	DBL dRad2 = outline->GetNode(n1)->GetPoint().Len(outline->GetNode(n2)->GetPoint()) / 2.0 + EPS;

	//! предотвращаем ошибки записи\окргулений, считанный радиус не может быть меньше стандартного, больше - может.
	//! есть зависимость от точности координат узлов
	if (dRad - dRad2 < EPS) {
		dRad = dRad2;
	}

	C2DCircleArc *arc = new C2DCircleArc(n1,n2,dRad);

	int nAdd = outline->AddCurve(arc);

	return nAdd;
}

} // namespace MM