#include "StdAfx.h"
#include "CContour.h"

namespace DXF {

bool CContour::LoadFromDXF (CReader & reader) {
	if (reader.SeekForEntities() == true) {
		while (reader.NextToken().IsSectionEnd() != true) {
			CCurve * curve = AllocateCurveForToken(reader.GetToken());
			
			if (curve != nullptr) {
				if (curve->LoadFromDXF(reader) == true)
					AddCurve(curve);
				else {
					logMessage("ошибка: не удалось загрузить примитив");
					delete curve;
					return false;
				}
			} else {
				logMessage("предупреждение: примитив не опознан и будет проигнорирован");
				reader.SkipEntity();
				continue;
			}
		}
	} else {
		logMessage("ошибка: секция примитивов не найдена");
		return false;
	}
	
	return true;
}

void CContour::RegisterCurvesInGraph (CGraph & graph) {
	for (CCurvePtrArray::iterator i = m_curves.begin(); i != m_curves.end(); ++i) {
		(*i)->RegisterPointsInGraph(graph);
		
		*i = graph.AddCurve(*i);
	}
}

int CContour::ExportToOutline(C2DOutline * outline) const {
	C2DContour * contour = (C2DContour *)IO::CreateObject(C2DCONTOUR);
	contour->SetParent(outline);

	for (CCurvePtrArray::const_iterator i = m_curves.begin(); i != m_curves.end(); ++i) {
		CCurve *pCurve = *i;
		int curve_index = pCurve->ExportToOutline(outline);

		contour->AddCurve(curve_index);
	}

	return outline->AddContour(contour);
}

bool CContour::operator == (const CContour & contour) const {
	if (m_curves == contour.m_curves) return true;
	// Возможно, нужно проходить по всему массиву и проверять соответствие не только указателей,
	// но и значений объектов.
	
	return false;
}

} // namespace DXF