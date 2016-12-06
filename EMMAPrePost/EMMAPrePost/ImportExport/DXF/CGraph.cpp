#include "StdAfx.h"
#include "CGraph.h"

namespace DXF {

CGraph::~CGraph () {
	deleteVector(m_points);
	deleteVector(m_curves);
}

bool CGraph::LoadFromDXF (const char * path) {
	CReader reader(path);
	
	return LoadFromDXF(reader);
}

bool CGraph::LoadFromDXF (CReader & reader) {
	CContour contour;
	
	if (contour.LoadFromDXF(reader) == true)
		AddContour(contour);
	else {
		logMessage("предупреждение: не загружено ни одного контура");
		return false;
	}
	
	return true;
}

void CGraph::AddContour (CContour & contour) {
	contour.RegisterCurvesInGraph(*this);
	
	for (CContourArray::iterator i = m_contours.begin(); i != m_contours.end(); ++i) {
		if (*i == contour) return;
	}
	
	m_contours.push_back(contour);
}

CCurve * CGraph::AddCurve (CCurve * curve) {
	for (CCurvePtrArray::iterator i = m_curves.begin(); i != m_curves.end(); ++i) {
		CCurve *pCurve = *i;
		if (pCurve->Compare(curve) == true) {
			delete curve;
			return pCurve;
		}
	}

	m_curves.push_back(curve);
	
	return curve;
}

CPoint * CGraph::AddPoint (CPoint * point) {
	for (CPointPtrArray::iterator i = m_points.begin(); i != m_points.end(); ++i) {	
		CPoint *pPoint = *i;
		if (*pPoint == *point) {
			delete point;
			return pPoint;
		}
	}
	
	m_points.push_back(point);

	return point;
}

C2DOutline * CGraph::MakeOutline () const {
	C2DOutline * outline = dynamic_cast<C2DOutline *>(IO::CreateObject(C2DOUTLINE));

	for (CContourArray::const_iterator i = m_contours.begin(); i != m_contours.end(); ++i)
		i->ExportToOutline(outline);
		
	return outline;
} //раскомментировать

} // namespace DXF