#include "StdAfx.h"
#include "utility.h"
#include "CPoint.h"

namespace DXF {

bool CPoint::LoadFromDXF (CReader & reader) {
	CToken token;

	token = reader.GetToken();
	if (isXCoordinate(token.m_code) == true) x = atof(token.m_command.c_str());
	else {
		logMessage("ошибка: отсутствует x-координата точки");
		return false;
	}
	
	token = reader.GetToken();
	if (isYCoordinate(token.m_code) == true) y = atof(token.m_command.c_str());
	else {
		logMessage("ошибка: отсутствует y-координата точки");
		return false;
	}
	
	token = reader.GetToken();
	if (isZCoordinate(token.m_code) == true) z = atof(token.m_command.c_str());
	else {
		logMessage("ошибка: отсутствует z-координата точки");
		return false;
	}
	
	return true;
}

bool CPoint::operator == (const CPoint & point) {
	if (
		fabs(x - point.x) < EPS &&
		fabs(y - point.y) < EPS &&
		fabs(z - point.z) < EPS
	) return true;
	else return false;
}

CPoint makeCPoint (double x, double y, double z) {
	CPoint point(x, y, z);
	
	return point;
}

}