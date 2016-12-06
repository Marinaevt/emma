#include "StdAfx.h"
#include "CCurve.h"

namespace DXF {

bool CLine::LoadFromDXF (CReader & reader) {
	
	DXF::GroupCode cline_code = reader.NextToken().m_code;
	
	while (cline_code != ELEMENT_START) {
		if (cline_code == FIRST_X_COORDINATE) {
			CPoint * point = new CPoint;
			
			if (point->LoadFromDXF(reader) == true) {
				m_startPoint = point;
				continue;
			} else {
				logMessage("ошибка: не удалось загрузить начальную точку отрезка");
				delete point;
				return false;
			}
		} else if (cline_code == SECOND_X_COORDINATE) {
			CPoint * point = new CPoint;
			
			if (point->LoadFromDXF(reader) == true) {
				m_endPoint = point;
				continue;
			} else {
				logMessage("ошибка: не удалось загрузить конечную точку отрезка");
				delete point;
				return false;
			}
		} else {
			logMessage("предупреждение: параметр не опознан и будет проигнорирован");
		}
		
		reader.GetToken();
	}
	
	return true;
}

void CLine::RegisterPointsInGraph (CGraph & graph) {
	m_startPoint = graph.AddPoint(m_startPoint);
	m_endPoint = graph.AddPoint(m_endPoint);
}

bool CLine::Compare (CCurve * curve) const {
	if (curve->IsLine()) {
		return (*this == *dynamic_cast<CLine *>(curve));
	} else return false;
}

bool CLine::operator == (const CLine & line) const {
	if (
		*m_startPoint == *line.m_startPoint &&
		*m_endPoint == *line.m_endPoint
	) return true;
	
	return false;
}

// ------------------------------------------------------------------------------------------------

bool CCircle::LoadFromDXF (CReader & reader) {
	DXF::CToken tNext = reader.NextToken();
	while (tNext.m_code != ELEMENT_START) {
		if (tNext.m_code == FIRST_X_COORDINATE) {
			CPoint * point = new CPoint;
			
			if (point->LoadFromDXF(reader) == true) {
				m_centerPoint = point;
				continue;
			} else {
				logMessage("ошибка: не удалось загрузить точку центра окружности");
				delete point;
				return false;
			}
		} else if (tNext.m_code == RADIUS) {
			m_radius = atof(tNext.m_command.c_str());
		} else {
			logMessage("предупреждение: параметр не опознан и будет проигнорирован");
		}
		
		reader.GetToken();
	}
	
	return true;
}
	
void CCircle::RegisterPointsInGraph (CGraph & graph) {
	m_centerPoint = graph.AddPoint(m_centerPoint);
}

bool CCircle::Compare (CCurve * curve) const {
	if (curve->IsCircle()) {
		return (*this == *dynamic_cast<CCircle *>(curve));
	} else return false;
}

bool CCircle::operator == (const CCircle & circle) const {
	if (
		*m_centerPoint == *circle.m_centerPoint &&
		fabs(m_radius - circle.m_radius) < EPS
	) return true;
	
	return false;
}

// ------------------------------------------------------------------------------------------------

bool CArc::LoadFromDXF (CReader & reader) {
	
	DXF::CToken nextToken = reader.NextToken();
	while (nextToken.m_code != ELEMENT_START) {
		if (nextToken.m_code == FIRST_X_COORDINATE) {
			CPoint * point = new CPoint;
			
			if (point->LoadFromDXF(reader) == true) {
				m_centerPoint = point;
				continue;
			} else {
				logMessage("ошибка: не удалось загрузить точку центра окружности");
				delete point;
				return false;
			}
		} else if (nextToken.m_code == RADIUS) {
			m_radius = atof(nextToken.m_command.c_str());
		} else if (nextToken.m_code == START_ANGLE) {
			m_startAngle = atof(nextToken.m_command.c_str());
		} else if (nextToken.m_code == END_ANGLE) {
			m_endAngle = atof(nextToken.m_command.c_str());
		} else {
			logMessage("предупреждение: параметр не опознан и будет проигнорирован");
		}
		
		nextToken = reader.GetToken();
	}
	
	return true;
}

void CArc::RegisterPointsInGraph (CGraph & graph) {
	m_centerPoint = graph.AddPoint(m_centerPoint);
}

bool CArc::Compare (CCurve * curve) const {
	if (curve->IsArc()) {
		return (*this == *dynamic_cast<CArc *>(curve));
	} else return false;
}

bool CArc::operator == (const CArc & arc) const {
	if (
		*m_centerPoint == *arc.m_centerPoint &&
		fabs(m_radius - arc.m_radius) < EPS &&
		fabs(m_startAngle - arc.m_startAngle) < EPS &&
		fabs(m_endAngle - arc.m_endAngle) < EPS
	) return true;
	
	return false;
}

// ------------------------------------------------------------------------------------------------

bool CVertex::LoadFromDXF (CReader & reader) {

	DXF::CToken nextToken = reader.NextToken();

	while (nextToken.m_code != ELEMENT_START) {
		if (nextToken.m_code == FIRST_X_COORDINATE) {
			CPoint * point = new CPoint;
			
			if (point->LoadFromDXF(reader) == true) {
				m_point = point;
				continue;
			} else {
				logMessage("ошибка: не удалось загрузить координаты вершины полилинии");
				delete point;
				return false;
			}
		} else if (nextToken.m_code == START_WIDTH) {
			m_startWidth = atof(nextToken.m_command.c_str());
		} else if (nextToken.m_code == END_WIDTH) {
			m_endWidth = atof(nextToken.m_command.c_str());
		} else if (nextToken.m_code == CURVATURE) {
			m_curvature = atof(nextToken.m_command.c_str());
		} else if (nextToken.m_code == TANGENT_DIRECTION) {
			m_tangentDirection = atof(nextToken.m_command.c_str());
		} else if (nextToken.m_code == FLAG70) {
			m_flag = atof(nextToken.m_command.c_str());
		} else {
			logMessage("предупреждение: параметр не опознан и будет проигнорирован");
		}
		
		nextToken = reader.GetToken();
	}
	
	return true;
}

void CVertex::RegisterPointInGraph (CGraph & graph) {
	m_point = graph.AddPoint(m_point);
}

bool CVertex::operator == (const CVertex & vertex) const {
	if (
		*m_point == *vertex.m_point &&
		fabs(m_startWidth - vertex.m_startWidth) < EPS &&
		fabs(m_endWidth - vertex.m_endWidth) < EPS &&
		fabs(m_curvature - vertex.m_curvature) < EPS &&
		fabs(m_tangentDirection - vertex.m_tangentDirection) < EPS &&
		fabs(m_flag - vertex.m_flag) < EPS
	) return true;

	return false;
}

bool CSequenceEnd::LoadFromDXF (CReader & reader) {
	while (reader.NextToken().m_code != ELEMENT_START) {
		reader.GetToken();
	}
	
	return true;
}

bool CPolyline::LoadFromDXF (CReader & reader) {
	// Чтение параметров полилинии.
	while (reader.NextToken().m_code != ELEMENT_START) {
		CToken token = reader.GetToken();
		
		if (token.m_code == START_WIDTH) {
			m_startWidth = atof(token.m_command.c_str());
		} else if (token.m_code == END_WIDTH) {
			m_endWidth = atof(token.m_command.c_str());
		} else {
			logMessage("предупреждение: параметр не опознан и будет проигнорирован");
		}
	}
	
	// Чтение вершин полилинии.
	while (reader.NextToken().IsVertex() == true) {
		CVertex vertex;
		
		// Одну команду необходимо пропустить, чтобы вершина начала считывать себя с нужного места.
		reader.GetToken();
		
		if (vertex.LoadFromDXF(reader) == true) m_vertexes.push_back(vertex);
		else {
			logMessage("ошибка: не удалось загрузить вершину полилинии");
			return false;
		}
	}
	
	// Чтение окончания последовательности.
	if (reader.NextToken().IsSequenceEnd() == true) {
		CSequenceEnd seqend;
		
		reader.GetToken();
		
		if (seqend.LoadFromDXF(reader) == false) {
			logMessage("ошибка: не может быть загружена команда завершения последовательности вершин полилинии");
			return false;
		}
	} else {
		logMessage("ошибка: нет команды завершения последовательности вершин полилинии");
		return false;
	}
	
	return true;
}

void CPolyline::RegisterPointsInGraph (CGraph & graph) {
	for (CVertexArray::iterator i = m_vertexes.begin(); i != m_vertexes.end(); ++i)
		(*i).RegisterPointInGraph(graph);
}

bool CPolyline::Compare (CCurve * curve) const {
	if (curve->IsPolyline()) {
		return (*this == *dynamic_cast<CPolyline *>(curve));
	} else return false;
}

bool CPolyline::operator == (const CPolyline & polyline) const {
	if (
		m_vertexes == polyline.m_vertexes &&
		fabs(m_startWidth - polyline.m_startWidth) < EPS &&
		fabs(m_endWidth - polyline.m_endWidth) < EPS
	) return true;
	
	return false;
}

// ------------------------------------------------------------------------------------------------

CCurve * AllocateCurveForToken (const CToken & token) {
	CCurve * curve = nullptr;

	if (token.IsLine()) curve = new CLine;
	else if (token.IsCircle()) curve = new CCircle;
	else if (token.IsArc()) curve = new CArc;
	else if (token.IsPolyline()) curve = new CPolyline;
	
	return curve;
}

} // namespace DXF