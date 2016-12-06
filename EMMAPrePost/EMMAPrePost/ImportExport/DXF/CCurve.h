#pragma once

#include <vector>
#include "../../../../Common/2DOutline.h"

namespace DXF {

class CCurve;
typedef std::vector<CCurve *> CCurvePtrArray;

}

#include "CPoint.h"
#include "CReader.h"
#include "utility.h"

namespace DXF {

class CGraph;

//! Абстрактная кривая.
/*!
	Базовый класс для всех типов кривых.
	\sa CLine, CCircle, CArc, CPolyline
*/
class CCurve {
public:
	CCurve () {}
	virtual ~CCurve () {}
	
	virtual bool LoadFromDXF (CReader & reader) = 0;
	
	//! Регистрация точек, содержащихся в кривой, в чертеже.
	/*!
		Проверяет каждую точку, которая есть в кривой, на наличие в чертеже "graph". Если такая
		точка в чертеже уже есть, то указатель на неё освобождается и переназначается на указатель
		на аналогичную точку из чертежа. Если же такой точки в чертеже нет, то она туда
		добавляется.
	*/
	virtual void RegisterPointsInGraph (CGraph & graph) = 0;
	
	//! Сравнение кривых.
	/*!
		Сравнивает кривые, независимо от их типа.
		Для этого сначала сравнивает типы, и если типы одинаковы, вызывает оператор равенства.
	*/
	virtual bool Compare (CCurve * curve) const = 0;
	
	virtual bool IsLine () const { return false; }
	virtual bool IsCircle () const { return false; }
	virtual bool IsArc () const { return false; }
	virtual bool IsPolyline () const { return false; }

	virtual int ExportToOutline (C2DOutline * outline) const {
		logMessage("Функция \"ExportToOutline\" у этого элемента пока не поддерживается.");
		return 0;
	} //раскомментировать
};

// ------------------------------------------------------------------------------------------------

//! Отрезок прямой.
/*!
	Примитив, содержащий указатели на точки своего начала и конца.
	\sa CCurve, CCircle, CArc, CPolyline
*/
class CLine : public CCurve {
public:
	CLine () {
		m_startPoint = nullptr;
		m_endPoint = nullptr;
	}
	~CLine () {
		delete m_startPoint;
		delete m_endPoint;
	}
	
	bool LoadFromDXF (CReader & reader);
	
	void RegisterPointsInGraph (CGraph & graph);
	
	const CPoint & StartPoint () const { return *m_startPoint; }
	const CPoint & EndPoint () const { return *m_endPoint; }
	
	bool Compare (CCurve * curve) const;
	bool operator == (const CLine & line) const;
	
	bool IsLine () const { return true; }
	
	int ExportToOutline (C2DOutline * outline) const {
		C2DCurve * line = dynamic_cast<C2DCurve *>(IO::CreateObject(C2DLINESEG));
		
		line->SetNodes(
			m_startPoint->ExportToOutline(outline),
			m_endPoint->ExportToOutline(outline)
		);
		
		return outline->AddCurve(line);
	} //раскомментировать
	
private:
	// Указатели на точки отрезка. Освобождать память не нужно, т.к. на самом деле точки лежат в
	// массиве точек в классе "CContour".
	CPoint * m_startPoint;
	CPoint * m_endPoint;
};

// ------------------------------------------------------------------------------------------------

//! Окружность.
/*!
	Примитив, содержащий указатель на точку своего центра, а так же свой радиус.
	\sa CCurve, CLine, CArc, CPolyline
*/
class CCircle : public CCurve {
public:
	CCircle () {
		m_centerPoint = nullptr;
		m_radius = 0.0;
	}
	~CCircle () {
		delete m_centerPoint;
	}
	
	bool LoadFromDXF (CReader & reader);
	
	void RegisterPointsInGraph (CGraph & graph);
	
	const CPoint & CenterPoint () const { return *m_centerPoint; }
	double Radius () const { return m_radius; }
	
	bool Compare (CCurve * curve) const;
	bool operator == (const CCircle & circle) const;
	
	bool IsCircle () const { return true; }

	int ExportToOutline (C2DOutline * outline) const {


		//пока пусто, потому что я не знаю, как работать с окружностями в C2DOutline.
		logMessage("Функция \"ExportToOutline\" у этого элемента пока не поддерживается.");
		return 0;

	} //раскомментировать
	
private:
	CPoint * m_centerPoint;
	
	double m_radius;
};

// ------------------------------------------------------------------------------------------------

//! Дуга окружности.
/*!
	Примитив, содержащий указатель на точку своего центра, а так же свои радиус, начальный и конечный углы.
	\sa CCurve, CLine, CCircle, CPolyline
*/
class CArc : public CCurve {
public:
	CArc () {
		m_centerPoint = nullptr;
		m_radius = 0.0;
		m_startAngle = 0.0;
		m_endAngle = 0.0;
	}
	~CArc () {
		delete m_centerPoint;
	}
	
	bool LoadFromDXF (CReader & reader);
	
	void RegisterPointsInGraph (CGraph & graph);
	
	const CPoint & CenterPoint () const { return *m_centerPoint; }
	double Radius () const { return m_radius; }
	double StartAngle () const { return m_startAngle; }
	double EndAngle () const { return m_endAngle; }
	
	bool Compare (CCurve * curve) const;
	bool operator == (const CArc & arc) const;
	
	bool IsArc () const { return true; }
	
	int ExportToOutline (C2DOutline * outline) const {
		C2DCurve * arc = dynamic_cast<C2DCurve *>(IO::CreateObject(C2DCIRCARC));
		CPoint startPoint, endPoint; // CSimplePoint startPoint, endPoint;

		startPoint.set(
			m_centerPoint->x + m_radius * cos(M_PI * m_startAngle / 180),
			m_centerPoint->y + m_radius * sin(M_PI * m_startAngle / 180),
			0
		);

		endPoint.set(
			m_centerPoint->x + m_radius * cos(M_PI * m_endAngle / 180),
			m_centerPoint->y + m_radius * sin(M_PI * m_endAngle / 180),
			0
		);

		arc->SetNodes(
			startPoint.ExportToOutline(outline),
			endPoint.ExportToOutline(outline)
		);
		
		dynamic_cast<C2DCircleArc *>(arc)->SetRadius(m_radius);
		
		return outline->AddCurve(arc);
	} //раскомментировать
	
private:
	CPoint * m_centerPoint;
	
	double m_radius;
	double m_startAngle;
	double m_endAngle;
};

// ------------------------------------------------------------------------------------------------

//! Вершина полилинии.
/*!
	Представляет собой точку с параметрами, такими, как начальная и конечная ширина линии, кривизна линии в этой точке и пр.
	\sa CPolyline
*/
class CVertex {
public:
	CVertex () {
		m_point = nullptr;
		m_startWidth = m_endWidth = 0.0;
		m_curvature = 0.0;
		m_tangentDirection = 0.0;
		m_flag = 0.0;
	}
	~CVertex () {
		delete m_point;
	}
	
	bool LoadFromDXF (CReader & reader);
	
	void RegisterPointInGraph (CGraph & graph);
	
	const CPoint & Point () const { return *m_point; }
	double StartWidth () const { return m_startWidth; }
	double EndWidth () const { return m_endWidth; }
	double Curvature () const { return m_curvature; }
	double TangentDirection () const { return m_tangentDirection; }
	double Flag () const { return m_flag; }
	
	bool operator == (const CVertex & vertex) const;

private:
	CPoint * m_point;
	
	double m_startWidth;
	double m_endWidth;
	double m_curvature;
	double m_tangentDirection;
	double m_flag;
};

typedef std::vector<CVertex> CVertexArray;

//! Конец последовательности примитивов типа ATTRIB или VERTEX.
/*!
	Бесполезный класс, необходимый только для того, чтобы пропустить пару команд в файле.
*/
class CSequenceEnd {
public:
	CSequenceEnd () {}
	~CSequenceEnd () {}

	bool LoadFromDXF (CReader & reader);
};

//! Полилиния.
/*!
	Представляет собой набор вершин.
	\sa CCurve, CLine, CCircle, CArc
*/
class CPolyline : public CCurve {
public:
	CPolyline () {
		m_startWidth = m_endWidth = 0.0;
	}
	~CPolyline () {
	}
	
	bool LoadFromDXF (CReader & reader);
	
	void RegisterPointsInGraph (CGraph & graph);
	
	CVertexArray::iterator VertexesBegin () { return m_vertexes.begin(); }
	CVertexArray::iterator VertexesEnd () { return m_vertexes.end(); }
	
	CVertexArray::const_iterator VertexesBegin () const { return m_vertexes.begin(); }
	CVertexArray::const_iterator VertexesEnd () const { return m_vertexes.end(); }
	
	double StartWidth () const { return m_startWidth; }
	double EndWidth () const { return m_endWidth; }
	
	bool Compare (CCurve * curve) const;
	bool operator == (const CPolyline & polyline) const;
	
	bool IsPolyline () const { return true; }
	
private:
	CVertexArray m_vertexes;
	
	double m_startWidth;
	double m_endWidth;
	
	// Ещё должны быть разные флаги, но они пока не поддерживаются.
};

//! Создание кривой в зависимости от команды.
/*!
	По команде "token" определяет, какой тип кривой нужно создать, а затем выделяет память под эту
	кривую и возвращает указатель на неё. Если команда не опознана, то возвращается нулевой
	указатель ("nullptr").
*/
CCurve * AllocateCurveForToken (const CToken & token);

} // namespace DXF

#include "CGraph.h"