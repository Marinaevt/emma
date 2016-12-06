#pragma once
#include "IOInterface.h"
#include "2DPoint.h"
#include "GLParam.h"
#include "./Logger/Logger.h"
#include "math.h"

// **************************** C2DNode ****************************

namespace IO {
typedef CTemplate<Math::C2DPoint, CIO2DPOINT> C2DPoint;
}

class C2DNode : public IO::CClass
{
	GETTYPE (C2DNODE)
protected:
	IO::C2DPoint m_point;
public:
	C2DNode(const Math::C2DPoint& p = Math::C2DPoint::Zero)
	{
		RegisterMember(&m_point);
		SetPoint(p);
	}

	void SetPoint(const Math::C2DPoint& p)
	{
		m_point() = p;
	}

	Math::C2DPoint GetPoint()
	{
		return m_point();
	}
	
	IO::C2DPoint GetIOPoint(){
		return m_point;
	} 

	virtual Math::C2DRect GetBoundingBox();
	virtual bool GetBoundingBox(Math::C2DRect &rect);

	virtual void DrawGL(GLParam &parameter);

	//virtual bool Copy(IO::CInterface* ptr);

	virtual ~C2DNode()
	{}
};

IODEF (C2DNode)

// **************************** C2DCurve ****************************
// C2DCURVE == C2DLINESEG

class C2DOutline;

class C2DCurve : public IO::CClass
{
	GETTYPE (C2DCURVE)
protected:
	C2DOutline* m_parent;
	IO::CSize_t m_start, m_end;
public:

	C2DCurve(int start = 0, int end = 0);

	void SetParent(C2DOutline* parent)
	{
		m_parent = parent;
	}

	size_t GetStart() const
	{
		return m_start();
	}

	size_t GetEnd() const
	{
		return m_end();
	}

	C2DNode* GetStartNode() const;
	C2DNode* GetEndNode() const;

	bool SetNodes(size_t start, size_t end);
	
	virtual void SwapNodes()
	{
		Math::swap(m_start(), m_end());
	}

	// проверка на "равенство", т.е. совпадение основных параметров
	virtual bool IsEqual(const C2DCurve* c);

	// есть ли общий узел?
	bool IsConnectedTo(const C2DCurve* c);

	//Сколькими узлами соединены кривые
	int HowConnectedTo(const C2DCurve* c);

	bool HasNode(int node) const
	{
		return (m_start() == node) || (m_end() == node);
	}

	// длина отрезков
	virtual DBL CalcLength();

	virtual bool GetPoint(DBL l, Math::C2DPoint& p, Math::C2DPoint& tau);
	//Находим точку на середине длины
	virtual Math::C2DPoint GetPointHalfLenght();

	// записывает точку пересечения отрезка p1p2 с кривой и возвращает:
	// -1 - ошибка
	//  0 - более одного пересечения
	//  1 - одно пересечение (возвращает точку пересечения в result)
	//  2 - нет пересечений
	virtual int GetIntersection(const Math::C2DPoint& p1, const Math::C2DPoint& p2, Math::C2DPoint& result);

	// записывает все точки пересечения прямой l с кривой в result
	virtual bool GetAllIntersections(const Math::C2DLine& l, std::vector<Math::C2DPoint>& result);
	
	// записывает в result точку на кривой, ближайшую к указанной и возвращает:
	// -1 - ошибка
	//  0 - ближайшая точка является начальным узлом кривой
	//  1 - ближайшая точка является конечным узлом кривой
	//  2 - ближайшая точка лежит где-то на кривой
	virtual int GetClosestPoint(const Math::C2DPoint& p, Math::C2DPoint& result);

	virtual Math::C2DRect GetBoundingBox();
	
	//Yura
	virtual C2DCurve* Clone();

	virtual void DrawGL(GLParam &parameter);

	//bool Copy(IO::CInterface* ptr);

	virtual ~C2DCurve()
	{}
};

IODEF(C2DCurve)

// **************************** C2DCircleArc ****************************

class C2DCircleArc : public C2DCurve
{
	GETTYPE (C2DCIRCARC)
protected:
	IO::CDouble m_rad; ///< радиус дуги (отрицательный радиус означает, что центр дуги находится с другой стороны)
	IO::CBool m_cw; ///< clockwise, параметр, позволяющий выбрать одну из двух дуг

	bool CalcCircleData(double& phi0, double& phi1, Math::C2DPoint& center) const;
public:

	C2DCircleArc();
	C2DCircleArc(int start, int end, DBL rad);

	void SetRadius(DBL rad)	{ m_rad() = rad; }
	DBL  GetRadius() const 	{ return m_rad(); }
	void SetCW(bool cw) { m_cw() = cw; }
	bool GetCW() const	{ return m_cw(); }

	//////ROMA
	bool CheckRadius(int acc);

	// проверка на "равенство", т.е. совпадение основных параметров
	virtual bool IsEqual(const C2DCurve* c);	

	// длина дуги
	virtual DBL CalcLength();

	virtual bool GetPoint(DBL l, Math::C2DPoint& p, Math::C2DPoint& tau);
	//Получение точки на середине длины дуги
	virtual Math::C2DPoint C2DCircleArc::GetPointHalfLenght();

	virtual int GetIntersection(const Math::C2DPoint& p1, const Math::C2DPoint& p2, Math::C2DPoint& result);
	virtual int GetIntersection2Points(const Math::C2DPoint& p1, const Math::C2DPoint& p2, Math::C2DPoint& result1, Math::C2DPoint& result2);
	
	// записывает все точки пересечения прямой l с кривой в result
	virtual bool GetAllIntersections(const Math::C2DLine& l, std::vector<Math::C2DPoint>& result);

	//Yura записывает точка(и) пересечения с дугами 
	virtual int GetIntersectionArcs(C2DCircleArc& c1, Math::C2DPoint& result);
	virtual int GetIntersectionArcs2Points(C2DCircleArc& c1, Math::C2DPoint& result, Math::C2DPoint& result2);
	virtual bool IfPoint(const Math::C2DPoint& p1, DBL eps, Math::C2DPoint& result);
	/////

	virtual void SwapNodes();

	virtual int GetClosestPoint(const Math::C2DPoint& p, Math::C2DPoint& result);

	virtual Math::C2DRect GetBoundingBox();
	
	virtual void DrawGL(GLParam &parameter);

	//Yura
	virtual C2DCircleArc* Clone();

	//bool Copy(IO::CInterface* ptr);
	
	virtual ~C2DCircleArc()
	{}
};

IODEF(C2DCircleArc)

// **************************** C2DContour ****************************

// здесь curve_i - индекс кривой в массиве кривых в C2DOutline
//       index   - индекс в массиве m_indexes в C2DContour

class C2DContour : public IO::CClass
{
	GETTYPE (C2DCONTOUR)
protected:
	C2DOutline* m_parent;
	IO::CSize_tArray m_indexes;
	IO::CBool m_closed;

	//--- contour approximation points cache ---
	IO::CBool m_bcache;
	std::vector<Math::C2DPoint> m_cache;
public:
	
	C2DContour(C2DOutline* parent = nullptr);

	void SetParent(C2DOutline* parent);

	//int operator [] (size_t i);

	size_t& At(size_t i)
	{
		return m_indexes()[i];
	}


	// найти локальный индекс кривой curve_i в контурном массиве m_indexes
	int GetCurveIndexByIndex(size_t curve_i);
	// найти локальный индекс кривой по её указателю
	int GetCurveIndexByCurve(C2DCurve* pCur);
	// получить указатель на кривую в C2DOutline по ее индексу в массиве m_indexes
	C2DCurve* GetCurve(size_t index);
	// очистить контур
	void Clear();

	bool AddCurve(int curve_i);

	// возвращает true, если:
	// 1. контур пустой
	// 2. контур состоит из одной кривой (nHead = nTail = 0)
	// 3. контур не имеет раздвоений
	// 4. контур замкнут (nHead = nTail = 0)
	bool CheckContour(size_t &nHead, size_t &nTail);

	bool OrderContour();

	bool DelCurve(size_t curve_i);

	// количество кривых в контуре
	/*size_t GetCount() const
	{
		return m_indexes().size();
	} //*/

	// количество кривых в контуре (название привычнее)
	size_t GetCurveCount() const
	{
		return m_indexes().size();
	}

	//контур замкнут или открыт?
	bool IsClosed() const
	{
		return m_closed();
	}

	//контур открыт
	void UnClose()
	{
		m_closed() = false;
	}

	//замкнуть контур
	void Close(){
		m_closed() = true;
	}

	virtual void DrawGL(GLParam &parameter);

	bool IsCache() const
	{
		return m_bcache();
	}

	void SetCache(bool set){
		m_bcache() = set;
	}

	void FillCache();

	const std::vector<Math::C2DPoint>& GetCache()
	{
		return m_cache;
	}

	//Находим среднюю точку для всего контура
	Math::C2DPoint GetCenterPoint();

	//Находит сумму длин всех кривых в контуре
	DBL CalcLength();

	virtual ~C2DContour()
	{}
};

IODEF (C2DContour)

// **************************** C2DOutline ****************************

class C2DOutline : public IO::CClass
{
	GETTYPE (C2DOUTLINE)
protected:
	IO::CPtrArray<CIO2DNODEPTRARRAY>  m_nodes;
	IO::CPtrArray<CIO2DCURVEPTRARRAY> m_curves;
	IO::CPtrArray<CIO2DCONTPTRARRAY>  m_contours;

public:
	C2DOutline();

	void SetParentToAll(); // установить себя в качестве parent для всех кривых и контуров

	virtual int  CreateNode(TYPEID type, const Math::C2DPoint& p = Math::C2DPoint::Zero);
	virtual int  AddNode(C2DNode* pNode);
	virtual bool DelNode(size_t index);
	virtual C2DNode* GetNode(size_t index);
	size_t GetNodeCount() const;
	void ClearAll();
	////Roma
	//// Объединяет узлы, стоящии рядом. acc - точность.
	void DelRepNodes(double acc);
	void DelRepNodes2(double acc);
	size_t FindNumNode(C2DNode *pNode);  ///ищет номер узла в outline 
	/////
	virtual int  CreateCurve(TYPEID type, int start, int end);
	virtual ptrdiff_t  AddCurve(C2DCurve* pCurve);
	virtual bool DelCurve(size_t index);
	//возвращает кривую по индексу Outline
	virtual C2DCurve* GetCurve(size_t index);
	// возвращает индекс кривой в Outline
	virtual int GetCurveIndexByCurve(C2DCurve* pCurve);
	size_t GetCurveCount() const;
	
	// Эта функция должна удалить "Висящие" узлы и линии, 
	// найти внешний контур и перенумеровать его по часовой стрелке
	bool PrerareToMeshing();

	virtual int  AddNodeToCurve(int curve_i, DBL length);

	virtual int  CreateContour(TYPEID type);
	virtual size_t  AddContour(C2DContour* pContour);
	virtual bool DelContour(size_t index);
	virtual C2DContour* GetContour(size_t index);
	size_t GetContourCount() const;
	bool FindContours(double ribLen);  // поиск контуров среди всех кривых с учетом их пересечения

	virtual bool AddCurveToContour(int curve_i, int contour_i);
	virtual bool DelCurveFromContour(size_t curve_i, size_t contour_i);	// удаляем кривую curve_i из контура contour_i
	virtual bool DelCurveFromContours(size_t curve_i);		// удаляем кривую curve_i из всех контуров

	bool GetPoint(int c, DBL l, Math::C2DPoint& p, Math::C2DPoint& tau); // выдает точку на контуре c, отстоящую от начального узла на длину l вдоль всей кривой, а также касательную в этой точке

	Math::C2DPoint GetCenterPoint();	// Получаем центральную точку чертежа (центр масс)

	// записывает точку пересечения отрезка p1p2 с кривой curve_num и возвращает:
	// -1 - ошибка
	//  0 - более одного пересечения
	//  1 - одно пересечение (возвращает точку пересечения в result)
	//  2 - нет пересечений
	virtual int GetIntersection(const Math::C2DPoint& p1, const Math::C2DPoint& p2, Math::C2DPoint& result, int& curve_num);

	// рассчитывает все пересечения прямой с контуром c
	virtual bool GetAllIntersections(size_t c, const Math::C2DLine& line, std::vector<Math::C2DPoint>& points);

	// записывает в result точку на кривой curve_num, ближайшую к указанной и возвращает:
	// -1 - ошибка
	//  0 - ближайшая точка является начальным узлом кривой
	//  1 - ближайшая точка является конечным узлом кривой
	//  2 - ближайшая точка лежит где-то на кривой
	virtual int GetClosestPoint(const Math::C2DPoint& p, Math::C2DPoint& result, int& curve_num);
	
	//Возвращает номер ближайшего узла или -1, если узлов нет
	//в len записывается расстояние до ближайшего узла, если таковой имеется
	virtual int GetClosestNode(const Math::C2DPoint& p, DBL& len);

	//Возвращает наименьшую длину кривой по всем контурам
	DBL GetSmallestCurveLenght();
	//Возвращает номер наименьшей по длине кривой по всем контурам или -1
	int GetSmallestCurveIndex();

	// проверяет, находится ли точка p внутри контура c (если он замкнут)
	// возвращает -2 - если ошибка, -1 - если точка снаружи, 0 - если на контуре, 1 - если внутри
	int IsInsideContour(const Math::C2DPoint& p, int c);

	// проверяет, находится ли точка p внутри хотя бы одного из контуров (= любого)
	int IsInside(const Math::C2DPoint& p);

	virtual Math::C2DRect GetBoundingBox();
	virtual bool GetBoundingBox(CRect2D &rect);
	
	virtual void DrawGL(GLParam &parameter);

	DBL GetSumLengthOfAllCurves();

	//сохраняем чертёж в файл с точностью до 12 знака
	CString SaveMetaToString(size_t precisizon = 12);
	virtual bool Load(CStorage& file);

	virtual bool Copy(IO::CInterface* ptr);

	virtual void WriteToLog();

	virtual ~C2DOutline();
};

IODEF(C2DOutline)

// **************************** C2DFacet_common ****************************

class C2DFacet_common : public C2DOutline
{
public:

	int m_nNode1, m_nNode2, m_nNode3;
	Math::C2DPoint m_pVect1, m_pVect2, m_ugl;
	DBL m_dRad;

	C2DFacet_common(){
		m_nNode1 = m_nNode2 = m_nNode3 = 0;
		m_dRad = 0.0;
	}

	bool Facet_ugl(Math::C2DPoint *vect1, Math::C2DPoint *vect2, Math::C2DPoint w1b, Math::C2DPoint w1e, Math::C2DPoint w2b, Math::C2DPoint w2e);
	Math::C2DPoint Facet_scale_rotate(Math::C2DPoint pbegin, Math::C2DPoint p1, Math::C2DPoint p2, Math::C2DPoint p3);
	int Facet_create_arc(Math::C2DPoint p1, Math::C2DPoint p2, Math::C2DPoint p3, int n1, int n2, DBL rad);
	void SetTempFacetPointsOrder(Math::C2DPoint s1, Math::C2DPoint e1, Math::C2DPoint s2, Math::C2DPoint e2);	
	//void CorrectArcs(short t1, short t2, Math::C2DPoint s1, Math::C2DPoint e1, Math::C2DPoint s2, Math::C2DPoint e2);
};

IODEF(C2DFacet_common)

// **************************** C2DFacet_2lines ****************************

class C2DFacet_2lines : public C2DFacet_common
{
public:	
	C2DFacet_2lines(Math::C2DPoint p, DBL r){
		m_ugl = p;
		m_dRad = r;
	}

	void CreateData();
	bool Create();
	bool CreateFacet(Math::C2DPoint tper1, Math::C2DPoint tper2);
	
	~C2DFacet_2lines(){
	};
};

IODEF(C2DFacet_2lines)

// **************************** C2DFacet_line_arc ****************************

class C2DFacet_line_arc : public C2DFacet_common
{
public:
	short m_Cur, m_Arc; //порядковые номера кривых
	Math::C2DPoint  m_P1up, m_P2up, 
					m_P1down, m_P2down;

	C2DFacet_line_arc(){
		m_Cur = m_Arc = 0;
	}

	C2DFacet_line_arc(Math::C2DPoint p, DBL r){
		m_ugl = p;
		m_dRad = r;
		m_Cur = m_Arc = 0;
	}

	void CreateData(short t1, short t2, bool cw, DBL r);	//Заполнение Facet_line_arc
	void CorrectArcs(short t1, short t2, int n1, int n2, C2DOutline* out);	//Корректировка дуг
	bool Create();
	bool Facet_line_arc_point(int up, int down, Math::C2DPoint *result);	//Поиск точки пересечения
	bool Facet_line_arc_2points(int up, int down, Math::C2DPoint *result);	//Поиск 2х точек пересечения
	bool CreateFacet(Math::C2DPoint tper1, Math::C2DPoint tper2, Math::C2DPoint o2);
	
	~C2DFacet_line_arc(){
	};
};

IODEF(C2DFacet_line_arc)

// **************************** C2DFacet_2arcs ****************************

class C2DFacet_2arcs : public C2DFacet_common
{
public:
	int m_nArc1_up, m_nArc1_down, 
		m_nArc2_up, m_nArc2_down;
	
	C2DFacet_2arcs(){
		m_nArc1_up = m_nArc1_down = 0; 
		m_nArc2_up = m_nArc2_down = 0;
	}

	C2DFacet_2arcs(Math::C2DPoint p, DBL r){
		m_ugl = p;
		m_dRad = r;
		m_nArc1_up = m_nArc1_down = 0; 
		m_nArc2_up = m_nArc2_down = 0;
	}

	void CreateData(bool cw1, DBL r1, bool cw2, DBL r2);	//Заполнение Facet_2arcs
	void CorrectArcs(int n1, int n2, C2DOutline* out);	//Корректировка дуг
	bool Create();
	bool Facet_2arcs_point(Math::C2DPoint l1,Math::C2DPoint l2, Math::C2DPoint *result);	//Поиск точки пересечения
	bool Facet_2arcs_2points(Math::C2DPoint l1,Math::C2DPoint l2, Math::C2DPoint* result);	//Поиск 2х точек пересечения
	bool CreateFacet(Math::C2DPoint tper1, Math::C2DPoint tper2, Math::C2DPoint l1,Math::C2DPoint l2, DBL r1, DBL r2);	

	~C2DFacet_2arcs(){
	}
};

IODEF(C2DFacet_2arcs)