#pragma once
#include "Globals.h"
#include "RibLine.h"
#include "../2DOutline.h"
#include "../IOInterface.h"
#include "../AllTypes.h"
#include "../ResourceCommon.h"

/////////////////////////////////////////////////////////////////////////////
//! Сетка

class CMesh : public IO::CClass
{
	GETTYPE (С2DMESH)
public:

	//TODO: IO
	CFinEl	*m_Elements; //массив элементов
	Math::C2DPoint	*m_Points; //массив узлов
	FRONTINF	*m_FrInf;

	size_t m_nElNum; //количество элементов
	size_t m_nElNumMax; //максимально возможное количество элементов
	size_t m_nPointsNum; //количество узлов
	size_t m_nPointsNumMax; //максимально возможное количество узлов
	
	IO::CUnInt m_nMeshType;		//тип КЭ сетки

	IO::C2DPoint m_MeshPoint;	//Точка сгущения КЭ сетки
	IO::CDouble m_dMeshRad;		//Радиус сгущения
	IO::CBool m_IsSetPoint;		//Установлена точка сгущения или нет

	IO::CDouble m_dRibLen;	//длина ребра сетки
	CRibLine m_RibLine;

public:
	CMesh();
	~CMesh();

	size_t GetNodesNum() const {return m_nPointsNum;}
	size_t GetElNum() const {return m_nElNum;}
	const Math::C2DPoint &GetNode(size_t nNode) const;

	const CFinEl &GetElement(size_t nEl) const;

	ptrdiff_t GetNodeInEl(size_t nEl, int nIndex)const{
		return m_Elements[nEl].m_Points[nIndex] - m_Points;
	}
	int GetStar(int nPoint, Math::C2DPoint ***Star);
	void Optimize(int n=3);
	void Renumber(Math::C2DPoint ***Stars, int *StarSize);
	void SetOutline();
	int	AddPoint(const Math::C2DPoint &DP, const FRONTINF &FI);
	int AddElement(const CFinEl &FinEl);

	//! Создание сетки (тип, длина разбиения, угол, чертёж)
	void SetMeshType(UNINT id);
	UNINT GetMeshType() { return m_nMeshType(); };
	bool GenerateMesh(double dRibLen, C2DOutline* Outline);
	bool GenerateMeshFront(double dRibLen, double Ang, C2DOutline* Outline);
	bool GenerateMeshMSW(double dRibLen, C2DOutline* Outline);
	void DeleteMesh();
	void DeleteMeshPoint();

	CFinEl *NewElement();
	void ConstructMap();
	int  FindNeighbor(size_t num, int rib, int *NeigRib = nullptr);
	bool SetMarks(int index);
	void DrawGl(GLParam &parameter) const;
	void WriteToLog() const;
	double GetRibLen(const Math::C2DPoint &Point);

};
IODEF (CMesh)
