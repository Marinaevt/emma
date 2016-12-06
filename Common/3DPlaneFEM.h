#pragma once

#include "i3dinteractiveobject.h"
#include "AllTypes.h"
#include "2DMesh.h"
#include "SLAE.h"
#include "3DBCAtom.h"

class C3DPlaneFEM :	public I3DInteractiveObject
{
	GETTYPE (T3DPLANEFEM)
protected:
	std::vector<int> m_BBN, m_TBN; // bottom and top border nodes


	void CalcFEM(DBL dt); ///< основная процедура, производящая вычисления МКЭ
	void CalcIntDef(DBL dt);  ///< накопление интенсивности деформации в элементах
	void CalcAvgDef(DBL dt); ///< расчет средней деформации
	void InitStep();
	void MakeSmooth(); ///< осреднение (сглаживание) полей
public:

	IO::CInterfaceAdapter m_mesh_adapt;				  ///< сетка
	IO::CMaterial m_mat;							  ///< свойства материала
	CSLAE m_slae;									  ///< СЛАУ (система линейных алгебраических уравнений)
	IO::CSimpleArray<C3DBCAtom, C3DBCATOMARRAY> m_bc; ///< ГУ (граничные условия)

public:
	C3DPlaneFEM();

	void SetMesh(C2DMeshInterface* newmesh)
	{
		m_mesh_adapt.SetObject(newmesh);
	}
	bool SetMaterial(const CString& path)
	{
		return m_mat.ReLoadMaterial(path);
	}

	//virtual bool Load(CStorage& file);
	//virtual bool Save(CStorage& file);

	virtual bool Init();
	virtual void Calc(DBL dt);
	virtual void Move(DBL dt);
	
	virtual void DrawGL(GLParam &parameter);
	virtual void DrawGL3D(GLParam &parameter);

	virtual CRect3D GetBoundingBox();
	virtual bool GetBoundingBox(CRect3D &rect);

	virtual bool IsInside(const Math::C3DPoint& point);

	virtual bool GetBC(const Math::C3DPoint& node, C3DBCAtom& bc);
	/*!	Возвращает граничное условие для указанной точки, если точка на границе.	!*/
	virtual bool GetBC(/*const C3DMeshInterface *pMesh, */ ptrdiff_t nBoundaryNode, C3DBCAtom& bc);
	/*!	Возвращает граничные условия для указанных точек, если точки на границе.	!*/
	virtual bool GetBC(/*const C3DMeshInterface *pMesh, */ std::vector<C3DBCAtom> *bc);
	
	// подготовка объектов для расчёта
	virtual void Preparations(const ITask *task);

	virtual ~C3DPlaneFEM();
};

IODEF (C3DPlaneFEM)
