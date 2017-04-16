#pragma once
#include "AllTypes.h"
#include "SLAE.h"
#include "2DBCAtom.h"
#include "2DMesh.h"
#include "2DRigid.h"

//#include"2DOutline.h"

#include "2DPoint.h"
#include "I2DInteractiveobject.h"
#include "ITask.h"
#include "./Logger/AdvancedLogger.h"
#include "DlgShowError.h"
#include "IOFiles.h"
#include "math.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//! Класс сборки сетки, ГУ, материалов
class C2DPlaneFEM :	public I2DInteractiveObject
{
	GETTYPE (T2DPLANEFEM)
protected:
	std::vector<int> BBN, TBN; // bottom and top border nodes


	void CalcFEM(DBL dt); ///< основная процедура, производящая вычисления МКЭ
	void CalcIntDef(DBL dt);  ///< накопление интенсивности деформации в элементах
	void CalcAvgDef(DBL dt); ///< расчет средней деформации (avg_d)
	void CalcAvgDefSpeed(); // расчёт средней скорости деформации (avg_ds)

	void InitStep();
	void MakeSmoothDef(); ///< осреднение (сглаживание) полей деформации (smo_d)
	void MakeSmoothDefSpeed(); // сглаживание полей скоростей деформации (smo_ds)
	void FieldsToNodes(DBL dt);	///< раскидываем значения полей по узлам
	
	void CalcForce(DBL dt);	///< считаем силу в элементах (основание заготовки)
	//void CalcPressure(DBL dt); // считаем давление в элементах (основание заготовки)

public:

	IO::CInterfaceAdapter m_mesh_adapt; // сетка
	IO::CMaterial m_mat;				// свойства материала
	CSLAE m_slae;						// СЛАУ (система линейных алгебраических уравнений)
	IO::CSimpleArray<C2DBCAtom, C2DBCATOMARRAY> m_bc;				// ГУ (граничные условия)
	IO::CDouble m_CurTime;				// текущее время расчёта

	DBL m_FirstVolSum;		//суммарная НАЧАЛЬНАЯ площадь заготовки
	DBL m_CurrentVolSum;		//суммарная площадь заготовки в момент t
	bool m_FirstVolSumStart;	//флаг, что сумма расчитана

public:
	C2DPlaneFEM();

	void InitMesh(C2DMeshInterface* newmesh)
	{
		m_mesh_adapt.SetObject(newmesh);
	}

	bool SetMaterial(const CString& path)
	{
		return m_mat.ReLoadMaterial(path);
	}

	//virtual bool Load(CStorage& file);
	//virtual bool Save(CStorage& file);
	C2DMeshInterface *GetMesh(){return dynamic_cast<C2DMeshInterface *>(m_mesh_adapt());}

	//Инициализация КЭ сетки (сохранение в 2DBuilder)
	bool InitMesh(const CMesh *pInitMesh);
	//Инициализация материала (сохранение в 2DBuilder)
	bool InitMaterial(const CString& path);
	void InitBC();

	virtual bool Init();			//! Инициализация в Srv
	virtual void Calc(DBL dt);
	virtual void Move(DBL dt);

	virtual void Preparations(const ITask *task);	// подготовка объектов для расчёта
	
	virtual void DrawGL(GLParam &parameter);
	virtual void DrawGL3D(GLParam &parameter);

	virtual Math::C2DRect GetBoundingBox();
	virtual bool GetBoundingBox(CRect2D &rect);

	virtual bool IsInside(const Math::C2DPoint& point)
	{ return false;	}

	virtual bool GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc)
	{ return false; }
	virtual bool GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom> *bc)
	{ return false; }

	void WriteToLogNewLines(int numNewLines);	//! делаем кол-во пустых строк
	virtual void WriteBCToLog();		//! Записываем ГУ в лог
	void WriteEFieldsToLog();	//! Записываем поля элементов в лог
	void WriteEFieldsToLog(DBL dt);	//! Записываем поля элементов в лог

	virtual ~C2DPlaneFEM();

};

IODEF (C2DPlaneFEM)
