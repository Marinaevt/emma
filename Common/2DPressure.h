#pragma once
#include "AllTypes.h"
#include "I2DInteractiveObject.h"
#include "2DBCAtom.h"
#include "ITask.h"
#include "Material.h"
#include "IOFiles.h"
#include "DlgShowError.h"

class C2DPressure : public I2DInteractiveObject
{
	IO::CDouble m_dValue;
	IO::CIntArray m_borderNodes;
	DBL m_dt0;//продолжительность набора давления


	GETTYPE (T2DPRESSURE)
public:

	C2DPressure();

	void SetValue(double dPressure){m_dValue = dPressure;}

	void initBorderNodes(size_t nNodesCount, int bVal = 0);

	void setBorderNode(size_t nNode, int bVal = 1);

	virtual bool Init() { return true; }
	virtual void Calc(double dt) {}
	virtual void Move(double dt) {}

	virtual void DrawGL(GLParam &parameter){}
	virtual void DrawGL3D(GLParam &parameter) {}

	//Задание ГУ для одного граничного узла
	virtual bool GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc);
	//Задание ГУ для вектора граничных узлов
	virtual bool GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom> *bc);
	// подготовка объектов для расчёта
	virtual void Preparations(const ITask *task){}

	virtual Math::C2DRect GetBoundingBox(){ return Math::C2DRect(); }
	virtual bool GetBoundingBox(CRect2D &rect){ return false; }

	virtual ~C2DPressure(void);
};

IODEF (C2DPressure)
