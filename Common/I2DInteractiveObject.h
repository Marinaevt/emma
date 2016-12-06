#pragma once

#include "InteractiveObject.h"
#include "2DPoint.h"
#include "2DBCAtom.h"
#include "2DMesh.h"

class I2DInteractiveObject : public IInteractiveObject
{
public:
	I2DInteractiveObject(/*int id*/) //: IInteractiveObject(id)
	{}

	virtual Math::C2DRect GetBoundingBox() = 0;
	virtual bool GetBoundingBox(CRect2D &rect) = 0;

	/*!	Возвращает граничное условие для указанной точки, если точка на границе.	!*/
	virtual bool GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc) = 0;
	/*!	Возвращает граничные условия для указанных точек, если точки на границе.	!*/
	virtual bool GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom> *bc) = 0;

	virtual void Preparations(const ITask *task) = 0;	// подготовка объектов для расчёта
	
	virtual bool Init() = 0;

	virtual ~I2DInteractiveObject() {}
};