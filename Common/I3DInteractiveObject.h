#pragma once

#include "InteractiveObject.h"
#include "3DPoint.h"
#include "3DBCAtom.h"
//#include "2DMesh.h" //TODO: нужно 3DMesh

class I3DInteractiveObject : public IInteractiveObject
{
public:
	I3DInteractiveObject(/*int id*/) //: IInteractiveObject(id)
	{}

	virtual Math::C3DRect GetBoundingBox() = 0;
	virtual bool GetBoundingBox(CRect3D &rect) = 0;

	/*!	Возвращает граничное условие для указанной точки, если точка на границе.	!*/
	virtual bool GetBC(/*const C3DMeshInterface *pMesh, */ ptrdiff_t nBoundaryNode, C3DBCAtom& bc ) = 0;
	/*!	Возвращает граничные условия для указанных точек, если точки на границе.	!*/
	virtual bool GetBC(/*const C3DMeshInterface *pMesh, */ std::vector<C3DBCAtom> *bc ) = 0;

	virtual void Preparations(const ITask *task) = 0;	// подготовка объектов для расчёта
	
	virtual bool Init() = 0;

	virtual ~I3DInteractiveObject() {}
};