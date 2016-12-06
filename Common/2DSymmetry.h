#pragma once

#include "I2DInteractiveObject.h"
#include "AllTypes.h"
#include "../Common/2DMotion.h"
#include "2DBCAtom.h"

class C2DSymmetry : public I2DInteractiveObject
{
	IO::C2DLine m_axle; ///< ось симметрии
	IO::CBool m_left; ///< левая или правая полуплоскость

	IO::CInt m_mode;


	GETTYPE (T2DSYMMETRY)
public:
	enum{
		modeX = 1,  //Относительно X
		modeY = 2,  //Относительно Y
		modeAny = 0 //Относительно заданной оси
	};

	C2DSymmetry(int nMode = modeX);

	virtual bool Init() { return true; }
	virtual void Calc(double dt) {}
	virtual void Move(double dt) {}

	virtual void DrawGL(GLParam &parameter){ /*shape.DrawGL();*/ }
	virtual void DrawGL3D(GLParam &parameter) {}
	virtual bool IsInside(const Math::C2DPoint& point)
	{ return false; }

	virtual bool GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc);
	virtual bool GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom> *bc);
	virtual void Preparations(const ITask *task);	// подготовка объектов для расчёта

	virtual Math::C2DRect GetBoundingBox(){ return Math::C2DRect(); }
	virtual bool GetBoundingBox(CRect2D &rect){ return false; };

	virtual ~C2DSymmetry(void);
};

IODEF (C2DSymmetry)
