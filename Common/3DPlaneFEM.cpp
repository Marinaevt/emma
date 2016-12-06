#include "stdafx.h"
#include "3DPlaneFEM.h"

IOIMPL(C3DPlaneFEM, T3DPLANEFEM)

C3DPlaneFEM::C3DPlaneFEM()
{
	RegisterMember(&m_mesh_adapt);
	//RegisterMember(&m_mat);
	RegisterMember(&m_bc);
}

bool C3DPlaneFEM::Init() {
	return false;
};
void C3DPlaneFEM::Calc(DBL dt) {};
void C3DPlaneFEM::Move(DBL dt) {};

void C3DPlaneFEM::DrawGL(GLParam &parameter) {};
void C3DPlaneFEM::DrawGL3D(GLParam &parameter) {};

CRect3D C3DPlaneFEM::GetBoundingBox() {
	return CRect3D(0, 0, 0, 0, 0, 0);
};
bool C3DPlaneFEM::GetBoundingBox(CRect3D &rect) {
	return false;
};

bool C3DPlaneFEM::IsInside(const Math::C3DPoint& point)
{
	return false;
}

bool C3DPlaneFEM::GetBC(const Math::C3DPoint& node, C3DBCAtom& bc)
{
	return false;
}
/*!	Возвращает граничное условие для указанной точки, если точка на границе.	!*/
bool C3DPlaneFEM::GetBC(/*const C3DMeshInterface *pMesh, */ ptrdiff_t nBoundaryNode, C3DBCAtom& bc) {
	return false;
};
/*!	Возвращает граничные условия для указанных точек, если точки на границе.	!*/
bool C3DPlaneFEM::GetBC(/*const C3DMeshInterface *pMesh, */ std::vector<C3DBCAtom> *bc) {
	return false;
};

// подготовка объектов для расчёта
void C3DPlaneFEM::Preparations(const ITask *task) {
	
};

C3DPlaneFEM::~C3DPlaneFEM() {
};