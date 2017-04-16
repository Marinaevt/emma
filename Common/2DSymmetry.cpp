#include "StdAfx.h"
#include "2DSymmetry.h"

IOIMPL (C2DSymmetry, T2DSYMMETRY)

C2DSymmetry::C2DSymmetry(int nMode){
	m_mode() = nMode;
	m_left() = false;
	RegisterMember(&m_axle);
	RegisterMember(&m_mode);
	RegisterMember(&m_left);
}

bool C2DSymmetry::GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc)
{
	if (pMesh == nullptr) {
		CDlgShowError cError(ID_ERROR_2DSYMMETRY_NO_FE_MESH); //_T("pMesh is null"));
		return false;
	}

	Math::C2DPoint node = pMesh->GetBorderNode(nBoundaryNode);
	switch(m_mode()){
	case modeY:
		if(node.x <= m_axle().GetX(0) + EPS && !m_left || node.x + EPS >= m_axle().GetX(0) && m_left()){
			bc.setSymY();
			return true;
		}
		break;
	case modeX:
		if(node.y <= m_axle().GetY(0) + EPS && !m_left || node.y + EPS >= m_axle().GetY(0) && m_left()){
			bc.setSymX();
			return true;
		}
		break;
	case modeAny:
		if (m_axle().IsLeft(node) && m_left())
		{

			//bc.setSymAny();
			CDlgShowError cError(ID_ERROR_2DSYMMETRY_NO_ANGLE);		// НАДО УСТАНОВИТЬ УГОЛ!
			return false;
		}
		break;
	}
	return false;
}

bool C2DSymmetry::GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom> *bc)
{
	if (pMesh == nullptr) {
		CDlgShowError cError(ID_ERROR_2DSYMMETRY_NO_FE_MESH);
		return false;
	}
	
	for(size_t nBoundaryNode=0; nBoundaryNode < pMesh->m_bordernodes.GetSize(); nBoundaryNode++){
		Math::C2DPoint node = pMesh->GetBorderNode(nBoundaryNode);
		switch(m_mode()){
		case modeY:
			if((node.x <= m_axle().GetX(0) + EPS && !m_left) || (node.x + EPS >= m_axle().GetX(0) && m_left())){
				bc->at(nBoundaryNode).setSymY();
			}
			break;
		case modeX:
			if((node.y <= m_axle().GetY(0) + EPS && !m_left) || (node.y + EPS >= m_axle().GetY(0) && m_left())){
				bc->at(nBoundaryNode).setSymX();
			}
			break;
		case modeAny:
			if (m_axle().IsLeft(node) && m_left()){
				CDlgShowError cError(ID_ERROR_2DSYMMETRY_NO_ANGLE); // НАДО УСТАНОВИТЬ УГОЛ!
				return false;
			}
			break;
		}
	}
	return true;
}

void C2DSymmetry::Preparations(const ITask *task){

}

C2DSymmetry::~C2DSymmetry(void)
{
}
