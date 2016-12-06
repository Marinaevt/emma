#include "StdAfx.h"
#include "2DPressure.h"


IOIMPL (C2DPressure, T2DPRESSURE)

double g_dPressure = 0.1;//Значение давления, которое читается из файла init (временно)

class InitPressure{
public:
	double m_dValue;	
	double m_dt0;
	InitPressure(){
		m_dValue = 0.1;
		m_dt0 = 0.01;

		CStorage file;

		if(file.Open(_T("init.txt."), CStorage::modeRead | CStorage::typeBinary)){ // Загрузка параметров из файла init  (временно)
			ReadValue(file, "<P>", &g_dPressure); //m_dValue здесь читать нельзя, поскольку при загрузке оно все равно изменится на то, которое было сохранено
			ReadValue(file, "<t0>", &m_dt0);
			file.Close();
		}else{
			CDlgShowError cError(ID_ERROR_2DPRESSURE_INIT_NOT_LOADED); //_T("2DPressure: init.txt not found"));
		}
	}
} g_IP;

C2DPressure::C2DPressure(){
	RegisterMember(&m_dValue);
	RegisterMember(&m_borderNodes);

	m_dValue = g_IP.m_dValue;	//MPa pressure
	m_dt0 = g_IP.m_dt0;

}

void C2DPressure::initBorderNodes(size_t nNodesCount, int bVal) {
	m_borderNodes.SetSize(nNodesCount);
	for (size_t i = 0; i < m_borderNodes.GetSize(); i++) {
		m_borderNodes[i] = bVal;
	}
}

void C2DPressure::setBorderNode(size_t nNode, int bVal) {

	if (nNode >= m_borderNodes.GetSize()) {
		CDlgShowError cError(ID_ERROR_2DPRESSURE_NNODE_WRONG); //_T("2DPressure: nNode is wrong"));
		return;
	}
	m_borderNodes[nNode] = bVal;
}

bool C2DPressure::GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc)
{
	m_dValue = g_dPressure;
	size_t nBoundaryNodesCount = pMesh->m_bordernodes.GetSize();

	//проверка на номер узла
	if (nBoundaryNode >= nBoundaryNodesCount) {
		CDlgShowError cError(ID_ERROR_2DPRESSURE_NNODE_WRONG);
		return false;
	}
	
	size_t nBoundaryNode2 = nBoundaryNode > 0 ? nBoundaryNode - 1 : nBoundaryNodesCount - 1;	//соседний узел

	if(m_borderNodes[nBoundaryNode] || m_borderNodes[nBoundaryNode2]){
		Math::C2DPoint normS(0,0);
		if(m_borderNodes[nBoundaryNode2]){
			normS += pMesh->GetNormalS(nBoundaryNode2, 1);
		}
		if(m_borderNodes[nBoundaryNode]){
			normS += pMesh->GetNormalS(nBoundaryNode, -1);
		}
		double elapsedTyme = (GetParent() ? GetParent()->GetElapsedTime() : 0.)/m_dt0;
		double dVal = -m_dValue*min(elapsedTyme, 1);

		bc.setLoad(dVal*normS.x, dVal*normS.y);

		//TRACE(CString(_T("[")) + AllToString(nBoundaryNode) + CString(_T("]:  Px = ")) + AllToString(bc.getQx()) + 
		//	  CString(_T("  Py = ")+ AllToString(bc.getQy()))+ CString(_T("\n\r")));
		return true;
	}
	return false;
}

bool C2DPressure::GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom>  *bc)
{
	size_t nBoundaryNodesCount = pMesh->m_bordernodes.GetSize();
	m_dValue = g_dPressure;
	bool bResultBC = false;	//проверяем, есть ли хоть в 1м узле ГУ

	for(size_t nBoundaryNode=0; nBoundaryNode < nBoundaryNodesCount; nBoundaryNode++){

		size_t nBoundaryNode2 = nBoundaryNode > 0 ? nBoundaryNode - 1 : nBoundaryNodesCount - 1;	//соседний узел

		if(m_borderNodes[nBoundaryNode] || m_borderNodes[nBoundaryNode2]){
			Math::C2DPoint normS(0,0);
			if(m_borderNodes[nBoundaryNode2]){
				normS += pMesh->GetNormalS(nBoundaryNode2, 1);
			}
			if(m_borderNodes[nBoundaryNode]){
				normS += pMesh->GetNormalS(nBoundaryNode, -1);
			}
			double elapsedTyme = (GetParent() ? GetParent()->GetElapsedTime() : 0.)/m_dt0;
			double dVal = -m_dValue*min(elapsedTyme, 1);

			bc->at(nBoundaryNode).setLoad(dVal*normS.x, dVal*normS.y);

			//TRACE(CString(_T("[")) + AllToString(nBoundaryNode) + CString(_T("]:  Px = ")) + AllToString(bc->at(nBoundaryNode).getQx()) + 
			//	  CString(_T("  Py = ")+ AllToString(bc->at(nBoundaryNode).getQy()))+ CString(_T("\n\r")));
			bResultBC = true;
		}
	}
	return bResultBC;
}

C2DPressure::~C2DPressure(void)
{
}
