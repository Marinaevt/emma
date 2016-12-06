#include "StdAfx.h"
#include "3DMesher.h"


IOIMPL (C3DMesher, C3DMESHER)


C3DMesher::C3DMesher(void) : CEMMADoc(RUNTIME_CLASS(CChild3DView), IDR_EMMAPrePostTYPE)
{
	m_flag() = 0;
}

C3DMesher::~C3DMesher(void)
{
}

void C3DMesher::DrawGL(GLParam &parameter)
{
//	dynamic_cast<C2DMesher*>(m_pParentDoc)->DrawGL();
	//вот это очень плохо v
	if(m_flag() == 0)
	{
	//	Construct3DMesh();
		GetTetra();
		m_flag() = 1;
	}
//тут надо закомментировать
	//glBegin(GL_POINTS);
	//glColor3d(0,1,0);
	//	for(int i = 0; i < m_3DMesh.m_MeshPoints().GetCount(); i++)
	//	{
	//		glVertex3dv(m_3DMesh.m_MeshPoints()[i]);
	//	}
	//glEnd();
	glBegin(GL_TRIANGLES);
		for(size_t i = 0; i < m_3DMesh.m_MeshElements().size(); i++)
		{
			for(int j = 0; j < 4; j++)
			{
				//glColor3d(0, 0, 1);
				glColor4d(i%3, (i + 1)%3, (i+2)%3, 0.2);
				glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[j]]);
				glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[(j + 1)%4]]);
				glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[(j + 2)%4]]);
			}
		}
	glEnd();
//тут конец комментов
	//glBegin(GL_LINES);
	//for(int i = 0; i < m_3DMesh.m_MeshElements().GetCount(); i++)
	//	{
	//		glEnable(GL_BLEND);
	//		glLineWidth(1);
	//		//glColor4d(i%3, (i + 1)%3, (i+2)%3, 0.2);
	//		glColor4d(1, 0, 0, 0.2);
	//		for(int j = 0; j < 3; j++)
	//		{
	//			glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[j]]);
	//			glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[j + 1]]);
	//		}
	//		//glColor4d(1, 0, 0, 0.5f);
	//			glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[0]]);
	//			glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[3]]);
	//			glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[0]]);
	//			glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[2]]);
	//		//glColor4d(0,0,1, 0.5f);
	//		glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[3]]);
	//		glVertex3dv(m_3DMesh.m_MeshPoints()[m_3DMesh.m_MeshElements()[i].n[1]]);

	//		glDisable(GL_BLEND);

	//	}
	//glEnd();
}

void C3DMesher::FillPropList(CMFCPropertyGridCtrl *pGrid)
{
	if(pGrid == 0) return;
	
	pGrid->RemoveAll();	
	pGrid->AddProperty(AddProp(_T("Layers"), m_3DMesh.GetLayersAmount() - 1,21));

	CMFCPropertyGridProperty *t_pProp = new CMFCPropertyGridProperty(_T("Length"));
	for(int i = 1; i < m_3DMesh.GetLayersAmount(); i++)
	{
		t_pProp->AddSubItem(AddProp(_T("Length"), m_3DMesh.GetHeight(i),(22 + i)));
	}

	pGrid->AddProperty(t_pProp);
	pGrid->Invalidate();
}

CMFCPropertyGridProperty * C3DMesher::AddProp(CStringW name, double eval, size_t pos)
{
	CMFCPropertyGridProperty *pProp;
	
	pProp = new CMFCPropertyGridProperty(name, eval, name);
	pProp->SetData(pos);

	if (pos == 5) pProp->AllowEdit(0);
	
	return pProp;
}

void C3DMesher::GetBoundingBox()
 {
	DBL mx = 0.0, Mx = 0.0, my = 0.0, My = 0.0;
	CMesh* pMesh = dynamic_cast<C2DMesher*>(m_pParentDoc)->GetMesh();
	for(size_t i = 0; i < pMesh->m_nPointsNum; i++)
	{
		if(mx > pMesh->m_Points[i].x)
			mx = pMesh->m_Points[i].x;
		if(Mx < pMesh->m_Points[i].x)
			Mx = pMesh->m_Points[i].x;

		if(my > pMesh->m_Points[i].y)
			my = pMesh->m_Points[i].y;
		if(My < pMesh->m_Points[i].y)
			My = pMesh->m_Points[i].y;
	}
	m_3DMesh.SetDraftRect(mx,Mx,my,My);
 }

void C3DMesher::Construct3DMesh()
{
	//LOGGER.Init(CString("..\\..\\Logs\\emma_3dmesher_constructedmesh.txt"));
	m_3DMesh.SetCardinalityHeightArray();
	m_3DMesh.SetHeights(4);
	CMesh* pMesh = dynamic_cast<C2DMesher*>(m_pParentDoc)->GetMesh();
	int t_p = (pMesh->m_nPointsNum)*(m_3DMesh.GetLayersAmount() + 1);
	m_3DMesh.m_MeshPoints().resize(t_p);

	for(int j = 0; j < m_3DMesh.GetLayersAmount(); j++)
	{
		for(size_t i = 0; i < pMesh->m_nPointsNum; i++)
		{
			CPoint3D& p = m_3DMesh.m_MeshPoints()[i + (j*pMesh->m_nPointsNum)];
			p = CPoint3D(pMesh->m_Points[i], m_3DMesh.GetHeight(j));
			
			//DLOG(AllToString(p.x) + CString(_T(",")) + AllToString(p.y) + CString(_T(",")) + AllToString(p.z), log_info);
		}
	}
}

bool C3DMesher::GetTetra()
{
	CMesh *pMesh = dynamic_cast<C2DMesher*>(m_pParentDoc)->GetMesh();
	int nElNum = pMesh->m_nElNum;
	m_3DMesh.PrepareElements(nElNum);

	int nPointsNumMax = pMesh->m_nPointsNum;
	Math::C2DPoint* nPointsNumDec = pMesh->m_Points;

	for(int e = 0; e < nElNum; e++)
	{
		int n1 = -1, n2 = -1;

		for(int i = 0; i < 3; i++)
		{ 
			if(!pMesh->m_Elements[e].m_status[i]) 
				return false;
			if(pMesh->m_Elements[e].m_mark[i])
			{
				if(n1 == -1)
					n1 = i;
				else
					n1 = -2;
			}
			else
			{
				if(n2 == -1)
					n2 = i;
				else
					n2 = -2;
			}
		}
		if(n1 * n2 > 0) 
			return false;
		if(n1 >= 0)
		{
			for(int z = 0; z < m_3DMesh.GetLayersAmount() - 1; z++)
			{
				m_3DMesh.SetElements(3*e + (z*nElNum),
					(int)(Get2DEPoint(e, n1) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,(n1+1)%3) - nPointsNumDec) + z*nPointsNumMax,
					(int)(Get2DEPoint(e,(n1+2)%3) - nPointsNumDec) + z*nPointsNumMax,
					(int)(Get2DEPoint(e,n1) - nPointsNumDec) + z*nPointsNumMax);

				m_3DMesh.SetElements(3*e + 1 + (z*nElNum),
					(int)(Get2DEPoint(e,n1) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,(n1+1)%3) - nPointsNumDec) + z*nPointsNumMax,
					(int)(Get2DEPoint(e,(n1+1)%3) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,(n1+2)%3) - nPointsNumDec) + (z + 1)*nPointsNumMax);

				m_3DMesh.SetElements(3*e + 2 + (z*nElNum),
					(int)(Get2DEPoint(e,n1) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,(n1+1)%3) - nPointsNumDec) + z*nPointsNumMax,
					(int)(Get2DEPoint(e,(n1+2)%3) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,(n1+2)%3) - nPointsNumDec) + z*nPointsNumMax);
			}
		}
		if(n2 >= 0)
		{
			for(int z = 0; z < m_3DMesh.GetLayersAmount() - 1; z++)
			{
				m_3DMesh.SetElements(3*e + (z*nElNum),
					(int)(Get2DEPoint(e,n2) - nPointsNumDec) + z*nPointsNumMax,
					(int)(Get2DEPoint(e,(n2 + 1)%3) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,(n2 + 2)%3) - nPointsNumDec) + z*nPointsNumMax,
					(int)(Get2DEPoint(e,(n2 + 1)%3) - nPointsNumDec) + z*nPointsNumMax);  
			
				m_3DMesh.SetElements(3*e + 1 + (z*nElNum),
					(int)(Get2DEPoint(e,n2) - nPointsNumDec) + z*nPointsNumMax,
					(int)(Get2DEPoint(e,(n2 + 1)%3) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,n2) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,(n2+2)%3) - nPointsNumDec) + (z + 1)*nPointsNumMax);

				m_3DMesh.SetElements(3*e + 2 + (z*nElNum),
					(int)(Get2DEPoint(e,(n2 + 2)%3) - nPointsNumDec) + (z + 1)*nPointsNumMax, 
					(int)(Get2DEPoint(e,(n2 + 1)%3) - nPointsNumDec) + (z + 1)*nPointsNumMax,
					(int)(Get2DEPoint(e,n2) - nPointsNumDec) + z*nPointsNumMax,
					(int)(Get2DEPoint(e,(n2 + 2)%3)- nPointsNumDec) + z*nPointsNumMax);
			}
		}
	}
	return true;
}
 
Math::C2DPoint* C3DMesher::Get2DEPoint(int e, int n)
{
	if (e < 0 || n < 0) return nullptr;

	Math::C2DPoint* pPoint = dynamic_cast<C2DMesher*>(m_pParentDoc)->GetMesh()->m_Elements[e].m_Points[n];
	return pPoint;
}

void C3DMesher::UpdateProp(double *pEval, UNLONG *pId)
{
	if( (pEval == 0) || (pId == 0)) return;
	/*switch(*pId)
	{
	case 21:*/
	int nAmount = static_cast<int>(*pEval);
		m_3DMesh.SetLayersAmount(nAmount);
	/*	break;
	case 22:
		m_3DMesh.SetHeight(1,*pEval);
		break;
	default:
		m_3DMesh.SetHeight(1,*pEval);
		break;B
	}*/
	m_3DMesh.ClearAll();
	Construct3DMesh();
	GetTetra();
}