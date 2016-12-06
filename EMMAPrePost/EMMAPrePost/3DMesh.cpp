#include "StdAfx.h"
#include "3DMesh.h"


IOIMPL (C3DMesh, C3DMESH)

C3DMesh::C3DMesh(void) {
	m_LayersAmount() = 3;
	RegisterMember(&m_MeshPoints);
	RegisterMember(&m_MeshElements);
	RegisterMember(&m_Z);
	RegisterMember(&m_LayersAmount);
	//RegisterMember(&m_BoundingBox);	//TODO: перевести в IO
};

void C3DMesh::ClearAll()
{
	m_Z().clear();
	m_MeshElements().clear();
	m_MeshPoints().clear();
}

int C3DMesh::GetLayersAmount()
{
	return m_LayersAmount - 1;
}

void C3DMesh::PrepareElements(int nEl)
{
	m_MeshElements().resize((m_LayersAmount - 1)*nEl*3);
}

void C3DMesh::SetCardinalityHeightArray()
{
	m_Z().resize(m_LayersAmount);
}

void C3DMesh::SetCardinalityHeightArray(int m)
{
	m_LayersAmount() = m;
	m_Z().resize(m);
}

void C3DMesh::SetHeights() //по-умолчанию
{
	m_Z()[0] = 0;
	m_Z()[1] = 4;
	//m_Z()[2] = 8;
	//m_Z()[3] = 12;
}

void C3DMesh::SetHeights(double m)
{
	m_Z()[0] = 0;
	for(int i = 1; i < m_LayersAmount - 1; i++)
	{
		m_Z()[i] = m_Z()[i - 1] + m;
	}
}

void C3DMesh::SetHeight(int n, double h)
{
	m_Z()[n] = h;
}

double C3DMesh::GetHeight(int n)
{
	return m_Z()[n]-m_Z()[n+1];
}

void C3DMesh::SetElements(int i, int v1, int v2, int v3, int v4)
{
	m_MeshElements()[i].n0 = v1;
	m_MeshElements()[i].n1 = v2;
	m_MeshElements()[i].n2 = v3;
	m_MeshElements()[i].n3 = v4;
}

C3DElement4& C3DMesh::GetElement(int n)
{
	return m_MeshElements()[n];
}

void C3DMesh::SetLayersAmount(int v)
{
	m_LayersAmount() = v;
}

CRect3D C3DMesh::GetDraftRect()
{
	return m_BoundingBox;
}

void C3DMesh::SetDraftRect(DBL mx,DBL Mx,DBL my,DBL My)
{
	m_BoundingBox = CRect3D(mx,Mx,my,My,0,(m_LayersAmount - 2)*4);
}

