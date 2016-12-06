#include "StdAfx.h"
#include "MaterialDoc.h"


IOIMPL(CMaterialDoc, CMATERIALDOC)

CMaterialDoc::CMaterialDoc(void) : CEMMADoc(RUNTIME_CLASS(CMaterialView), IDR_EMMAPrePostTYPE)
{
	RegisterMember(&m_Material);
	//RegisterMember(&m_surfaces);	//TODO:
	SetIconName(IDS_ICON_MATERIALDOC);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CMaterialDoc::~CMaterialDoc(void)
{
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CMaterialDoc::LoadMaterial(void){
	CFileDialog fileDlg(1, _T("spf"), 0, 4|2, _T("material files (*.spf)|*.spf|"));	//1 - open, 0 - save
	if(fileDlg.DoModal() == IDOK)
	{
		CString sPathName = fileDlg.GetPathName();
		bool bLoad = m_Material.ReLoadMaterial(sPathName);
		SetName(fileDlg.GetFileName());
		//MatToPoints();	// СОЗДАЕМ КЛЮЧЕВЫЕ ТОЧКИ 
		return bLoad;
	}
	return false;
}

void CMaterialDoc::FillPropList(CMFCPropertyGridCtrl *pGrid){
}

void CMaterialDoc::UpdateProp(double *pEval, UNLONG *pId){
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CMaterialDoc::MatToPoints()
{

	/*
	// ВРЕМЕННО //  
	DBL minT = 800;		// диапазон 800 -> 1400
	DBL maxT = 1400;
	DBL deltaT = 100;

	DBL minE = 0;		// диапазон 0 -> 1
	DBL maxE = 3;
	DBL deltaE = 0.15;
				 
	DBL minV = 0;		// диапазон 0 -> 1
	DBL maxV = 100;
	DBL deltaV = 0.5;
	// END //

	for(DBL T = minT; T < maxT; T += deltaT)
	{
		std::vector<std::vector<CPoint3D>> kp;
		for(DBL e = minE; e < maxE; e += deltaE) 	
		{
			kp.SetSize(kp.GetSize() + 1);
			for(DBL v = minV; v < maxV; v += deltaV) 
			{
				kp[kp.GetSize() - 1].Add(CPoint3D(e, v/100, m_Material.Material().Si(e, v, T)/100));
			}
		}
		m_surfaces().SetSize(m_surfaces().GetSize() + 1);
		m_surfaces()[m_surfaces().GetSize() - 1].setEvenData(kp);
	}

	//*/
}

void CMaterialDoc::PointsToMat()
{
	/*
	for(INT i = 0; i < m_surfaces().GetSize(); i++)
	{
		std::vector<std::vector<CPoint3D>> kp;
		m_surfaces()[i].getEvenData(kp);

		//....
	}
	//*/
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMaterialDoc::DrawGL()
{ 
	/*
	if (m_surfaces().GetSize() <= 0) {
		//CDlgShowError cError(_T("Размер m_surfaces <= 0"));
		return;
	}

	//for(INT_PTR i = 0; i < m_surfaces().GetSize(); i++)
	{
		m_surfaces()[i].outGL(false, true, true, true);
	}

	//..
	//*/
}
/*for(int T = 0; T < kp.sz(); T++)
	{
		for(int e = 0; e < kp[T].sz(); e++) 		
		{
			for(int v = 0; v < kp[T][e].sz(); v++) 
			{
				if(surface && e > 0 && v > 0)
				{
					glColor4d(0, 1, (1 - T/kp.sz()), 0.3);
					glBegin(GL_POLYGON);
					glVertex3d(kp[T][e-1][v-1]);
					glVertex3d(kp[T][e-1][ v ]);
					glVertex3d(kp[T][ e ][ v ]);
					glVertex3d(kp[T][ e ][v-1]);
					glEnd();
				}
				if(joinLines && e > 0)
				{
					glColor4d(0, 1, (1 - T/kp.sz()), 1);
					glLineWidth((GLfloat)1.7);
					glBegin(GL_LINES);
					glVertex3d(kp[T][e-1][v]);		
					glVertex3d(kp[T][ e ][v]);
					glEnd();
				}
				if(basicLines && v > 0)
				{
					glColor4d(0, 1, (1 - T/kp.sz()), 1);
					glLineWidth((GLfloat)1.7);
					glBegin(GL_LINES);
					glVertex3d(kp[T][e][ v ]);		
					glVertex3d(kp[T][e][v-1]);
					glEnd();
				}
				if(allPoints)
				{
					glColor4d(1, 0.3, 0, 1);
					glBegin(GL_POINTS);	
					glVertex3d(kp[T][e][v]);
					glEnd();
				}

			}		
		}
	}*/

/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CMaterialDoc::createKeyPoint()
{
	std::vector<std::vector<C3DPoint>> &kp = m_keyPoints;
	//kp.~std::vector();
	kp.SetSize(0);

	// ВРЕМЕННО //
	DBL T = 1100;		

	DBL minE = -1;	// деапозон -1 -> 0
	DBL maxE = 0;	//
	DBL m_deltaE;
				 
	DBL minV = -1;	// деапозон -1 -> 1
	DBL maxV = 1;	//
	DBL deltaV;
	// END //

	for(DBL e = minE; e < maxE; e += deltaE) 	
	{
		kp.SetSize(kp.GetSize() + 1);
		for(DBL v = minV; v < maxV; v += deltaE) 
		{
			kp[kp.GetSize() - 1].Add(C3DPoint(e, v, m_Material.Si(e, v, T)));
		}
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ВНЕМАНИЕ!!! толко для случая, когда частота точек для каждого 'e' одинаковая !!! 
void CMaterialDoc::DrawGL()
{ 
	std::vector<std::vector<C3DPoint>> &kp = m_keyPoints;
	//ASSERT(kp[0].GetSize() == kp[kp.GetSize() - 1].GetSize());

	bool basicLines = true;	   // is draw	(мужду разными v)
	bool joinLines = true;	  // is draw	(между разными e)
	bool surface = true;	 // is draw		(сама поверхность)
	bool allPoints = true;  // is draw		(сами кл. точки)

	glEnable(GL_POLYGON);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POINTS);			// чтобы точка была круглая	
	glEnable(GL_POINT_SMOOTH);	
	glPointSize(7.0);

	for(int e = 0; e < kp.GetSize(); e++) 		
	{
		for(int v = 0; v < kp[e].GetSize(); v++) 
		{
			if(surface && e > 0 && v > 0)
			{
				glColor4d(0, 1, 0.2, 0.3);
				glBegin(GL_POLYGON);
				glVertex3d(kp[e-1][v-1]);
				glVertex3d(kp[e-1][ v ]);
				glVertex3d(kp[ e ][ v ]);
				glVertex3d(kp[ e ][v-1]);
				glEnd();
			}
			if(joinLines && e > 0)
			{
				glColor4d(0, 1, 0.2, 1);
				glLineWidth((GLfloat)1.7);
				glBegin(GL_LINES);
				glVertex3d(kp[e-1][v]);		
				glVertex3d(kp[ e ][v]);
				glEnd();
			}
			if(basicLines && v > 0)
			{
				glColor4d(0, 1, 0.2, 1);
				glLineWidth((GLfloat)1.7);
				glBegin(GL_LINES);
				glVertex3d(kp[e][ v ]);		
				glVertex3d(kp[e][v-1]);
				glEnd();
			}
			if(allPoints)
			{
				glColor4d(1, 0.3, 0, 1);
				glBegin(GL_POINTS);	
				glVertex3d(kp[e][v]);
				glEnd();
			}

		}		
	}

	glColor4d(1, 0, 0, 1);
	glBegin(GL_POINTS);	
	if(m_beActKeyPoint) glVertex3d(m_actKeyPoint);
	glEnd();
	//..
}
*/
