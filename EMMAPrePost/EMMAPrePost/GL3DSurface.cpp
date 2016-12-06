#include "StdAfx.h"
#include "GL3DSurface.h"

CGL3DSurface::CGL3DSurface(bool isPon, bool isSurf, bool isBL,  bool isJL)
{
	m_isPon = isPon;
	m_isSurf = isSurf;
	m_isBL = isBL;
	m_isJL = isJL;
	initData();
}

CGL3DSurface::~CGL3DSurface(void)
{
	delData();
}

void CGL3DSurface::outGL()
{
	// включаем ружимы
	glEnable(GL_POLYGON);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POINTS);			// чтобы точка была круглая	
	glEnable(GL_POINT_SMOOTH);
	glEnableClientState(GL_VERTEX_ARRAY);	// для работы с массивом
	glEnableClientState(GL_COLOR_ARRAY);
	glPointSize(7.0);

	if (m_isSurf) drawArray(GL_QUADS, m_surface, m_surfSz);
	if (m_isJL) drawArray(GL_LINES, m_joinLines, m_jlSz);
	if (m_isBL) drawArray(GL_LINES, m_basicLines, m_blSz);
	if (m_isPon) drawArray(GL_POINTS, m_points, m_ponSz);

	// возвращаем в обычное состояние OpenGL
	glDisable(GL_POLYGON);
	glDisable(GL_POLYGON_SMOOTH);
	glDisable(GL_POINTS);
	glDisable(GL_POINT_SMOOTH);
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
}

void CGL3DSurface::outGL(bool isPon, bool isSurf, bool isBL, bool isJL)
{
	m_isPon = isPon;
	m_isSurf = isSurf;
	m_isBL = isBL;
	m_isJL = isJL;
	outGL();
}

void CGL3DSurface::drawArray(GLenum mode, DBL *comb_ar, UNINT value_tops)const
{
	glColorPointer(4, GL_DOUBLE, 7 * sizeof(DBL), comb_ar);		// 3тий параметр - смещение между вершинами
	glVertexPointer(3, GL_DOUBLE, 7 * sizeof(DBL), &comb_ar[4]);
	glDrawArrays(mode, 0, value_tops);
}


void CGL3DSurface::initData()
{
	m_ponSz = 0;
	m_surfSz = 0;
	m_blSz = 0;
	m_jlSz = 0;
	m_surface = nullptr;
	m_basicLines = nullptr;
	m_joinLines = nullptr;
	m_points = nullptr;
}

void CGL3DSurface::delData()
{
	//if(m_surface) 
	delete[] m_surface;
	//if(m_basicLines) 
	delete[] m_basicLines;
	//if(m_joinLines) 
	delete[] m_joinLines;
	//if(m_points) 
	delete[] m_points;
	initData();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CGL3DSurface::setEvenData(const std::vector<std::vector<CPoint3D> > &ps)
{

	if (ps.size() <= 0 || ps[0].size() <= 0) {
		CDlgShowError cError(ID_ERROR_GL3DSURFACE_PS_WRONG); //_T("ps is wrong"));
		return;
	}

	delData();
	
	m_ponSz = ps.size() * ps[0].size();
	m_blSz = ps.size() * (ps[0].size() - 1) * 2;
	m_jlSz = (ps.size() - 1) * ps[0].size() * 2;
	m_surfSz = (ps.size() - 1) * (ps[0].size() - 1) * 6;//4;
	
	m_points = new DBL[m_ponSz * (4+3)];	// 4 - это RGBA; 3 - это xyz;
	m_basicLines = new DBL[m_blSz * (4+3)]; 
	m_joinLines = new DBL[m_jlSz * (4+3)];
	m_surface = new DBL[m_surfSz * (4+3)];

	UNLONG p = 0;
	UNLONG bl = 0;
	UNLONG jl = 0;
	UNLONG sf = 0;
	for(size_t i = 0; i < ps.size(); i++)
	{
		for(size_t i2 = 0; i2 < ps[i].size(); i2++)
		{
			if (ps[i].size() != ps[0].size()) {
				CDlgShowError cError(ID_ERROR_GL3DSURFACE_PS_WRONG); //_T("ps.GetSize() is wrong"));
				return;
			}

			// POINTS
			m_points[p++] = 1.0;	// сами точки красным
			m_points[p++] = 0.0;
			m_points[p++] = 0.0;
			m_points[p++] = 1.0;	// alpha 100%
			for(ptrdiff_t c = 0; c < 3; c++) m_points[p++] = ps[i][i2][c];	// cords

			// BASIC LINES
			if(i2 > 0) 
			{
				m_basicLines[bl++] = 0.0;	
				m_basicLines[bl++] = 0.0;
				m_basicLines[bl++] = 0.0;
				m_basicLines[bl++] = 1.0;	
				for(ptrdiff_t c = 0; c < 3; c++) m_basicLines[bl++] = ps[i][i2-1][c];

				m_basicLines[bl++] = 0.0;	
				m_basicLines[bl++] = 0.0;
				m_basicLines[bl++] = 0.0;
				m_basicLines[bl++] = 1.0;	
				for(ptrdiff_t c = 0; c < 3; c++) m_basicLines[bl++] = ps[i][i2][c];
			}
		
			// JOIN LINES
			if(i > 0) 
			{
				m_joinLines[jl++] = 0.0;	
				m_joinLines[jl++] = 0.0;
				m_joinLines[jl++] = 0.0;
				m_joinLines[jl++] = 1.0;	
				for(ptrdiff_t c = 0; c < 3; c++) m_joinLines[jl++] = ps[i-1][i2][c];

				m_joinLines[jl++] = 0.0;	
				m_joinLines[jl++] = 0.0;
				m_joinLines[jl++] = 0.0;
				m_joinLines[jl++] = 1.0;	
				for(ptrdiff_t c = 0; c < 3; c++) m_joinLines[jl++] = ps[i][i2][c];
			}

			// SURFACE
			if(i > 0 && i2 > 0)
			{

				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 1.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for (ptrdiff_t c = 0; c < 3; c++) {
					m_surface[sf++] = ps[i - 1][i2 - 1][c];
				}
				
				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 0.6;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.6;	
				for (ptrdiff_t c = 0; c < 3; c++) {
					m_surface[sf++] = ps[i - 1][i2][c];
				}
				
				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 1.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for (ptrdiff_t c = 0; c < 3; c++) {
					m_surface[sf++] = ps[i][i2][c];
				}				
				
				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 1.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for (ptrdiff_t c = 0; c < 3; c++) {
					m_surface[sf++] = ps[i][i2 - 1][c];
				}


				/*m_surface[sf++] = 0.0;	
				m_surface[sf++] = 1.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for(UNINT c = 0; c < 3; c++) m_surface[sf++] = ps[i-1][i2-1][c];
				
				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for(UNINT c = 0; c < 3; c++) m_surface[sf++] = ps[i-1][i2][c];
				
				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 1.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for(UNINT c = 0; c < 3; c++) m_surface[sf++] = ps[i][i2][c];
				
				
				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 1.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for(UNINT c = 0; c < 3; c++) m_surface[sf++] = ps[i][i2][c];

				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 1.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for(UNINT c = 0; c < 3; c++) m_surface[sf++] = ps[i][i2-1][c];
				
				m_surface[sf++] = 0.0;	
				m_surface[sf++] = 1.0;
				m_surface[sf++] = 0.0;
				m_surface[sf++] = 0.3;	
				for(UNINT c = 0; c < 3; c++) m_surface[sf++] = ps[i-1][i2-1][c];*/
			}

		}
	}

}
