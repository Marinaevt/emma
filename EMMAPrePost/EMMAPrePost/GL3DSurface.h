#pragma once

#include "../../Common/3DPoint.h"
#include "../../Common/IOFiles.h"
#include "resource.h"

//TODO: преобразовать в IO из-за MaterialDoc

/****************************************************************************************************/

// реализовано пока-что только для одинакового количества ключевых точек во входящем массиве!
class CGL3DSurface	// Surface - плоскость
{
public:
	CGL3DSurface(bool isPon = false, bool isSurf = true, bool isBL = true,  bool isJL = true);
	~CGL3DSurface(void);
	void setEvenData(const std::vector<std::vector<CPoint3D>> &ps);	// Even - равномерный



	void outGL(bool isPon, bool isSurf, bool isBL,  bool isJL);
	void outGL();
	virtual bool Load(CStorage& file){ return true; };
	virtual bool Save(CStorage& file){ return true; };

public:
			// .. is Draw
	bool m_isSurf;
	bool m_isJL;
	bool m_isBL;
	bool m_isPon;

protected:	
			// combaine arrays (col[4] + vert[3])
	DBL *m_surface;
	DBL *m_basicLines;
	DBL *m_joinLines;
	DBL *m_points;
			// НЕ реальные размеры, а количество вершин
	UNLONG m_ponSz;
	UNLONG m_surfSz;
	UNLONG m_blSz;
	UNLONG m_jlSz;

protected:
	void initData();
	void delData();
	inline void drawArray(GLenum mode, DBL *comb_ar, UNINT value_tops)const;
};
