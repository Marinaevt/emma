#pragma once
#include "afxwin.h"
#include <gl/GL.h>
#include <gl/GLU.h>
#include "../../../Common/2DOutline.h"
#include "../GL2DView.h"
///////////////////////////////////////////////////
/*!  
	CGLModelPreview - класс превью загружаемого чертежа 
*/
class CGLModelPreview : public CGL2DView
{
	DECLARE_DYNCREATE(CGLModelPreview);
public:
	//конструктор 
	CGLModelPreview(void);
	//установка модели 
	bool SetModel(C2DOutline *Outline);

	//очистка окна превью и сохранённой модели
	void ClearModel();
	//перегруженый метод вывода сцены 
	virtual void DrawScene(void);
	//определение области отрисовки 
	afx_msg LRESULT OnDraftRect(WPARAM wParam, LPARAM lParam);
	//перед созданием окна 
	virtual int PreCreateWindow(CREATESTRUCT& cs);
	//во время движения колеса мыши (отменяет зум)
	afx_msg int OnMouseWheel(UNINT nFlags, short zDelta, CPoint pt){return 1;};
	//деструктор ~CGLModelPreview(void);
	~CGLModelPreview(void);
protected:
	DECLARE_MESSAGE_MAP()
public:
	//чертеж
	C2DOutline *m_pModelPreviewOutline;
};

