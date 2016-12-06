#pragma once
#include "2DScale.h"
#include <gl/GL.h>
#include "GLFont.h"

// Класс, преставляющий параметры для Viewport. Позволяет получать данные в виде
// CRect и GLint[]. Так же позволяет менять размеры и/или сдвигать координаты. 
class GLViewport{

protected:
	GLint m_X;
	GLint m_Y;
	GLsizei m_Width;
	GLsizei m_Height;

public:
	GLViewport(){
		Set(0, 0, 0, 0);
	}

	GLViewport(GLint x, GLint y, GLsizei width, GLsizei height){
		Set(x, y, width, height);
	}

	GLViewport(const GLint *viewport){
		Set(viewport);
	}

	GLViewport(const CRect &rect){ 
		Set(rect);
	}

	void Set(GLint x, GLint y, GLsizei width, GLsizei height);

	void Set(const GLint *viewport);

	void Set(const CRect &rect);

	void ShiftBy(GLint x, GLint y);

	void ShiftBy(CPoint point);

	void ResizeLeftBottom(GLint x, GLint y);

	void ResizeLeftBottom(CPoint point){
		ResizeLeftBottom(point.x, point.y);
	}

	void ResizeRightTop(GLint width, GLint height);

	void ResizeRightTop(CPoint point){
		ResizeRightTop(point.x, point.y);
	}

	CRect ToCRect() const{
		return CRect(m_X, m_Y - m_Height, m_X + m_Width, m_Y);
	}

	void FillArray(GLint out[4]) const;

};

//------------------- class ScaleVisualizer --------------

class CScaleVisualizer{

protected:
	
	GLuint m_DisplayListIndex;
	GLint m_ClientGLintArray[4] = {0, 0, 0, 0};
	
	GLdouble m_Trasparency;
	GLdouble m_ScaleTop;
	GLdouble m_ScaleBottom;
	GLdouble m_ScaleLeft;
	GLdouble m_ScaleRight;
	GLdouble m_FontHeight;
	
	CGLFont *m_Font;

public:
	bool m_IsMoved;
	GLViewport m_ScaleViewport;

public:
	CScaleVisualizer();

	~CScaleVisualizer()
	{
		ClearDisplayList();
	}

	void DrawScale(const C2DScale &scale);

	bool IsScaleVisible(){
		return fabs(m_Trasparency) > EPS;
	}

	void ClearDisplayList();

	void SetFont(CGLFont *font){
		m_Font = font;
	}

	void SetTrasparency(int percent);

protected:

	void SetProjection();
	void RestoreProjection();

	void DrawViewportBorder();
	void DrawScaleBorder();

	std::vector<size_t> GetDrawableMarks(const C2DScale &scale);

	void DrawScaleMarks(const C2DScale &scale);
	void DrawScaleMarkValue(double y, const CString &value);

};