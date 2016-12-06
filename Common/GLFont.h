#pragma once
#include <gl/GL.h>
#include "SimpleMath.h"
#include "DlgShowError.h"
#include "ResourceCommon.h"

class CGLFont
{
public:
	CGLFont(void);
	~CGLFont(void);
	void BuildFont(LOGFONT *lf, HDC hDC);
	void KillFont();

	void OutText2D(wchar_t *text, double x, double y, double H, double angle = 0);
	double GetTextWidth(wchar_t *text, double H);
	void TranslateFromUnicod(wchar_t *dst, wchar_t* src);
	void TranslateToUnicod(wchar_t *dst, wchar_t* src);
protected:
	CString m_strBase;
	bool m_bCreated;
	GLuint m_uBase;
	GLYPHMETRICSFLOAT *m_GMF;
};
