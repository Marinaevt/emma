#pragma once
#include "StdAfx.h"
#include "GLFont.h"

CGLFont::CGLFont(void) : m_bCreated(false) {
	//int nNum = m_strBase.LoadString(IDS_FONTBASE);	//TODO: временно отказался, т.к. иначе подключение из /Common/ в проект resource.h (23.06 - есть ResourceCommon.h)
	m_strBase = _T("? -0123456789+=.,!qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNMйцукенгшщзхъфывапролджэячсмитьбюёЙЦУКЕНГШЩЗХЪФЫВАПРОЛДЖЭЯЧСМИТЬБЮЁěščřžýáíéů§úĚŠČŘŽÝÁÍÉŮÚ[]{}()/%^&*@#:;""'<>");
	m_GMF = new GLYPHMETRICSFLOAT[m_strBase.GetLength()];
	m_uBase = 0;
}

CGLFont::~CGLFont(void){
	if(m_bCreated)KillFont();
	delete[] m_GMF;
}

//////////////////////////////////////////////////////////////
void CGLFont::BuildFont(LOGFONT *lf, HDC hDC){
	HFONT	font;										// Windows Font ID

	m_uBase = glGenLists(m_strBase.GetLength());								// Storage For 256 Characters

	font = CreateFontIndirect(lf);

    HFONT oldfont = static_cast<HFONT>(SelectObject(hDC, font));        // Selects The Font We Want

	for(int i=0, uBase=m_uBase; i<m_strBase.GetLength(); i++, uBase++){
		wglUseFontOutlines(	hDC,										// Select The Current DC
						m_strBase.GetAt(i),											// Starting Character
						1,										// Number Of Display Lists To Build
						uBase,										// Starting Display Lists
						0.0f,										// Deviation From The True Outlines
						0.1f,										// Font Thickness In The Z Direction
						WGL_FONT_POLYGONS,							// Use Polygons, Not Lines
						//WGL_FONT_LINES,							// Use Polygons, Not Lines
						&m_GMF[i]);										// Address Of Buffer To Recieve Data
	} //*/

	//	if(!wglUseFontBitmaps(hDC, L'у', 256, m_uBase)){				// Builds 96 Characters Starting At Character 32
/*	if(!wglUseFontOutlines(	hDC,										// Select The Current DC
						0,											// Starting Character
						255,										// Number Of Display Lists To Build
						m_uBase,										// Starting Display Lists
						0.0f,										// Deviation From The True Outlines
						0.1f,										// Font Thickness In The Z Direction
						WGL_FONT_POLYGONS,							// Use Polygons, Not Lines
						m_GMF)){										// Address Of Buffer To Recieve Data*/
/*		UNLONG error =  GetLastError();
		
	}//*/
	SelectObject(hDC, oldfont);							// Selects The Font We Want
	DeleteObject(font);	
	m_bCreated = true;
};

//////////////////////////////////////////////////////////////
void CGLFont::KillFont(){

	if (!m_bCreated) {
		CDlgShowError cError(ID_ERROR_GLFONT_MBCREATED_FALSE); //CDlgShowError cError(_T("m_bCreated is false"));
		return;
	}
	glDeleteLists(m_uBase, m_strBase.GetLength());	
	m_bCreated = false;
}

//////////////////////////////////////////////////////////////
/*void CGLFont::OutText2D(wchar_t *text, double x, double y){
	UNINT uLen = wcslen(text);
	wchar_t *List = new wchar_t[uLen];
	TranslateFromUnicod(List, text);
	glPushMatrix();
	glRasterPos2d(x,y);
	glTranslated(x , y, 0.0);					// Center Our Text On The Screen
	glPushAttrib(GL_LIST_BIT);							// Pushes The Display List Bits
	glListBase(m_uBase);									// Sets The Base Character to 0
	glCallLists(uLen, GL_UNSIGNED_SHORT, List);	// Draws The Display List Text
	glPopAttrib();										// Pops The Display List Bits
	glTranslated(-x , -y, 0.0);					// Center Our Text On The Screen
	glPopMatrix();
	delete[] List;
}*/

//////////////////////////////////////////////////////////////
void CGLFont::TranslateFromUnicod(wchar_t *dst, wchar_t* src){
	size_t len = wcslen(src);
	for(size_t i=0; i < len; i++){
		dst[i] = max(0, m_strBase.Find(src[i]));
	}
}

//////////////////////////////////////////////////////////////
void CGLFont::TranslateToUnicod(wchar_t *dst, wchar_t* src){
	size_t len = wcslen(src);
	for(size_t i=0; i < len; i++){
		dst[i] = m_strBase.GetAt(src[i]);
	}
}

//////////////////////////////////////////////////////////////
void CGLFont::OutText2D(wchar_t *text, double x, double y, double H, double angle){
	size_t uLen = wcslen(text);
	if(!uLen)return;
	
	wchar_t *List = new wchar_t[uLen];
	TranslateFromUnicod(List, text);
	double H0 = m_GMF[0].gmfBlackBoxY;
	glPushMatrix();
	glTranslated(x , y, 0.0);					// Center Our Text On The Screen
	glScaled(H/H0, H/H0, 0);
	
	if(fabs(angle) > EPS){
		glRotated(angle, 0, 0, 1);
	}
	
	glPushAttrib(GL_LIST_BIT);							// Pushes The Display List Bits
	glListBase(m_uBase);									// Sets The Base Character to 0
	glCallLists(uLen, GL_UNSIGNED_SHORT, List);	// Draws The Display List Text
	glPopAttrib();										// Pops The Display List Bits
	glPopMatrix();
	delete[] List;
}

//////////////////////////////////////////////////////////////
double CGLFont::GetTextWidth(wchar_t *text, double H){
	size_t uLen = wcslen(text);
	if(!uLen)return 0;
	
	wchar_t *List = new wchar_t[uLen];
	TranslateFromUnicod(List, text);
	double H0 = m_GMF[0].gmfBlackBoxY, W0 = m_GMF[List[0]].gmfCellIncX;
	for(size_t i=1; i<uLen; i++){
		W0 += m_GMF[List[i]].gmfCellIncX;
	}
	delete[] List;
	return H*W0/H0;
}