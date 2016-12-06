#pragma once
#include "BasicView.h"
#include "../../Common/GLFont.h"
#include "EMMAPrePost.h"

// CGLBasicView
#define	CGLBasicView_Parent CBasicView
////////////////////////////////////////////////////////////////////////////////
class CGLBasicView : public CGLBasicView_Parent{
public:
	////////////////////////////////////////////////////////////////////////////////
	CGLBasicView();// Конструктор
	virtual ~CGLBasicView();// Деструктор
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	//Пре-деструктор
	afx_msg void OnDestroy();
	// Отрисовка контента окна
	afx_msg void OnPaint();
	// Обработка изменения размеров окна
	afx_msg void OnSize(UNINT nType, int cx, int cy);
	////////////////////////////////
	int GetCursorType() { return m_nCurType; };
	void SetCursorType(int nType) { m_nCurType = nType; };

protected:
	int m_nCurType; // типы курсора (2D и 3D вид): 
					// Вид: 0 - стрелка, 
					// Вид: 1 - лупа, 2 - перемещение (рука), 3 - полет (3D), 4 - поворот (3D)
					
					// Редактор: 100 - добавление узла, 101 - отрезка, 102 - дуги, 103 - добавление скругления,
					// Редактор: 200 - привязка инструмента, 
					// КЭ сетка: 300 - точка сгущения КЭ

protected:
	CDC		*m_pDC; //контекст вывода окна
	HGLRC	m_hrc;  //контекст воспроизведения
	GLclampf	m_BGcolor[4];	//цвет фона
	CGLFont	m_Font;				//фонт

	// Базовые функции
	int InitGL();// Создание контента OpenGL
	int SetPixelFormat(HDC hdc);// Создание конвы OpenGL
	virtual void SetProjection(void){};
	virtual void DrawScene(void);// Отрисовка всей сцены
	bool SetCurrentContext(void);// Переназначить контент OpenGL
	void SetViewport(void);
	void InitFont(void);

protected:
	DECLARE_MESSAGE_MAP()
public:

};
