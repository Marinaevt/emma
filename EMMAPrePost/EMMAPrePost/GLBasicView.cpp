// GLBasicView.cpp : implementation file
//
#include "stdafx.h"
#include "GLBasicView.h"


// CGLBasicView
BEGIN_MESSAGE_MAP(CGLBasicView, CGLBasicView_Parent)
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_PAINT()
	ON_WM_SIZE()
END_MESSAGE_MAP()
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
CGLBasicView::CGLBasicView()
{// Конструктор
	m_pDC = nullptr;
	m_hrc = nullptr;
	m_BGcolor[0]=m_BGcolor[1]=m_BGcolor[2]=m_BGcolor[3]=1;
	m_nCurType = 0;	//курсор-стрелка
}
////////////////////////////////////////////////////////////////////////////////
CGLBasicView::~CGLBasicView()
{// Деструктор
}
////////////////////////////////////////////////////////////////////////////////
// CGLBasicView message handlers
int CGLBasicView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{// Пост-конструктор
	//
	if (CGLBasicView_Parent::OnCreate(lpCreateStruct) == -1)return -1;
	if(!InitGL()) return -1;// Пытаемся создать контент OpenGL
	return 0;
}

// Пре-деструктор
void CGLBasicView::OnDestroy()
{
	//TRACE(_T("CGLBasicView::OnDestroy() Start\n"));
	if (m_pDC == nullptr || m_hrc == nullptr) {
		CDlgShowError cError(ID_ERROR_GLBASICVIEW_MPDC_NULL); //_T("No m_pDC or m_hrc !"));
	}
	
	// Перед удалением он не должен быть текущим
	::wglMakeCurrent(nullptr, nullptr);
	
	// Удаляем контекст воспроизведения
	if(m_hrc)::wglDeleteContext(m_hrc);
	
	// Удаляем контекст рабочей области	
	delete m_pDC; 

	CGLBasicView_Parent::OnDestroy();	
	//TRACE(_T("CGLBasicView::OnDestroy() End\n"));
}
////////////////////////////////////////////////////////////////////////////////
void CGLBasicView::OnPaint(){
	// Отрисовка контента окна
	CGLBasicView_Parent::OnPaint();
	//
	if(!m_hrc || !m_pDC->GetSafeHdc())return;
	if(::wglMakeCurrent(m_pDC->GetSafeHdc(), m_hrc) == 0) return;
	glClearColor(m_BGcolor[0], m_BGcolor[1], m_BGcolor[2], m_BGcolor[3]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearDepth(1.0);
	glPushMatrix();
	DrawScene();
	glPopMatrix();
	glFinish();
	SwapBuffers(wglGetCurrentDC());

	::wglMakeCurrent(nullptr, nullptr);
}
////////////////////////////////////////////////////////////////////////////////
void CGLBasicView::OnSize(UNINT nType, int cx, int cy){
	// Обработка изменения размеров окна
	if(cy > 0 && cx>0){
		SetProjection();
	}
	CGLBasicView_Parent::OnSize(nType, cx, cy);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Базовые функции конвы
int CGLBasicView::InitGL(){
	
	// Создание контента OpenGL
	m_pDC = new CClientDC(this);
	
	if (SetPixelFormat(m_pDC->GetSafeHdc())==0){
		CDlgShowError cError(ID_ERROR_GLBASICVIEW_FAIL_SETPIXELFORMAT); //_T("SetPixelFormat fail"));
		return 0;
	}

	m_hrc = ::wglCreateContext(m_pDC->GetSafeHdc()); // Создаем контекст воспроизведения
	if (m_hrc==nullptr){
		CDlgShowError cError(ID_ERROR_GLBASICVIEW_FAIL_CREATEGL); //_T("CreateGLContext fail"));
		return 0;
	}
	
	InitFont();
	return 1;
}

////////////////////////////////////////////////////////////////////////////////
void CGLBasicView::InitFont(){
	
	if (m_pDC == nullptr) {
		CDlgShowError cError(ID_ERROR_GLBASICVIEW_NO_DOCUMENT); //_T("Документ пуст"));
		return;
	}

	if(!m_hrc || !m_pDC->GetSafeHdc()) return;
	if(::wglMakeCurrent(m_pDC->GetSafeHdc(), m_hrc) == 0) return;

	LOGFONT lf = {  -12,        // Высота фонта
		0,        // Ширина фонта
		0,        // Угол отношения
		0,        // Угол наклона
		FW_THIN,      // Ширина шрифта
		0,        // Курсив
		0,        // Подчеркивание
		0,        // Перечеркивание
		DEFAULT_CHARSET,      // Идентификатор набора символов
		OUT_TT_PRECIS,      // Точность вывода
		CLIP_DEFAULT_PRECIS,    // Точность отсечения
		ANTIALIASED_QUALITY,    // Качество вывода
		FF_DONTCARE|DEFAULT_PITCH,  // Семейство и шаг
		_T("Arial")} ;//*/
//	memset(&lf, 0, sizeof(LOGFONT));       // zero out structure
//	lf.lfHeight = -20;                      // request a 12-pixel-height font
//	_tcsncpy_s(lf.lfFaceName, LF_FACESIZE, _T("Arial"), 7);                    // request a face name "Arial"
	m_Font.BuildFont(&lf, m_pDC->GetSafeHdc());
	::wglMakeCurrent(nullptr, nullptr);
}

////////////////////////////////////////////////////////////////////////////////
int CGLBasicView::SetPixelFormat(HDC hdc){
	// Создание конвы OpenGL
	// Заполняем поля структуры необходимыми значениями
	static PIXELFORMATDESCRIPTOR pfd = 	{
			sizeof(PIXELFORMATDESCRIPTOR),	// размер структуры
			1,                              // номер версии
			PFD_DRAW_TO_WINDOW   |          // поддержка вывода в окно
			PFD_SUPPORT_OPENGL |          // поддержка OpenGL
			PFD_DOUBLEBUFFER,             // двойная буферизация
			0, 0, 0, 0,            // биты аккумулятора игнорируются
			32,                    // 32-разрядный буфер глубины
			0,                     // буфер трафарета не используется
			0,                     // вспомогательный буфер не используется
			PFD_MAIN_PLANE,        // основной слой
			0,                     // зарезервирован
			0, 0, 0                // маски слоя игнорируются
	};

	int pixelFormat;
	// Поддерживает ли система необходимый формат пикселей?
	if((pixelFormat = ::ChoosePixelFormat(hdc, &pfd)) == 0){
		CDlgShowError cError(ID_ERROR_GLBASICVIEW_PIXELFORMAT_WRONG); //_T("С заданным форматом пикселей работать нельзя"));
		return 0;
	}

	if (::SetPixelFormat(hdc, pixelFormat, &pfd) == 0)
    {
		CDlgShowError cError(ID_ERROR_GLBASICVIEW_FAIL_SETPIXELFORMAT); //_T("SetPixelFormat failed"));
        return 0;
    }
	//
	return 1;
}

// Переназначать контент OpenGL
bool CGLBasicView::SetCurrentContext(void){
	if(!m_hrc)return false;
	if(::wglMakeCurrent(m_pDC->GetSafeHdc(), m_hrc) == 0) return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////
void CGLBasicView::SetViewport(void){
	CRect clientRect;
	GetClientRect(&clientRect);
	glViewport(0, 0, clientRect.Width(), clientRect.Height());
}

////////////////////////////////////////////////////////////////////////////////
void CGLBasicView::DrawScene(void){
	// Отрисовка всей сцены

}

