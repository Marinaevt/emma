#include "stdafx.h"
#include "DlgNewScale.h"


CDlgNewScale::CDlgNewScale(CWnd* pParent) : CDialog(CDlgNewScale::IDD, pParent)
{
	m_Scale = nullptr;
	m_dlg_path.Empty();
}

CDlgNewScale::CDlgNewScale(CString path, CWnd* pParent) : CDialog(CDlgNewScale::IDD, pParent)
{
	m_Scale = nullptr;
	m_dlg_path = path;
}

CDlgNewScale::~CDlgNewScale(void)
{
	delete m_Scale;
}

BEGIN_MESSAGE_MAP(CDlgNewScale, CDialog)

	ON_COMMAND(IDOK, &CDlgNewScale::OnButtonSaveClick)
	ON_COMMAND(IDC_GENERATE_SCALE, &CDlgNewScale::GenerateScale)

END_MESSAGE_MAP()

void CDlgNewScale::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EDIT_SCALE_NAME, m_ScaleName);
	DDX_Control(pDX, IDC_EDIT_SCALE_COUNT, m_ScaleCount);
	DDX_Control(pDX, IDC_CHECKDISCRETE, m_IsDiscreteChk);
	DDX_Control(pDX, IDC_CHECKLOG, m_IsLogChk);
}

int CDlgNewScale::OnInitDialog() {
	CDialog::OnInitDialog();

	CRect rect;

	GetDlgItem(IDC_OPENGL_WINDOW)->GetWindowRect(rect);

	ScreenToClient(rect);

	return 1;	// return 1 unless you set the focus to a control
				// EXCEPTION: OCX Property Pages should return 0
}

void CDlgNewScale::OnButtonSaveClick()
{

	if (!m_Scale) {
		GenerateScale();
	}

	bool bSaveScale = SaveScale();

	EndDialog(IDOK);
}

bool CDlgNewScale::SaveScale()
{
	bool bSave = m_Scale->SaveToPath(m_dlg_path);
	return bSave;

}

void CDlgNewScale::GenerateScale()
{
	delete m_Scale;

	m_Scale = new C2DScale();
	m_Scale->m_Name() = GetString(m_ScaleName);
	int count = _wtoi(GetString(m_ScaleCount));
	m_Scale->Set(count, m_IsDiscreteChk.GetState() == BST_CHECKED, m_IsLogChk.GetState() == BST_CHECKED);

	DrawScale();
}

CString CDlgNewScale::GetString(CEdit &control)
{
	CString str;

	int len = control.LineLength(control.LineIndex(0));
	control.GetLine(0, str.GetBuffer(len), len);
	str.ReleaseBuffer(len);

	return str;
}

CString CDlgNewScale::GetScaleName()
{
	return m_Scale->m_Name();
}

void CDlgNewScale::DrawScale()
{
	CClientDC cdc(this);
	CRect rect;
	GetDlgItem(IDC_OPENGL_WINDOW)->GetWindowRect(rect);
	ScreenToClient(rect);
	cdc.FillSolidRect(rect, GetSysColor(COLOR_3DFACE));								// Заливка области рисования стандартным цветом

	CPen pen(PS_SOLID, 1, RGB(0, 0, 0));
	cdc.SelectObject(&pen);

	DBL dDistribution = 0.0; 
	if (m_Scale->GetCount() > 1) {
		//не делим на 0
		dDistribution = static_cast<double>(rect.Width()) / (m_Scale->GetCount() - 1) + 0.5;	// + 0.5 для правильного округления
	}
	
	int cellWidth = m_Scale->IsLog() ? 0 : static_cast<int>(dDistribution);
	int right = 0;
	int left = rect.left;

	for (size_t i = 0; i < m_Scale->GetCount() - 1; i++)
	{
		if (m_Scale->IsLog())
		{
			left += cellWidth;
			cellWidth = static_cast<int>((m_Scale->GetValue(i + 1) - m_Scale->GetValue(i)) * rect.Width());
		}
		else {
			left = rect.left + i * cellWidth;
		}

		if (left >= rect.right)
			break;

		right = (left + cellWidth > rect.right) ? rect.right : left + cellWidth;

		CRect cell = CRect(left, rect.top, right, rect.bottom);						// Ячейка шкалы

		if (!m_Scale->IsDescrete())
			DrawGradientRect(cell, m_Scale->GetColor(m_Scale->GetValue(i)), m_Scale->GetColor(m_Scale->GetValue(i + 1)));
		else
		{
			C2DColor fillColor = m_Scale->GetColor(m_Scale->GetValue(i));
			cdc.FillSolidRect(cell, RGB(fillColor.GetR() * 255, fillColor.GetG() * 255, fillColor.GetB() * 255));
		}

		cdc.MoveTo(left, cell.top);													// Чёрная линия слева ячейки
		cdc.LineTo(left, cell.bottom);
	}

	cdc.MoveTo(right - 1, rect.top);												// Чёрная линия справа последней ячейки
	cdc.LineTo(right - 1, rect.bottom);
}

// Да-да, корпорация зла не потрудилась научить MFC заливать градиентом
void CDlgNewScale::DrawGradientRect(CRect &r, C2DColor &cLeft, C2DColor &cRight)
{
	CClientDC cdc(this);
	CRect stepR;									// Маленький прямоугольник [n X 2]
	C2DColor fillColor;								// Цвет его зарисовки
	C2DColor col_i, col_i1;							// Вспомогательные цвета

	for (int i = 0; i < r.Width(); i++)
	{
		// Устанавливаем координаты прямоугольника
		SetRect(&stepR, r.left + i, r.top, r.left + i + 1, r.bottom);

		// Вычисляем цвета левой и правой стороны прямоугольника и берем средний цвет
		col_i = cLeft * (1.0 - (static_cast<double>(i) / r.Width())) + cRight * (static_cast<double>(i) / r.Width());
		col_i1 = cLeft * (1.0 - (static_cast<double>(i + 1) / r.Width())) + cRight * (static_cast<double>(i + 1) / r.Width());
		fillColor = (col_i + col_i1) / 2;

		// Заливаем прямоугольник этим цветом
		cdc.FillSolidRect(stepR, RGB(fillColor.GetR() * 255, fillColor.GetG() * 255, fillColor.GetB() * 255));
	}
}