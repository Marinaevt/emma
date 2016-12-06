#include "StdAfx.h"
#include "ColorScale.h"

C2DScaleColor::C2DScaleColor() {
	m_Value = 0;
	m_SecondValue = 0;
	m_Color = C2DColor::WHITE;
}

C2DScaleColor::C2DScaleColor(double value, C2DColor color) {
	m_Value = value;
	m_SecondValue = 0;
	m_Color = color;
}

C2DScaleColor::C2DScaleColor(const C2DScaleColor& scaleValue) {
	m_Value = scaleValue.m_Value;
	m_SecondValue = scaleValue.m_SecondValue;
	m_Color = scaleValue.m_Color;
}

void C2DScaleColor::operator = (const C2DScaleColor &scaleColor) {
	m_Color = scaleColor.m_Color;
	m_Value = scaleColor.m_Value;
	m_SecondValue = scaleColor.m_SecondValue;
}

bool C2DScaleColor::Exists(double value) const {
	if (fabs(m_SecondValue) < EPS) {
		return fabs(value - m_Value) < EPS;
	}

	if (fabs(m_SecondValue - 1.0) < EPS) {
		return (m_Value <= value && value <= m_SecondValue);
	}

	return m_Value <= value && value < m_SecondValue;
}

CString C2DScaleColor::GetStringValue(double minValue, double maxValue, bool onlyFirstValue) const {

	CString result;

	if (fabs(m_SecondValue) < EPS || onlyFirstValue) {
		result.Format(_T("%g"), minValue + m_Value * (maxValue - minValue));
	}
	else {
		CString templateString = _T("[%g - %g");
		templateString += (fabs(m_SecondValue - 1.0) < EPS) ? _T("]") : _T(")");
		result.Format(templateString, minValue + m_Value * (maxValue - minValue), minValue + m_SecondValue * (maxValue - minValue));
	}

	return result;
}

