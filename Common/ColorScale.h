#pragma once

#include "SimpleMath.h"
#include "Color.h"

//---------- class C2DScaleColor ------------------
class C2DScaleColor
{
private:
	double m_Value;			//начальное значение
	double m_SecondValue;	//конечное значение
	C2DColor m_Color;		//цвет

public:
	C2DScaleColor();

	C2DScaleColor(double value, C2DColor color);

	C2DScaleColor(const C2DScaleColor& scaleValue);

	void operator = (const C2DScaleColor &scaleColor);

	double GetValue() const{
		return m_Value;
	}

	void SetValue(double value){
		m_Value = value;
	}

	double GetSecondValue() const{
		return m_SecondValue;
	}

	void SetSecondValue(double secondValue){
		m_SecondValue = secondValue;
	}

	C2DColor GetColor() const{
		return m_Color;
	}

	void SetColor(const C2DColor &color){
		m_Color = color;
	}

	// Ищет в себе значение value
	bool Exists(double value) const;

	// Возвращает строку вида i для непрерывной шкалы, "[i; i + 1)" или
	// если это последний сплит, то вида "[n - 1; n]" для дискретной
	CString GetStringValue(double minValue, double maxValue, bool onlyFirstValue) const;

};
