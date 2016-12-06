#pragma once
#include "AllTypes.h"

//---------- class C2DColor ----------------
// Не имеет Alpha. Почему? - Потому что нигде не используется
class C2DColor
{
private:
	double m_R;
	double m_G;
	double m_B;

public:
	C2DColor()
	{
		SetColors(0, 0, 0);
	}

	C2DColor(double r, double g, double b)
	{
		SetColors(r, g, b);
	}

	C2DColor(UNLONG color);

	C2DColor(const C2DColor& color)
	{
		SetColors(color.m_R, color.m_G, color.m_B);
	}

	void SetColors(double r, double g, double b);

	const double GetR()
	{
		return m_R;
	}

	void SetR(double r)
	{
		m_R = r;
	}

	const double GetG()
	{
		return m_G;
	}

	void SetG(double g)
	{
		m_G = g;
	}

	const double GetB()
	{
		return m_B;
	}

	void SetB(double b)
	{
		m_B = b;
	}

	UNLONG ToDword()
	{
		return RGB(static_cast<UNLONG>(m_R * 255), static_cast<UNLONG>(m_G * 255), static_cast<UNLONG>(m_B * 255));
		//return ((long)(m_R * 255) << 16) | ((long)(m_G * 255) << 8) | (long)(m_B * 255);
	}

	double ToDouble()
	{
		return m_R + m_G*0.1 + m_B*0.01;
	}

	void operator = (const C2DColor &color);

	C2DColor operator - (const C2DColor &color) const
	{
		return C2DColor(m_R - color.m_R, m_G - color.m_G, m_B - color.m_B);
	}

	C2DColor operator + (const C2DColor &color) const
	{
		return C2DColor(m_R + color.m_R, m_G + color.m_G, m_B + color.m_B);
	}

	C2DColor operator / (int value) const
	{
		return C2DColor(m_R / value, m_G / value, m_B / value);
	}

	C2DColor operator * (int value) const
	{
		return C2DColor(m_R * value, m_G * value, m_B * value);
	}

	C2DColor operator / (double value) const
	{
		return C2DColor(m_R / value, m_G / value, m_B / value);
	}

	C2DColor operator * (double value) const
	{
		return C2DColor(m_R * value, m_G * value, m_B * value);
	}

	static const C2DColor RED;
	static const C2DColor GREEN;
	static const C2DColor BLUE;
	static const C2DColor BLACK;
	static const C2DColor WHITE;
};
