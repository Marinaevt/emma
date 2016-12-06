#include "StdAfx.h"
#include "Color.h"

/////////	C2DColor

const C2DColor C2DColor::RED = C2DColor(1, 0, 0);
const C2DColor C2DColor::GREEN = C2DColor(0, 1, 0);
const C2DColor C2DColor::BLUE = C2DColor(0, 0, 1);
const C2DColor C2DColor::BLACK = C2DColor(0, 0, 0);
const C2DColor C2DColor::WHITE = C2DColor(1, 1, 1);

C2DColor::C2DColor(UNLONG color)
{
	double b = ((color << 8) >> 24) / 255.0;
	double g = ((color << 16) >> 24) / 255.0;
	double r = ((color << 24) >> 24) / 255.0;

	SetColors(r, g, b);
}

void C2DColor::SetColors(double r, double g, double b)
{
	m_R = r;
	m_G = g;
	m_B = b;
}

void C2DColor::operator = (const C2DColor &color)
{
	m_R = color.m_R;
	m_G = color.m_G;
	m_B = color.m_B;
}