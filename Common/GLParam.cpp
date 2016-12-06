#include "StdAfx.h"
#include "GLParam.h"

//---------------------------- class GLParam -------------

GLParam::GLParam(void)
{
	m_eField = 0;
	m_Scale = nullptr;
	m_DrawInstrument = true;
	m_FillColor = true;

	SetAngle(10);
	SetStep3D(16);
}

GLParam::GLParam(const GLParam &parameter)
{
	m_Scale = parameter.m_Scale;
	m_PolygonBuilder = parameter.m_PolygonBuilder;
	m_ScaleVisualizer = parameter.m_ScaleVisualizer;
	m_EmptyScale = parameter.m_EmptyScale;
	m_eField = parameter.m_eField;
	m_DrawInstrument = parameter.m_DrawInstrument;
	m_FillColor = parameter.m_FillColor;
	SetAngle(10);
	SetStep3D(16);
}

GLParam& GLParam::operator = (const GLParam &parameter) {
	m_Scale = parameter.m_Scale;
	m_PolygonBuilder = parameter.m_PolygonBuilder;
	m_ScaleVisualizer = parameter.m_ScaleVisualizer;
	m_EmptyScale = parameter.m_EmptyScale;
	m_eField = parameter.m_eField;
	m_DrawInstrument = parameter.m_DrawInstrument;
	m_FillColor = parameter.m_FillColor;
	return *this;
}

C2DScale* GLParam::GetScale(){
	if (m_Scale == nullptr) {
		return &m_EmptyScale;
	}

	return m_Scale;
}

void GLParam::SetScale(C2DScale &scale){
	m_Scale = &scale;

	m_PolygonBuilder.ClearDisplayList();
	m_ScaleVisualizer.ClearDisplayList();
}

void GLParam::SetStep3D(int step3d)
{
	m_Step3D = step3d /* *3.6 */;
	if (m_Step3D <= 3)
		m_Step3D = 3;
}

void GLParam::SetAngle(DBL angle)
{
	const DBL minAngle = 1;
	m_Angle = angle * 3.6;

	if (m_Angle <= minAngle)
		m_Angle = minAngle;
	if (m_Angle > 360 - minAngle)
		m_Angle = 360;/*- minAngle;*/
}