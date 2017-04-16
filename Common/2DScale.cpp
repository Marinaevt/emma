#include "stdafx.h"
#include "2DScale.h"


IOIMPL(C2DScale, C2DSCALE)

C2DScale::C2DScale(void) : m_IsInit(false)
{
	Set(2, false, false);
	RegisterMembers();
}

C2DScale::C2DScale(size_t count, bool isDescrete, bool isLog) : m_IsInit(true)
{
	Set(count, isDescrete, isLog);
	RegisterMembers();
}

C2DScale::C2DScale(const C2DScale &scale)
{
	m_Name() = scale.m_Name();
	m_Colors = scale.m_Colors;
	m_IsInit = scale.m_IsInit;
	m_IsDescrete = scale.m_IsDescrete;
	m_IsLog = scale.m_IsLog;
	m_DrawLines = scale.m_DrawLines;

	m_MinValue = scale.m_MinValue;
	m_MaxValue = scale.m_MaxValue;	

	RegisterMembers();
}

void C2DScale::RegisterMembers()
{
	RegisterMember(&m_Name);
	RegisterMember(&m_Colors);
	RegisterMember(&m_IsDescrete);
	RegisterMember(&m_IsLog);
}

void C2DScale::operator = (const C2DScale &scale)
{
	m_Name = scale.m_Name;
	m_Colors = scale.m_Colors;
	m_IsDescrete = scale.m_IsDescrete;
	m_IsLog = scale.m_IsLog;
	m_IsInit = scale.m_IsInit;

	m_MinValue = scale.m_MinValue;
	m_MaxValue = scale.m_MaxValue;
	m_DrawLines = scale.m_DrawLines;
	m_Callbackable = scale.m_Callbackable;
}

bool C2DScale::Load(CStorage& file)
{
	bool result = IO::CClass::Load(file);
	m_IsInit = true;

	return result;	
}

bool C2DScale::SaveToPath(CString &path)
{
	CString fileName;
	fileName += path; //Переданный путь до папки со шкалами (расчётная папка)
	fileName += m_Name();
	fileName += _T(".emma_scale");

	CFileException e;
	CStorage scaleFile;
	if (!scaleFile.Open(fileName, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeBinary | CStorage::shareDenyWrite, &e))
	{	
		CString cErrorStr;
		int nLoadStr = cErrorStr.LoadStringW(ID_ERROR_2DSCALE_FILE_OPEN_FAIL);	//"File could not be opened: "
		CDlgShowError cError(cErrorStr + AllToString(e.m_cause));
		return false;
	}
	bool bResult = IO::CClass::Save(scaleFile);
	scaleFile.Close();

	return bResult;
}

void C2DScale::Set(size_t count, bool isDescrete, bool isLog)
{
	m_IsDescrete = isDescrete;
	m_IsLog = isLog;
	m_DrawLines = true;
	m_Callbackable = nullptr;

	m_Colors().clear();

	if (count > 2){
		if (IsLog())		
			for (size_t i = 1; i < count - 1; i++)
			{
				double value = log10((double)(count + 1 + 10 * i) / (count + 1));
				m_Colors().push_back(C2DScaleColor(value, CalculateDefaultColor(value)));
			}		
		else		
			for (size_t i = 0; i < count; i++)
			{
				double value = (1. / (count - 1)) * i;
				m_Colors().push_back(C2DScaleColor(value, CalculateDefaultColor(value)));
			}
	}else{
		m_Colors().insert(m_Colors().begin(), C2DScaleColor(0, C2DColor::BLUE));	//вставляем в 0 позицию
		m_Colors().push_back(C2DScaleColor(1, C2DColor::RED));
	}

	ArrangeColors();
	
	m_MinValue = 0;
	m_MaxValue = 1;
}

bool C2DScale::IsDescrete() const
{
	return m_IsDescrete();
}

bool C2DScale::IsLog() const
{
	return m_IsLog();
}

bool C2DScale::IsInit() const
{
	return m_IsInit;
}

bool C2DScale::IsDrawLines() const
{
	return m_DrawLines;
}

void C2DScale::SetDrawLines(bool draw)
{
	m_DrawLines = draw;
}

void C2DScale::AddMark(double value)
{
	if (value < 0 || value > 1)
		return;

	for (size_t i = 0; i < m_Colors().size(); i++){
		if (m_Colors[i].GetValue() < value)
		{
			m_Colors().insert(m_Colors().begin() + i + 1, C2DScaleColor(value, C2DColor::WHITE));
			return;
		}
		else if (fabs(m_Colors[i].GetValue() - value) < EPS) {
			return;
		}
	}
}

void C2DScale::SetCallback(ICallbackable *callbackable)
{
	m_Callbackable = callbackable;
}

void C2DScale::SetBounds(double minValue, double maxValue)
{
	m_MinValue = minValue;
	m_MaxValue = maxValue;

	//TODO: разобраться, зачем коллбэки тут
	if (m_Callbackable != nullptr) {
		m_Callbackable->Callback();
	}
}

DBL C2DScale::GetValue(size_t index) const
{ 
	return m_Colors[index].GetValue(); 
}

size_t C2DScale::GetCount() const
{
	return m_Colors().size(); 
}

C2DColor C2DScale::GetColor(double value) const
{
	if (value < 0) {
		return m_Colors[0].GetColor();
	}
	if (value > 1) {
		return m_Colors[m_Colors.GetSize() - 1].GetColor();
	}

	return GetAppropriateColor(value);
}

C2DColor C2DScale::GetColor(size_t index) const
{
	C2DColor tmpColor = m_Colors[index].GetColor();
	return tmpColor;
}

void C2DScale::SetColor(size_t index, const C2DColor &color)
{
	m_Colors[index].SetColor(color);
}

C2DColor C2DScale::GetAppropriateColor(double value) const
{
	if (IsDescrete()) {
		for (size_t i = 0; i < m_Colors().size(); i++)
		{
			if (m_Colors[i].Exists(value)) {
					return m_Colors[i].GetColor();			
			}
		}
	}
	else {
		for (size_t i = 0; i < m_Colors().size() - 1; i++)
		{
			if (m_Colors[i].GetValue() <= value && value <= m_Colors[i + 1].GetValue()){
				double firstValue = m_Colors[i].GetValue();
				double secondValue = m_Colors[i + 1].GetValue();
				double percentOfSecond = (value - firstValue) / (secondValue - firstValue);
				C2DColor firstColor = m_Colors[i].GetColor();
				C2DColor secondColor = m_Colors[i + 1].GetColor();

				double r = ((1 - percentOfSecond) * firstColor.GetR() + percentOfSecond * secondColor.GetR());
				double g = ((1 - percentOfSecond) * firstColor.GetG() + percentOfSecond * secondColor.GetG());
				double b = ((1 - percentOfSecond) * firstColor.GetB() + percentOfSecond * secondColor.GetB());

				return C2DColor(r, g, b);
			}
		}
	}

	return CalculateDefaultColor(value);
}

C2DColor C2DScale::CalculateDefaultColor(double value) const
{
	double N = 4;
	double gr = 0.9;
	double r, g, b;

	r = g = b = gr*(1 - min(value*N, 1));

	if(value<1./N){
		b += (value)*N;
	}else if(value<2./N){
		b += 1 - (value-1./N)*N;
		g += (value-1./N)*N;
	}else if(value<3./N){
		g += 1 - (value-2./N)*N;
		r += (value-2./N)*N;
	}else if(value<4./N){
		r += 1;
		g += (value-3./N)*N;
	}
	
	return C2DColor(r, g, b);
}

bool C2DScale::IsScaleMark(double value) const
{
	if (fabs(value) < EPS || fabs(value - 1) < EPS) //0 или 1, тогда...
	{
		return true;
	}

	for (size_t i = 0; i < m_Colors().size(); i++){
		if (fabs(m_Colors[i].GetValue() - value) < EPS){
			return true;
		}
	}

	return false;
}

void C2DScale::ArrangeColors()
{
	if (!m_IsDescrete) return;

	for (size_t i = 0; i < m_Colors().size() - 1; i++){
		m_Colors[i].SetSecondValue(m_Colors[i + 1].GetValue());
	}

	m_Colors[m_Colors().size() - 1].SetSecondValue(1);
}

CString C2DScale::GetStringValue(size_t i, bool onlyFirstValue) const
{
	CString str = m_Colors[i].GetStringValue(m_MinValue, m_MaxValue, onlyFirstValue);
	return str;
}