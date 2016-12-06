#include "stdafx.h"
#include "Concentration.h"


/////////////////////////////////////////////////////////////////////////////
//CConcentration 
/////////////////////////////////////////////////////////////////////////////
CConcentration::CConcentration(){
	m_Weights = nullptr;
	m_nWSize = 0;
	m_TgAng = 0.0;
	m_dCeiling = 0;
	m_nType = 0;
	m_dVal = 0.0;
}

CConcentration::~CConcentration(){
	delete[] m_Weights;
}

/////////////////////////////////////////////////////////////////////////////
void CConcentration::LoadFromFront(CFrontier *pFront, double dTgAng, double dKff, double dCeiling){

	m_nWSize = pFront->m_nPointsNum;
	//if(m_Weights) 
	delete[] m_Weights;
	m_Weights = new double[m_nWSize];

	for(int i=0; i<m_nWSize; i++){
		m_Weights[i] = pFront->GetWeight(i)*dKff;
	}
	m_TgAng = dTgAng;
	m_dCeiling = dCeiling;
}
/////////////////////////////////////////////////////////////////////////////
double CConcentration::GetValue(Math::C2DPoint *pDP, FRONTINF inf, CFrontier *pFront){
	if(m_nType == 1)return m_dVal;

	double dVal = m_dCeiling, dTemp, dDistance;
	double dW1;
	Math::C2DPoint DP1, DPRes;

	int i = m_nWSize;
	while(i--){
		DP1 = pFront->m_Points[i];
		dW1 = m_Weights[i];
		dDistance = pDP->Len(DP1);
		dTemp = dW1 + dDistance*m_TgAng;
		if( dTemp < dVal)
		{
			dVal = dTemp;
		}
	}
	return dVal;
}
/////////////////////////////////////////////////////////////////////////////

void CConcentration::SetConstVal(double val)
{
	m_nType = 1;
	m_dVal = val;
}
