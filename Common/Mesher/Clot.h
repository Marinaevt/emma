#include "Globals.h"

class CClot
{
protected:
	double *m_Weights;
	int		m_nWSize;
	int		m_nType;
	
	double m_TgAng;
	double m_dCeiling;
	double	m_dVal;
public:
	void SetConstVal(double val);
	CClot();
	~CClot();

	void LoadFromFront(CFrontier *pFront, double dTgAng, double dKff, double dCeiling);
	double GetValue(Math::C2DPoint *pDP, FRONTINF inf, CFrontier *pFront);
};
