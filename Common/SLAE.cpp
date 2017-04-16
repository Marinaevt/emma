#include "StdAfx.h"
#include "SLAE.h"


CSLAE::CSLAE(){
	m_is_null = true;
}

bool CSLAE::Init(size_t eqns, size_t band)
{
	if (!m_matr.resize(eqns, band, 0))
		return false;
	
	try {
		m_rp.resize(eqns);
		m_sol.resize(eqns);
		m_flg.resize(eqns);
	}
	catch (CException* pEx) {
		CDlgShowError cError(ID_ERROR_SLAE_INIT);	//_T("SLAE Init error"));
		pEx->Delete();
		return false;
	}
	//TRY & END_CATCH_ALL	//избавляемся от MFC макросов и вообще от MFC по максимуму

	ZeroAll();

	return true;
}

void CSLAE::ZeroAll()
{
	m_matr.reset(0.0);
	for (size_t i = 0; i < m_rp.size(); i++){
		m_rp[i] = m_sol[i] = 0;
	}
}
/*
void CSLAE::RotateMatrixLCS(size_t k, DBL alpha)
{
	DBL sina = sin(alpha),
		cosa = cos(alpha),
		cos2a = cosa * cosa,
		sin2a = sina * sina,
		sinacosa = sina * cosa;

	for (size_t i = k - m_matr.band(); i < k; i++) // rows() было band()
	{
		DBL ai1, ai2;
		ai1 = m_matr.cell(i, k) * cosa + m_matr.cell(i, k + 1) * sina;
		ai2 = m_matr.cell(i, k) * (-sina) + m_matr.cell(i, k + 1) * cosa;
		m_matr.cell(i, k) = ai1;
		m_matr.cell(i, k + 1) = ai2;
	}

	for (size_t j = k + 2; j < k + m_matr.band() + 1; j++) // rows() было band()
	{
		DBL a1j, a2j;
		a1j = m_matr.cell(k, j) * cosa + m_matr.cell(k + 1, j) * sina;
		a2j = m_matr.cell(k, j) * (-sina) + m_matr.cell(k + 1, j) * cosa;
		m_matr.cell(k, j) = a1j;
		m_matr.cell(k + 1, j) = a2j;
	}

	DBL akk = m_matr.cell(k, k) * cos2a + 2 * sinacosa * m_matr.cell(k, k + 1) + sin2a * m_matr.cell(k + 1, k + 1), 
		akk1 = (m_matr.cell(k + 1, k + 1) - m_matr.cell(k, k)) * sinacosa + m_matr.cell(k, k + 1) * (cos2a - sin2a), 
		ak1k1 = m_matr.cell(k, k) * sin2a + 2 * sinacosa * m_matr.cell(k, k + 1) + cos2a * m_matr.cell(k + 1, k + 1);

	m_matr.cell(k, k) = akk;
	m_matr.cell(k, k + 1) = akk1;
	m_matr.cell(k + 1, k + 1) = ak1k1;
}
*/

void CSLAE::RotateMatrixLCS(size_t k, DBL alpha) //k - номер узла
{
	DBL sina = sin(alpha), cosa = cos(alpha);
	for (size_t r = max(0,k- m_matr.band()); r < k; r++) // перебор элементов столбца k
	{
		size_t c = k-r; // индекс столбца k
		DBL ai1,ai2;
		ai1 = m_matr.direct_cell(r, c) * cosa - m_matr.direct_cell(r, c + 1) * sina;
		ai2 = m_matr.direct_cell(r, c) * sina + m_matr.direct_cell(r, c + 1) * cosa;
		
		m_matr.direct_cell(r, c) = ai1;
		m_matr.direct_cell(r, c + 1) = ai2;
	}
	for (size_t c = 2; c < m_matr.band(); c++)
	{
		DBL a1j, a2j;
		
		a1j = m_matr.direct_cell(k, c) * cosa - m_matr.direct_cell(k + 1, c-1) * sina;
		a2j = m_matr.direct_cell(k, c) * sina + m_matr.direct_cell(k + 1, c-1) * cosa;

		m_matr.direct_cell(k, c) = a1j;
		m_matr.direct_cell(k + 1, c-1) = a2j;
	}
	DBL akk = m_matr.cell(k, k) * cosa * cosa
		+ 2 * m_matr.cell(k, k + 1) * sina * cosa + m_matr.cell(k + 1, k + 1)* sina * sina;
	DBL akk1 = (m_matr.cell(k + 1, k + 1) - m_matr.cell(k, k))*sina * cosa
		+ m_matr.cell(k, k + 1) * (cosa * cosa - sina * sina);
	DBL ak1k1 = m_matr.cell(k + 1, k + 1) * cosa * cosa
		+ 2 * m_matr.cell(k, k + 1) * sina * cosa + m_matr.cell(k, k)* sina * sina;
	m_matr.cell(k, k) = akk;
	m_matr.cell(k, k + 1) = akk1;
	m_matr.cell(k + 1, k + 1) = ak1k1;
}
/*
void CSLAE::RotateMatrixLCS(size_t k, DBL alpha) //k - номер узла
{
	DBL sina = sin(alpha),
		cosa = cos(alpha);
	for (size_t i = 0; i < m_matr.band(); i++)
	{
		DBL ai1, ai2;
		if (i != k || i != k + 1)
		{
			continue;
		}
		ai1 = m_matr.cell(i, k) * cosa - m_matr.cell(i, k + 1) * sina;
		m_matr.cell(i, k) = ai1;
		ai2 = m_matr.cell(i, k) * sina + m_matr.cell(i, k + 1) * cosa;
		m_matr.cell(i, k + 1) = ai2;

	}
	for (size_t c = 1; c < m_matr.band(); c++)
	{
		DBL a1j, a2j;
		if (j != k || j != k + 1)
		{
			continue;
		}
		a1j = m_matr.cell(k, j) * cosa - m_matr.cell(k + 1, j) * sina;
		m_matr.cell(k, j) = a1j;
		a2j = m_matr.cell(k, j) * sina + m_matr.cell(k + 1, j) * cosa;
		m_matr.cell(k + 1, j) = a2j;
	}
	DBL akk = m_matr.cell(k, k) * cosa * cosa
		+ 2 * m_matr.cell(k, k + 1) * sina * cosa + m_matr.cell(k + 1, k + 1)* sina * sina;
	DBL akk1 = (m_matr.cell(k + 1, k + 1) - m_matr.cell(k, k))*sina * cosa
		+ m_matr.cell(k, k + 1) * (cosa * cosa - sina * sina);
	DBL ak1k1 = m_matr.cell(k + 1, k + 1) * cosa * cosa
		+ 2 * m_matr.cell(k, k + 1) * sina * cosa + m_matr.cell(k, k)* sina * sina;
	m_matr.cell(k, k) = akk;
	m_matr.cell(k, k + 1) = akk1;
	m_matr.cell(k + 1, k + 1) = ak1k1;
}
*/
void CSLAE::RotateRPLCS(size_t k, DBL alpha)
{
	DBL sina = sin(alpha),
		cosa = cos(alpha);

	DBL rp1 = m_rp[k] * cosa + m_rp[k + 1] * sina,
		rp2 = m_rp[k] * (-sina) + m_rp[k + 1] * cosa;

	m_rp[k] = rp1;
	m_rp[k + 1] = rp2;
}

//! внесение граничных условий для одного атомарного ГУ
void CSLAE::SetBC(size_t k, const C2DBCAtom& bc)
{
	// bc.type
	// [1] Qx = Px, Qy = Vy - симметрия Y (скорость по одной оси, 0 давление по другой)
	// [2] Qx = Vx, Qy = Py - симметрия X (скорость по одной оси, 0 давление по другой)
	// [3] Qx = Vx, Qy = Vy - прилипание (т.е. движется вместе с границей)

	switch (bc.getType())
	{
		case C2DBCAtom::symX: //Симметрия относительно оси X (y=0)
		{
			m_rp[k] += bc.getQx(); 

			m_matr.cell(k + 1, k + 1) = 1;

			for (size_t j = k + 2; j < m_matr.band(); j++)
			{
				m_rp[j] -= m_matr.cell(k + 1, j) * m_rp[k + 1];
				m_matr.cell(k + 1, j) = 0;
			}

			for (size_t i = 0; i < k + 1; i++)
			{
				m_rp[i] -= m_matr.cell(i, k + 1) * m_rp[k + 1];
				m_matr.cell(i, k + 1) = 0;
			}
			m_rp[k + 1] = bc.getQy();
			break;
		}

		case C2DBCAtom::symY: //Симметрия относительно оси Y (x=0)
		{
			m_rp[k] = bc.getQx(); 
			m_rp[k + 1] += bc.getQy(); 

			m_matr.cell(k, k) = 1;

			for (size_t j = k + 1; j < m_matr.band(); j++) // k + m_matr.band() // обнуляем строку k
			{
				m_rp[j] -= m_matr.cell(k, j) * m_rp[k];
				m_matr.cell(k, j) = 0;
			}

			for (size_t i = 0; i < k; i++) // k + m_matr.band() // правильно обнуляем столбец k
			{
				m_rp[i] -= m_matr.cell(i, k) * m_rp[k];
				m_matr.cell(i, k) = 0;
			}

			break;
		}

		case C2DBCAtom::kinematic:
		{
			m_rp[k] = bc.getQx(); //-= m_matr.cell(k, k) * bc.Qx;

			m_matr.cell(k, k) = 1;
			
			//Блок для координаты X
			
			//min(k+m_matr.band(), m_matr.rows())
			for (size_t j = k + 1; j < m_matr.band(); j++) // k + m_matr.band() // обнуляем строку k
			{
				m_rp[j] -= m_matr.cell(k, j) * m_rp[k];
				m_matr.cell(k, j) = 0;
			}

			for (size_t i = 0; i < k; i++) // k + m_matr.band() // правильно обнуляем столбец k
			{
				m_rp[i] -= m_matr.cell(i, k) * m_rp[k];
				m_matr.cell(i, k) = 0;
			}//*/

			////////////
			//Блок для координаты Y
			m_rp[k + 1] = bc.getQy(); //-= m_matr.cell(k + 1, k + 1) * bc.Qy;
			m_matr.cell(k + 1, k + 1) = 1;

			size_t band = m_matr.band();

			for (size_t j = k + 2; j < band; j++) // k + 1 + m_matr.band() // обнуляем строку k + 1
			{
				m_rp[j] -= m_matr.cell(k + 1, j) * m_rp[k + 1];	//чтобы матрица осталась такой же, вычитаем из правой части
				m_matr.cell(k + 1, j) = 0;
			}

			//сверху по строкам
			for (size_t i = 0; i < k + 1; i++) // k + 1 + m_matr.band() // правильно обнуляем столбец k + 1
			{
				m_rp[i] -= m_matr.cell(i, k + 1) * m_rp[k + 1]; //чтобы матрица осталась такой же, вычитаем из правой части
				m_matr.cell(i, k + 1) = 0;
			}
			
			break;
		}
		case C2DBCAtom::load:
		{	
			
			m_rp[k] += bc.getQx(); 
			m_rp[k + 1] += bc.getQy(); 
			break;
		}
	}
	
	//m_matr.WriteToLog(false);
}

//андер конструкшн бай Борхес
void CSLAE::Set3DBC(size_t k, const C3DBCAtom& bc)
{
	// bc.type
	// [1] Qx = Px, Qy = Vy, Qz = Pz
	// [2] Qx = Px, Qy = Py, Qz = Pz
	// [3] Qx = Vx, Qy = Vy, Qz = Vz - прилипание (т.е. движется вместе с границей)
	
	switch (bc.type)
	{
		case 1:
		{
			
			//m_rp[k] -= m_matr.cell(k, k) * bc.Qx;
			m_rp[k + 1] = bc.Qy; //-= m_matr.cell(k + 1, k + 1) * bc.Qy;

			//m_matr.cell(k, k) = 1;
			m_matr.cell(k + 1, k + 1) = 1;

			//for (int j = k + 1; j < m_matr.rows(); j++) // k + m_matr.band()
			//	m_matr.cell(k, j) = 0;

			for (size_t j = k + 2; j < m_matr.band(); j++) // k + 1 + m_matr.band()
			{
				m_rp[j] -= m_matr.cell(k + 1, j) * m_rp[k + 1];
				m_matr.cell(k + 1, j) = 0;
			}

			//for (int i = 0; i < k; i++) // k + m_matr.band()
			//	m_matr.cell(i, k) = 0;

			for (size_t i = 0; i < k + 1; i++) // k + 1 + m_matr.band()
			{
				m_rp[i] -= m_matr.cell(i, k + 1) * m_rp[k + 1];
				m_matr.cell(i, k + 1) = 0;
			}
			

			break;
		}

		case 2:
		{
			m_rp[k] = bc.Qx; //-= m_matr.cell(k, k) * bc.Qx;

			m_matr.cell(k, k) = 1;

			for (size_t j = k + 1; j < m_matr.band(); j++) // k + m_matr.band() // обнуляем строку k
			{
				m_rp[j] -= m_matr.cell(k, j) * m_rp[k];
				m_matr.cell(k, j) = 0;
			}

			for (size_t i = 0; i < k; i++) // k + m_matr.band() // правильно обнуляем столбец k
			{
				m_rp[i] -= m_matr.cell(i, k) * m_rp[k];
				m_matr.cell(i, k) = 0;
			}

			break;
		}

		case 3:
		{
			m_rp[k] = bc.Qx; //-= m_matr.cell(k, k) * bc.Qx;
			m_rp[k + 1] = bc.Qy; //-= m_matr.cell(k + 1, k + 1) * bc.Qy;

			m_matr.cell(k, k) = 1;
			m_matr.cell(k + 1, k + 1) = 1;

			for (size_t j = k + 1; j < m_matr.band(); j++) // k + m_matr.band() // обнуляем строку k
			{
				m_rp[j] -= m_matr.cell(k, j) * m_rp[k];
				m_matr.cell(k, j) = 0;
			}

			for (size_t j = k + 2; j < m_matr.band(); j++) // k + 1 + m_matr.band() // обнуляем строку k + 1
			{
				m_rp[j] -= m_matr.cell(k + 1, j) * m_rp[k + 1];
				m_matr.cell(k + 1, j) = 0;
			}

			for (size_t i = 0; i < k; i++) // k + m_matr.band() // правильно обнуляем столбец k
			{
				m_rp[i] -= m_matr.cell(i, k) * m_rp[k];
				m_matr.cell(i, k) = 0;
			}

			for (size_t i = 0; i < k + 1; i++) // k + 1 + m_matr.band() // правильно обнуляем столбец k + 1
			{
				m_rp[i] -= m_matr.cell(i, k + 1) * m_rp[k + 1];
				m_matr.cell(i, k + 1) = 0;
			}
			
			break;
		}
	}
}


void CSLAE::Gauss(/*int nxy2, int isl, bool bZZ*/)
{
	// nxy2 - кол-во уравнений
	// isl - ширина ленты (половины)

	size_t r, s, m, n, j;	// индексы
	double zn, anul;

	size_t nxy2 = m_matr.rows();
	size_t isl = m_matr.band();	//половина, так половина

	for( r = 0; r < nxy2; r++ )
	{
		m_rp[r] /= m_matr.direct_cell(r, 0);

		if( r == nxy2 - 1 )		//?????
			break;

		zn = m_matr.direct_cell(r, 0);	// ERROR_&_CRASH
		if (fabs(zn) < EPS) return;		// IF CRASH

		for(s = 1; s < isl; s++)
		{
			m_sol[s] = m_matr.direct_cell(r, s);
			
			if( fabs(m_sol[s]) < EPS )
				continue;
			
			m_matr.direct_cell(r, s) = m_sol[s] / zn;
		}

		for( m = 1; m < isl; m++ )
		{
			zn = m_sol[m];
			
			if( fabs(zn) < EPS )
				continue;

			n = r + m;
			
			if( n > nxy2 - 1 )
				continue;
			
			j = 0;
			for( s = m; s < isl; s++ )
			{
				anul = m_matr.direct_cell(n, j);
				m_matr.direct_cell(n, j) -= zn * m_matr.direct_cell(r, s);
				if( fabs(m_matr.direct_cell(n, 0)) < EPS ) {
					m_matr.direct_cell(n, 0) = EPS; //10^(-18)
				}
				j++;
			}
			m_rp[n] -= zn * m_rp[r];
		}
	}
	
	// цикл вычисления решения
	for( r = nxy2 - 2; r != 0; r-- ) {
		for( s = 1; s < isl; s++ ) {
			m = r + s;
			
			if( m > nxy2 - 1 )
				continue;
			
			m_rp[r] -= m_matr.direct_cell(r, s) * m_rp[m];
		}
	}

	// сохраняем решение в m_sol
	for( r = 0; r < nxy2; r++) {
		//Math::swap(m_sol[r], m_rp[r]);
		m_sol[r] = m_rp[r];
	}
}

//! Записывает в лог информацию о матрице ПРИ ОТЛАДКЕ (Debug)
void CSLAE::WriteToLogOnDebug(){

	DLOG(CString(_T("CSLAE LOG START ->")) ,log_info);
	LOGGER.IncShift();
		
		WriteToLogMatrixOnDebug();
		WriteToLogRightPartOnDebug();
		WriteToLogSolutionOnDebug();

	LOGGER.DecShift();
	DLOG(CString(_T("CSLAE LOG END   <-")) ,log_info);
}

//! Записывает в лог информацию о матрице ПРИ ОТЛАДКЕ
void CSLAE::WriteToLogMatrixOnDebug()
{
	DLOG(CString(_T("MATRIX LOG START ->")), log_info);
	LOGGER.IncShift();
	//m_matr.WriteToLog(true);			//only band matrix
	m_matr.WriteToLogFullMatrix(true);	//full 2d matrix
	LOGGER.DecShift();
	DLOG(CString(_T("MATRIX LOG END   <-")), log_info);
}

//! Записывает в лог информацию о правой части ПРИ ОТЛАДКЕ
void CSLAE::WriteToLogRightPartOnDebug()
{
	CString rpLog = _T("");

	DLOG(CString(_T("RIGHT PART LOG START ->")), log_info);
	LOGGER.IncShift();	//сдвиг блока вправо
	size_t rp_size = m_rp.size();
	for (size_t i = 0; i < rp_size; i++) {
		rpLog += AllToString(m_rp[i]) + _T(", ");
	}
	DLOG(rpLog, log_info);
	LOGGER.DecShift();	//убираем сдвиг
	DLOG(CString(_T("RIGHT PART LOG END   <-")), log_info);
}


//! Записывает в лог информацию о решении ПРИ ОТЛАДКЕ
void CSLAE::WriteToLogSolutionOnDebug()
{
	CString solLog = _T("");

	DLOG(CString(_T("SOLUTION LOG START ->")), log_info);
	LOGGER.IncShift();
	size_t sol_size = m_sol.size();
	for (size_t i = 0; i < sol_size; i++) {
		solLog += AllToString(m_sol[i]) + _T(", ");
	}
	DLOG(solLog, log_info);
	LOGGER.DecShift();
	DLOG(CString(_T("SOLUTION LOG END   <-")), log_info);
}

//! Записывает в лог информацию о матрице СЛАУ
void CSLAE::WriteToLog() {

	LOG(CString(_T("CSLAE LOG START ->")), log_info);
	LOGGER.IncShift();

	WriteToLogMatrix();
	WriteToLogRightPart();
	WriteToLogSolution();

	LOGGER.DecShift();
	LOG(CString(_T("CSLAE LOG END   <-")), log_info);
}

//! Записывает в лог информацию о матрице
void CSLAE::WriteToLogMatrix()
{
	LOG(CString(_T("MATRIX LOG START ->")), log_info);
	LOGGER.IncShift();
	//m_matr.WriteToLog(true);			//only band matrix
	m_matr.WriteToLogFullMatrix(true);	//full 2d matrix
	LOGGER.DecShift();
	LOG(CString(_T("MATRIX LOG END   <-")), log_info);
}

//! Записывает в лог информацию о правой части
void CSLAE::WriteToLogRightPart()
{
	CString rpLog = _T("");

	LOG(CString(_T("RIGHT PART LOG START ->")), log_info);
	LOGGER.IncShift();	//сдвиг блока вправо
	size_t rp_size = m_rp.size();
	for (size_t i = 0; i < rp_size; i++) {
		rpLog += AllToString(m_rp[i]) + _T(", ");
	}
	LOG(rpLog, log_info);
	LOGGER.DecShift();	//убираем сдвиг
	LOG(CString(_T("RIGHT PART LOG END   <-")), log_info);
}


//! Записывает в лог информацию о решении
void CSLAE::WriteToLogSolution()
{
	CString solLog = _T("");

	LOG(CString(_T("SOLUTION LOG START ->")), log_info);
	LOGGER.IncShift();
	size_t sol_size = m_sol.size();
	for (size_t i = 0; i < sol_size; i++) {
		solLog += AllToString(m_sol[i]) + _T(", ");
	}
	LOG(solLog, log_info);
	LOGGER.DecShift();
	LOG(CString(_T("SOLUTION LOG END   <-")), log_info);
}

CSLAE::~CSLAE()
{
}

DBL Norm2(const std::vector<DBL>& arr)
{
	DBL tmp = 0;

	for (size_t i = 0; i < arr.size(); i++) {
		tmp += arr[i] * arr[i];
	}

	return tmp;
}

DBL Length2(const std::vector<DBL>& arr1, const std::vector<DBL>& arr2)
{	
	DBL tmp0, tmp1 = 0.0;
	size_t n = min(arr1.size(), arr2.size());

	for (size_t i = 0; i < n; i++)
	{
		tmp0 = arr1[i] - arr2[i];
		tmp1 += tmp0 * tmp0;
	}

	return tmp1;
}

DBL MaxErrorP(const std::vector<DBL>& arr1, const std::vector<DBL>& arr2)
{	
	DBL maxabs = 0.0, avg = 0.0;
	size_t n = min(arr1.size(), arr2.size());

	for (size_t i = 0; i < n; i++){
		if (maxabs < fabs(arr1[i] - arr2[i])) {
			maxabs = fabs(arr1[i] - arr2[i]);
		}
		avg += fabs(arr2[i]);
	}
	avg /= n;
	if( fabs(avg) < EPS ) return 0;	//if(avg == 0) return 0;

	return maxabs/avg;
}