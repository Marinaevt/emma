// Material.cpp: implementation of the Material class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Material.h"

#ifdef _DEBUG
#define new DEBUG_NEW	//DEBUG_NEW keeps track of filename and line for each object it allocates
#endif

//IOIMPL(CMaterial, CMATERIAL)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMaterial::CMaterial() {
	m_strFName.Empty();
	n_e = n_v = n_t = 0;
	M_t = nullptr;
	M_e = nullptr;
	M_v = nullptr;
	Mat = nullptr;
	m_A = 0.0;
	m_m1 = m_m2 = m_m3 = m_m4 = 0.0;
	m_RX = { 0.0, 0.0, 0.0, 
			0.0, 0.0, 
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0,	//20 doubles
		0.0,
		0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,	//32 doubles
		0.0,
		0.0,0.0,0.0,0.0,
		0.0,0.0,	//39
		0.0,0.0,0.0,0.0,0.0
	};
	m_nType = 0;

	m_dK = 10000;
	m_dMaxMju = 1.0E+10;
	CStorage file;
	if (file.Open(_T("init.txt"), CStorage::modeRead | CStorage::typeBinary)) {
		ReadValue(file, "<K>", &m_dK);
		ReadValue(file, "<MaxMju>", &m_dMaxMju);
		file.Close();
	}
}

///////////////////////////////////////////////////////////
CMaterial::~CMaterial() {
	DeleteArrays();
}

CMaterial::CMaterial(const CMaterial& Src) {

	n_e = Src.n_e;
	n_v = Src.n_v;
	n_t = Src.n_t;

	M_e = new double[n_e];
	M_t = new double[n_t];
	M_v = new double *[n_t];
	Mat = new double **[n_t];

	for (int i = 0; i<n_t; i++) {
		M_v[i] = new double[n_v];
		Mat[i] = new double*[n_v];
		for (int j = 0; j<n_v; j++) {
			Mat[i][j] = new double[n_e];
		}
	}

	for (int t = 0; t<n_t; t++) {
		M_t[t] = Src.M_t[t];
		for (int v = 0; v<n_v; v++) {
			M_v[t][v] = Src.M_v[t][v];
			for (int e = 0; e<n_e; e++) {
				Mat[t][v][e] = Src.Mat[t][v][e];
			}
		}
	}

	for (int e = 0; e<n_e; e++) {
		M_e[e] = Src.M_e[e];
	}

	m_nType = Src.m_nType;
	m_strFName = Src.m_strFName;
	m_dK = Src.m_dK;
	m_dMaxMju = Src.m_dMaxMju;
	m_A = Src.m_A;
	m_m1 = Src.m_m1;
	m_m2 = Src.m_m2;
	m_m3 = Src.m_m3;
	m_m4 = Src.m_m4;

	m_Lambda = Src.m_Lambda;
	m_MD = Src.m_MD;
	m_C = Src.m_C;
	m_Emiss = Src.m_Emiss;
	m_RX = Src.m_RX;
}

///////////////////////////////////////////////////////////
const CMaterial &CMaterial::operator = (const CMaterial &Src) {
	if (&Src == this) return *this;

	if (n_t != Src.n_t || n_e != Src.n_e || n_v != Src.n_v) {
		DeleteArrays();

		n_e = Src.n_e;
		n_v = Src.n_v;
		n_t = Src.n_t;

		M_e = new double[n_e];
		M_t = new double[n_t];
		M_v = new double *[n_t];
		Mat = new double **[n_t];

		for (int i = 0; i<n_t; i++) {
			M_v[i] = new double[n_v];
			Mat[i] = new double*[n_v];
			for (int j = 0; j<n_v; j++) {
				Mat[i][j] = new double[n_e];
			}
		}
	}
	for (int t = 0; t<n_t; t++) {
		M_t[t] = Src.M_t[t];
		for (int v = 0; v<n_v; v++) {
			M_v[t][v] = Src.M_v[t][v];
			for (int e = 0; e<n_e; e++) {
				Mat[t][v][e] = Src.Mat[t][v][e];
			}
		}
	}
	for (int e = 0; e<n_e; e++) {
		M_e[e] = Src.M_e[e];
	}

	m_nType = Src.m_nType;
	m_strFName = Src.m_strFName;
	m_dK = Src.m_dK;
	m_dMaxMju = Src.m_dMaxMju;
	m_A = Src.m_A;
	m_m1 = Src.m_m1;
	m_m2 = Src.m_m2;
	m_m3 = Src.m_m3;
	m_m4 = Src.m_m4;

	m_Lambda = Src.m_Lambda;
	m_MD = Src.m_MD;
	m_C = Src.m_C;
	m_Emiss = Src.m_Emiss;
	m_RX = Src.m_RX;
	return *this;
}

double CMaterial::Si_t1(double e_, double v_, double t_)  const
{

	//проверка на тип материала
	if (m_nType != 1) {
		CDlgShowError cError(ID_ERROR_MATERIAL_TYPE1_WRONG); //_T("m_nType 1 is wrong"));
		return 0.0;
	}

	double Se;
	double e0 = m_m4 / (m_m2 - 1);
	if (e_ >= e0) {
		Se = pow(e_, m_m2)*exp(m_m4 / e_);
	}
	else {
		Se = e_*pow(e0, m_m2)*exp(m_m4 / e0) / e0;
	}
	double dSi = m_A*exp(m_m1*t_)*Se*pow(v_, m_m3);	//K*[exp^(m1*t)]*[Se]*[ee^m]
	return dSi;
}

///////////////////////////////////////////////////////////
double CMaterial::Si_t2(double e_, double v_, double t_)  const
{
	//проверка на тип материала
	if (m_nType != 2) {
		CDlgShowError cError(ID_ERROR_MATERIAL_TYPE2_WRONG);
		return 0.0;
	}

	int i1, j1, ki, k, i, j, i_tm, i_v, i_e;
	double si_, v0, vv, ee, temper, v_min0, v_max0, sig0, sig1;
	double Sig00[2], Sig000[2];

	si_ = 0.0;

	ee = e_;
	if (ee < M_e[0]) ee = M_e[0];
	if (ee > M_e[n_e - 1]) ee = M_e[n_e - 1];

	temper = t_;
	if (temper < M_t[0]) temper = M_t[0];
	if (temper > M_t[n_t - 1]) temper = M_t[n_t - 1];

	i1 = 1;
	if (n_t == 1) i1 = 0;
	j1 = 1;
	if (n_v == 1) j1 = 0;
	ki = 1;
	if (n_e == 1) ki = 0;

	//позиция температуры
	i_tm = -1;
	for (k = 0; k<n_t - 1; k++) {
		if (temper < M_t[k]) {
			i_tm = k;
			break;
		}
	}
	if (i_tm == -1) {
		i_tm = n_t - 1;
		i1 = 0;
	}

	//min, max скорости деформации при Т, T+1
	v_min0 = Math::APPR(M_t[i_tm], M_v[i_tm][0], M_t[i_tm + i1], M_v[i_tm + i1][0], temper);
	v_max0 = Math::APPR(M_t[i_tm], M_v[i_tm][n_v - 1], M_t[i_tm + i1], M_v[i_tm + i1][n_v - 1], temper);

	v0 = v_;
	if (v0 < v_min0) v0 = v_min0;
	else  if (v0 > v_max0)  v0 = v_max0;

	for (i = 0; i <= i1; i++)   // i_tm, i_tm+i1
	{
		vv = Math::APPR(M_v[i_tm + i][0], v_min0, M_v[i_tm + i][n_v - 1], v_max0, v0);
		//позиция скорости
		i_v = 0;
		for (j = 0; j < n_v - 1; j++) {
			if (vv >= M_v[i_tm + i][j]) i_v = j;
			else break;
		}
		//позиция деформации
		i_e = 0;
		for (k = 0; k < n_e - 1; k++) {
			if (ee >= M_e[k]) i_e = k;
			else break;
		}

		for (j = 0; j <= j1; j++) {
			sig0 = Mat[i_tm + i][i_v + j][i_e];
			sig1 = Mat[i_tm + i][i_v + j][i_e + ki];
			Sig00[j] = Math::APPR(M_e[i_e], sig0, M_e[i_e + ki], sig1, ee);
		}
		Sig000[i] = Math::APPR(M_v[i_tm + i][i_v], Sig00[0], M_v[i_tm + i][i_v + j1], Sig00[j1], vv);
	} // for i

	si_ = Math::APPR(M_t[i_tm], Sig000[0], M_t[i_tm + i1], Sig000[i1], temper);
	return si_;		//*/
}

///////////////////////////////////////////////////////////
double CMaterial::Si_t34(double e, double ee, double T, double Xd, double X0)  const
{
	//S  = (1-X)*Srv*(1-exp(-a_srv*e))^n_srv + X*Srx
	//X  = 1-exp(-ad*max(0, (e-e0)/e0)^nd)
	//Srv= A_rv*arcsinh(B_rv*[ee*exp(Q_rv/(R*T))]^m_rv)
	//Srx= A_rx*arcsinh(B_rx*[ee*exp(Q_rx/(R*T))]^m_rx)
	//e0 = A_e*[ee*exp(Q_e/(R*T))]^m_ee

	if (m_nType != 3 && m_nType != 4) {
		CDlgShowError cError(ID_ERROR_MATERIAL_TYPE34_WRONG);
		return 0.0;
	}

	//Xd = -1 by default
	//X0 = 0 by default
	if (X0 < 0 || X0 > 1) {
		CDlgShowError cError(ID_ERROR_MATERIAL_X0_WRONG);
		return 0.0;
	}

	double Srx, Srv;
	if (m_nType == 3) {
		Srx = m_RX.A_rx*Math::arcsinh(m_RX.B_rx*pow(ee*exp(m_RX.Q_rx / (m_RX.R*T)), m_RX.m_rx));
		Srv = m_RX.A_rv*Math::arcsinh(m_RX.B_rv*pow(ee*exp(m_RX.Q_rv / (m_RX.R*T)), m_RX.m_rv));
	}
	else {//type 4 - A*Z^n
		Srx = m_RX.A_rx*pow(ee*exp(m_RX.Q_rx / (m_RX.R*T)), m_RX.m_rx);
		Srv = m_RX.A_rv*pow(ee*exp(m_RX.Q_rv / (m_RX.R*T)), m_RX.m_rv);
	}

	double X;	//X_D = 1-exp(-a*((e-ec)/ec)^nd)
	if (fabs(Xd + 1) > EPS) {
		X = Xd;
	}
	else {
		double e0 = 0.0;
		/* если ee==0, то и e0==0
		if (fabs(ee) < EPS) {
		e0 = 0;
		}else{

		} //*/
		e0 = m_RX.A_e*pow(ee*exp(m_RX.Q_e / (m_RX.R*T)), m_RX.m_ee);	// e_critic
		if (e >= e0) {
			X = 1 - exp(-m_RX.a_drx*pow((e - e0) / e0, m_RX.n_drx));
		}
		else {
			X = 0.0;
		}
	}
	/* old one
	double S = Srv*pow(1-exp(-m_RX.a_drv*e), m_RX.n_drv);	//Sigma_drv = Sigma_ss_drv*(1-exp...)
	S = min(S, (1-X)*S + X*Srx);

	S = Srx*X0 + (1-X0)*S;
	//*/

	double Sigma_drv = Srv*pow(1 - exp(-m_RX.a_drv*e), m_RX.n_drv);	//Sigma_drv = Sigma_ss_drv*(1-exp...)

	double S = Sigma_drv*(1 - X) + X*Srx;

	return S;
	//return (1-X)*Srv*pow(1-exp(-m_RX.a_drv*e), m_RX.n_drv) + X*Srx;

}

///////////////////////////////////////////////////////////
double CMaterial::Si_t3(double e, double ee, double T, double Xd, double X0)  const
{
	//S  = (1-X)*Srv*(1-exp(-a_srv*e))^n_srv + X*Srx
	//X  = 1-exp(-ad*max(0, (e-e0)/e0)^nd)
	//Srv= A_rv*arcsinh(B_rv*[ee*exp(Q_rv/(R*T))]^m_rv)
	//Srx= A_rx*arcsinh(B_rx*[ee*exp(Q_rx/(R*T))]^m_rx)
	//e0 = A_e*[ee*exp(Q_e/(R*T))]^m_ee

	if (m_nType != 3) {
		CDlgShowError cError(ID_ERROR_MATERIAL_TYPE3_WRONG);
		return 0.0;
	}

	//Xd = -1 by default
	//X0 = 0 by default
	if (X0 < 0 || X0 > 1) {
		CDlgShowError cError(ID_ERROR_MATERIAL_X0_WRONG);
		return 0.0;
	}

	double Srx, Srv;
	// type 3 - 1\alpha * asinh(A*Z^n)
	Srx = m_RX.A_rx*Math::arcsinh(m_RX.B_rx*pow(ee*exp(m_RX.Q_rx / (m_RX.R*T)), m_RX.m_rx));
	Srv = m_RX.A_rv*Math::arcsinh(m_RX.B_rv*pow(ee*exp(m_RX.Q_rv / (m_RX.R*T)), m_RX.m_rv));

	double X;	//X_D = 1-exp(-a*((e-ec)/ec)^nd)
	if (fabs(Xd + 1) > EPS) {
		X = Xd;
	}
	else {

		double e0 = m_RX.A_e*pow(ee*exp(m_RX.Q_e / (m_RX.R*T)), m_RX.m_ee);	// e_critic
		if (e >= e0) {
			X = 1 - exp(-m_RX.a_drx*pow((e - e0) / e0, m_RX.n_drx));
		}
		else {
			X = 0.0;
		}
	}

	double Sigma_drv = Srv*pow(1 - exp(-m_RX.a_drv*e), m_RX.n_drv);	//Sigma_drv = Sigma_ss_drv*(1-exp...)
	double S = Sigma_drv*(1 - X) + X*Srx;

	return S;
}

///////////////////////////////////////////////////////////
double CMaterial::Si_t4(double e, double ee, double T, double Xd, double X0)  const
{
	//S  = (1-X)*Srv*(1-exp(-a_srv*e))^n_srv + X*Srx
	//X  = 1-exp(-ad*max(0, (e-e0)/e0)^nd)
	//Srv= A_rv*arcsinh(B_rv*[ee*exp(Q_rv/(R*T))]^m_rv)
	//Srx= A_rx*arcsinh(B_rx*[ee*exp(Q_rx/(R*T))]^m_rx)
	//e0 = A_e*[ee*exp(Q_e/(R*T))]^m_ee

	if (m_nType != 4) {
		CDlgShowError cError(ID_ERROR_MATERIAL_TYPE4_WRONG);
		return 0.0;
	}

	//Xd = -1 by default
	//X0 = 0 by default
	if (X0 < 0 || X0 > 1) {
		CDlgShowError cError(ID_ERROR_MATERIAL_X0_WRONG);
		return 0.0;
	}

	double Srx, Srv;
	//type 4 - A*Z^n
	Srx = m_RX.A_rx*pow(ee*exp(m_RX.Q_rx / (m_RX.R*T)), m_RX.m_rx);
	Srv = m_RX.A_rv*pow(ee*exp(m_RX.Q_rv / (m_RX.R*T)), m_RX.m_rv);

	double X;	//X_D = 1-exp(-a*((e-ec)/ec)^nd)
	if (fabs(Xd + 1) > EPS) {
		X = Xd;
	}
	else {
		double e0 = m_RX.A_e*pow(ee*exp(m_RX.Q_e / (m_RX.R*T)), m_RX.m_ee);	// e_critic
		if (e >= e0) {
			X = 1 - exp(-m_RX.a_drx*pow((e - e0) / e0, m_RX.n_drx));
		}
		else {
			X = 0.0;
		}
	}

	double Sigma_drv = Srv*pow(1 - exp(-m_RX.a_drv*e), m_RX.n_drv);	//Sigma_drv = Sigma_ss_drv*(1-exp...)
	double S = Sigma_drv*(1 - X) + X*Srx;

	return S;
}

///////////////////////////////////////////////////////////
double CMaterial::Si_t5(double e_, double v_, double t_)  const
{
	
	if (m_nType != 5) {
		CDlgShowError cError(ID_ERROR_MATERIAL_TYPE5_WRONG);
		return 0.0;
	}
	
	DBL sigma;	//Результирующее напряжение
	DBL sig_s = m_A, //Переименовываем для понятности
		sig_0 = m_m1,
		k = m_m2,
		m = m_m3,
		h = m_m4;
	//sigma = sig_s*(sig_0 + k*pow(v_,m))/(sig_s+k*pow(v_,m)) + h*e_;	//Уравнение состояния (АМг-6, Иван)
	sigma = pow(e_, h)*sig_s*(sig_0 + k*pow(v_, m)) / (sig_s + k*pow(v_, m));
	return sigma;
}

///////////////////////////////////////////////////////////
double CMaterial::Si_t6(double e_, double v_, double t_, double dGrainSize)  const
{
	//m_A = K;
	//m_m1 = m
	//m_m2 = p1
	//m_m3 = alpha1
	//m_m4 = alpha2
	//m_m5 = fa0
	//m_m6 = psi

	if (m_nType != 6) {
		CDlgShowError cError(ID_ERROR_MATERIAL_TYPE5_WRONG);
		return 0.0;
	}

	DBL sigma;	//Результирующее напряжение
	DBL  K = m_A, //Переименовываем для понятности
		m = m_m1,
		p1 = m_m2,
		alpha1 = m_m3,
		alpha2 = m_m4,
		fa0 = m_m5,
		psi = m_m6;
	//sigma = sig_s*(sig_0 + k*pow(v_,m))/(sig_s+k*pow(v_,m)) + h*e_;	//Уравнение состояния (АМг-6, Иван)
	sigma = K*pow(v_, m)*pow(dGrainSize, p1)*(1 - fa0*exp(psi*e_));
	return sigma;
}

///////////////////////////////////////////////////////////


double CMaterial::Si(double e_, double v_, double t_, double Xd, double X0) const
{
	e_ = max(1.0E-6, e_);
	v_ = max(1.0E-6, v_);
	//return 68*pow(v_, 0.175);
	switch (m_nType) {
	case 1: return Si_t1(e_, v_, t_);break;
	case 2: return Si_t2(e_, v_, t_); break;
	case 3: return Si_t3(e_, v_, t_, Xd, X0); break;	//3 - A*arcsinh(B*Z^N), see C-Mn.spf
	case 4: return Si_t4(e_, v_, t_, Xd, X0); break;	//4 - A*Z^n
	case 5: return Si_t5(e_, v_, t_); break;
	case 6: return Si_t6(e_, v_, t_, Xd); break; //Albakri
	default: CDlgShowError cError(ID_ERROR_MATERIAL_TYPE_WRONG); return 0;
	}
}

double CMaterial::CalcK(DBL dt, DBL mju) const {

	//return 0.1 * m_dK + (0.81 * m_dK * m_dK * dt) / (2.0 / 3.0 * mju + 0.9 * m_dK * dt) + 2.0 / 3.0 * mju / dt;

	//return max(0.01*m_dK*dt, m_dK*dt - (2.0/3.0)*mju);
	//return max(0.0, GetK(dt,mju)*dt - (2.0/3.0)*mju);
	return max(0.01*GetK(dt, mju)*dt, GetK(dt, mju)*dt - (2.0 / 3.0)*mju);
	/*if(mju <= 0.75*GetK(dt,mju)*dt)	return GetK(dt, mju)*dt - (2.0/3.0)*mju;
	return 0.375*GetK(dt,mju)*GetK(dt,mju)*dt*dt/mju;//*/
	//return GetK(dt, mju)*dt - (2.0/3.0)*mju;
	//return m_dMaxMju - (2.0/3.0)*mju;
}

double CMaterial::GetK(double dt, double mju) const {

	return m_dK;
	//return 2.1/3.0*m_dMaxMju/dt;
	//return max(0.7*mju/dt, 5000);
	//return 2./3.*m_dK*mju/dt;// ПЛОХО. Растет толщина
	/*if(mju <= 0.75*m_dK*dt)	return m_dK;
	return 0.375*m_dK*m_dK*dt/mju + 2./3.*mju/dt;//*/
	/*if(mju <= m_dK*dt)	return m_dK;
	return m_dK/mju*m_dK/mju*m_dK*dt*dt/3. + 2./3.*mju/dt;//*/ //Плохая сходимость
}

double CMaterial::Mu(double e, double v, double T, double Xd, double X0) const
{
	double dMju = Si(e, v, T, Xd, X0) / (3.*v);
	if (v <= 0 || dMju > m_dMaxMju) {
		return m_dMaxMju;
	}

	return dMju;
}

///////////////////////////////////////////////////////////
void CMaterial::DeleteArrays() {
	int i, j;
	//if(M_e) 
	delete[] M_e;
	if (Mat) {
		for (i = 0; i<n_t; i++) {
			for (j = 0; j<n_v; j++) {
				delete[] Mat[i][j];
			}
			delete[] Mat[i];
		}
		delete[] Mat;
	}
	if (M_v) {
		for (i = 0; i<n_t; i++) delete[] M_v[i];
		delete[] M_v;
	}
}


///////////////////////////////////////////////////////////
int CMaterial::Load_t1(CStorage &File) {
	const int wrong_file = 2;

	if (!ReadValue(File, "<A>", &m_A))	return wrong_file;
	if (!ReadValue(File, "<m1>", &m_m1))	return wrong_file;
	if (!ReadValue(File, "<m2>", &m_m2))	return wrong_file;
	if (!ReadValue(File, "<m3>", &m_m3))	return wrong_file;
	if (!ReadValue(File, "<m4>", &m_m4))	return wrong_file;

	return 0;
	//return LoadThermo(File);
}

///////////////////////////////////////////////////////////
// Разбор на основе VT6.spf
int CMaterial::Load_t2(CStorage &File) {
	int n = 0, i, j, k;
	char c;
	while (n < 5) {
		if (!File.Read(&c, 1)) break;
		if (c == '\n') n++;
	}//пропускаем первые 6 строк: <type> + (5 пустых)

	 //строка вида 10 6 6
	n_e = (int)read_d(&File, nullptr);		//Кол-во точек деформаций (10)
	n_v = (int)read_d(&File, nullptr);		//Кол-во значений скоростей деформации (6) (кол-во экспериментов)
	n_t = (int)read_d(&File, nullptr);		//Кол-во значений температур (6)

											//только положительные
	if (n_e <= 0 && n_v <= 0 && n_t <= 0) {
		return 2;
	}

	M_e = new double[n_e];
	M_t = new double[n_t];
	M_v = new double *[n_t];
	Mat = new double **[n_t];

	//"расширяем" массивы
	for (i = 0; i<n_t; i++) {
		M_v[i] = new double[n_v];
		Mat[i] = new double*[n_v];
		for (j = 0; j<n_v; j++) {
			Mat[i][j] = new double[n_e];
		}
	}

	//считываем строку деформаций
	for (j = 0; j < n_e; j++) {
		M_e[j] = read_d(&File, nullptr);
	}

	//считываем матрицу скоростей ??? из различных экспериментов
	for (j = 0; j < n_t; j++) {
		for (i = 0; i < n_v; i++) {
			M_v[j][i] = read_d(&File, nullptr);
		}
	}

	//считываем строку температуры
	for (j = 0; j < n_t; j++) {
		M_t[j] = read_d(&File, nullptr);
	}
	//конец первого (начального) блока

	//Сохраняем значений напряжений в Mat(температуры, скорости, деформ)
	//кол-во деформаций определяет кол-во "матриц"\"блоков"
	//скорость и температура определяют кол-во строк и столбцов (может наоборот) в "матрице"\"блоке"
	for (k = 0; k < n_e; k++) {
		for (j = 0; j < n_v; j++) {
			for (i = 0; i < n_t; i++) {
				Mat[i][j][k] = read_d(&File, nullptr);
			}
		}
	}
	return 0;
}

///////////////////////////////////////////////////////////
int CMaterial::Load_t3(CStorage &File) {
	if (!ReadRX(File, m_RX)) return 2;	// 1/alpha * asinh(A*Z^n)
	return 0; //LoadThermo(File);
}

///////////////////////////////////////////////////////////
int CMaterial::Load_t4(CStorage &File) {
	if (!ReadRX(File, m_RX)) return 2;	// A*Z^n
	LoadThermo(File);
	return 0;
}

///////////////////////////////////////////////////////////
int CMaterial::LoadThermo(CStorage &File) {
	int wrong_file = 2;

	if (!read_array(File, "<C>", m_C))	return wrong_file;
	if (!read_array(File, "<MD>", m_MD))	return wrong_file;
	if (!read_array(File, "<Lambda>", m_Lambda))	return wrong_file;
	if (!read_array(File, "<E>", m_Emiss))return wrong_file;

	for (size_t i = 0; i < m_C.size(); i++) { m_C[i].x += 273; }
	for (size_t i = 0; i < m_MD.size(); i++) { m_MD[i].x += 273; }
	for (size_t i = 0; i < m_Lambda.size(); i++) { m_Lambda[i].x += 273; }
	for (size_t i = 0; i < m_Emiss.size(); i++) { m_Emiss[i].x += 273; }

	return 0;
}

///////////////////////////////////////////////////////////
int CMaterial::Load_t5(CStorage &File) {
	const int wrong_file = 2;

	if (!ReadValue(File, "<Sig_s>", &m_A))	return wrong_file;
	if (!ReadValue(File, "<Sig_0>", &m_m1))	return wrong_file;
	if (!ReadValue(File, "<Kv>", &m_m2))	return wrong_file;
	if (!ReadValue(File, "<Mv>", &m_m3))	return wrong_file;
	if (!ReadValue(File, "<h>", &m_m4))	return wrong_file;

	return 0;
}
///////////////////////////////////////////////////////////
int CMaterial::Load_t6(CStorage &File) {
	const int wrong_file = 2;

	if (!ReadValue(File, "<K>", &m_A))	return wrong_file;
	if (!ReadValue(File, "<m>", &m_m1))	return wrong_file;
	if (!ReadValue(File, "<p1>", &m_m2))	return wrong_file;
	if (!ReadValue(File, "<alpha1>", &m_m3))	return wrong_file;
	if (!ReadValue(File, "<alpha2>", &m_m4))	return wrong_file;
	if (!ReadValue(File, "<fa0>", &m_m5))	return wrong_file;
	if (!ReadValue(File, "<psi>", &m_m6))	return wrong_file;

	return 0;
}
// Загружаем файл материала, возвращаем 0 -> успешно
int CMaterial::Load(const CString &strFile) {
	CStorage File;
	CFileException e;
	if (!File.Open(strFile, CStorage::modeRead | CStorage::typeBinary, &e)) {
		CString cErrorStr;
		int nLoadStr = cErrorStr.LoadStringW(ID_ERROR_MATERIAL_OPEN_FILE_FAIL);
		CDlgShowError cError(cErrorStr + AllToString(e.m_cause));
		return -1;
	}

	double dType;
	if (!ReadValue(File, "<TYPE>", &dType)) {
		CDlgShowError cError(ID_ERROR_MATERIAL_READ_TYPE_FAIL); //_T("Couldn't read Material TYPE"));
		return -1;
	};

	int ret = 0, nType = static_cast<int>(dType);
	switch (nType) {
	case 1: ret = Load_t1(File); break;		//1 - m*exp(m1*t)*S*ee^m3 (often K*(ee^m))
	case 2: ret = Load_t2(File); break;		//2 - table of values (t, e, ee, sigma)
	case 3: ret = Load_t3(File); break;		//3 - A*arcsinh(B*Z^N), see C-Mn.spf
	case 4: ret = Load_t4(File); break;		//4 - A*Z^n
	case 5: ret = Load_t5(File); break;		//5 - (e^h)*sig_s*(sig_0 + k*(ee^m))/(sig_s+k*(ee^m))
	case 6: ret = Load_t6(File); break;		//6 - Albakri
	default: ret = -1;
	}
	if (ret == 0) {
		m_nType = nType;
	}
	File.Close();
	return ret;
}

bool CMaterial::LoadFromString(const CString &strContent) {

	if(strContent.IsEmpty()) {
		return false;
	}
	m_strFName = _T("./Properties.spf");
	CStorage filetemp;

	if (filetemp.Open(m_strFName, CStorage::modeCreate | CStorage::modeWrite)){
		filetemp.WriteString(strContent);
	}
	else {
		return false;
	}

	return true;
}

CString CMaterial::LoadToString(const CString &strFile) {
	CStorage File;
	CFileException e;
	if (!File.Open(strFile, CStorage::modeRead | CStorage::typeText, &e)) {
		
		CString cErrorStr;
		int nLoadStr = cErrorStr.LoadStringW(ID_ERROR_MATERIAL_OPEN_FILE_FAIL);
		CDlgShowError cError(cErrorStr + AllToString(e.m_cause));
		
		/* Виды ошибок (types of errors)
		//CFileException::m_cause
		//CFileException::genericException   1: An unspecified error occurred.
		//CFileException::fileNotFound   2: The file could not be located.
		//CFileException::badPath   3 : All or part of the path is invalid.
		*/

		//если материал был сохранён, но путь к файлу утерян, то пересоздаём и пересохраняем новый файл ./Properties.spf в тот же каталог
		
		return _T("");
	}

	CString sRes, sLine;	//итоговая и текущая строка
	
	File.SeekToBegin();
	//Get each lines from the file
	//The CString version of this function removes the '\n' if present; the LPTSTR version does not.
	while (File.ReadString(sLine))
	{
		sRes += sLine +_T("\n");		//TODO: файлы страдают избытком \r, но он нужен для CEdit
	}

	File.Close();
	return sRes;
}

double CMaterial::Lambda(double T) {
	return Math::APPR_ex(m_Lambda, T);
}
double CMaterial::MD(double T) {
	return Math::APPR_ex(m_MD, T);
}
double CMaterial::SH(double T) {
	return Math::APPR_ex(m_C, T);
}
double CMaterial::Emiss(double T) {
	return Math::APPR_ex(m_Emiss, T);
}


CString CMaterial::GetName(void) {
	return m_strFName;
}

bool CMaterial::LoadMaterial(const CString &strDir, const CString &strFile) {
	CString StrPath = strDir + _T("\\") + strFile;
	if (Load(StrPath) == 0) {
		m_strFName = strFile;
		return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////
void CMaterial::UnLoad() {
	DeleteArrays();
	n_e = n_v = n_t = 0;
	M_t = nullptr;
	M_e = nullptr;
	M_v = nullptr;
	Mat = nullptr;
	m_nType = 0;
	m_strFName.Empty();
}
///////////////////////////////////////////////////////////
void CMaterial::GetDrxParameters(double e, double t, double T, double ee, double dt, double D0, double &Xd, double &Z, double &D) {
	if (m_nType != 3 && m_nType != 4) {
		CDlgShowError cError(ID_ERROR_MATERIAL_TYPE34_WRONG);
		return;
	}

	double a = m_RX.a_drx;
	double n = m_RX.n_drx;
	double dZ = ee*exp(m_RX.Q_mz / (T*m_RX.R));
	Z = (Z*t + dZ*dt) / (t + dt);
	double _ee = (e + ee*dt) / (t + dt);
	double ec = m_RX.A_e*pow(D0, m_RX.m_ed)*pow(_ee*exp(m_RX.Q_e / (m_RX.R*T)), m_RX.m_ee);
	double de = ee*dt;
	double dXd = (e + 0.5*de <= ec) ? 0 : (1 - Xd)*a*n*pow(e + 0.5*de - ec, n - 1) / pow(ec, n)*de;
	Xd = min(1, Xd + dXd);

	if (Xd < 0.0 || Xd > 1.0) {
		CDlgShowError cError(ID_ERROR_MATERIAL_XD_WRONG); //_T("Xd is wrong"));
		return;
	}
	
	D = m_RX.A_drx*pow(Z, m_RX.m_drx);
	D = pow(Xd, 4 / 3)*D + pow(1 - Xd, 2)*D0;
}

///////////////////////////////////////////////////////////
void CMaterial::GetSrxParameters(double e, double ee, double T, double Tdef, double Z, double D0, double Xd, double t, double dt, double &Xm, double &Xs, double &D, double &e0) {
	if (m_nType != 3 && m_nType != 4) return;
	//SRX
	double t50s = m_RX.A_s*pow(D0, m_RX.p1)*pow(e, m_RX.p2)*pow(ee, m_RX.p3)*exp(m_RX.Q_srx / (m_RX.R*T));
	double ns = m_RX.n_srx;
	double Bs = 0.693 / pow(t50s, ns);
	double dXs = (1 - Xs)*Bs*ns*pow(t + 0.5*dt, ns - 1)*dt;
	Xs = min(1, Xs + dXs);
	double Ds = m_RX.C0 + m_RX.C1*pow(D0, m_RX.q1)*pow(e, m_RX.q2)*pow(ee, m_RX.q3)*exp(-m_RX.Q_ds / (m_RX.R*Tdef));
	//MRX
	double t50m = m_RX.A_m*pow(Z, m_RX.p_m)*exp(m_RX.Q_mrx / (m_RX.R*T));
	double nm = m_RX.n_mrx;
	double Bm = 0.693 / pow(t50m, nm);
	double dXm = (1 - Xm)*Bm*nm*pow(t + 0.5*dt, nm - 1)*dt;
	Xm = min(1, Xm + dXm);
	double Dd = m_RX.A_drx*pow(Z, m_RX.m_drx);
	
	//recristallized grain size
	double Dm = m_RX.C_m*pow(Z, m_RX.m_mrx);
	double Xrex = (1 - Xd)*Xs + Xd;
	double Drex = (Xrex <= 0) ? Dm : ((1 - Xd)*Xs*Ds + Xd*(Dd + (Dm - Dd)*Xm)) / Xrex;

	if (Xs < 0.0 || Xs > 1.0) {
		CDlgShowError cError(ID_ERROR_MATERIAL_XS_WRONG);
		return;
	}

	if (Xd < 0.0 || Xd > 1.0) {
		CDlgShowError cError(ID_ERROR_MATERIAL_XD_WRONG);
		return;
	}

	if (Xm < 0.0 || Xm > 1.0) {
		CDlgShowError cError(ID_ERROR_MATERIAL_XM_WRONG);
		return;
	}

	if (Xrex < 0.0 || Xrex > 1.0) {
		CDlgShowError cError(ID_ERROR_MATERIAL_XREX_WRONG);
		return;
	}

	D = Drex*pow(Xrex, 4 / 3) + pow(1 - Xrex, 2)*D0;

	double Srv_S = pow(1 - exp(-m_RX.a_drv*e), m_RX.n_drv);
	double b = 1 - pow(Srv_S *(1 - Xs)*(1 - Xd) / (1 - Xd*(1 - Xm)), 1 / m_RX.n_drv);

	if (b <= 0) {
		CDlgShowError cError(ID_ERROR_MATERIAL_B_WRONG);
		return;
	}

	e0 = -1 / m_RX.a_drv*log(b);
	
	//grain growth
	if (Xrex>0.8) {
		double Delta;
		if (t<m_RX.t_gg1) {
			Delta = dt*pow(10, m_RX.A1 - m_RX.B1 / T);
		}
		else {
			Delta = dt*pow(10, m_RX.A2 - m_RX.B2 / T);
		}
		D = sqrt(D*D + Delta);
	}
}
//////////////////////////////////////////////////////////////////////
double CMaterial::GetGrainSize(double dGrainSize, double dStrainRate, double dTimeStep) {
	if (m_nType != 6)return dGrainSize;
	return dGrainSize + dTimeStep*dStrainRate*m_m3*exp(-m_m4*dStrainRate);
}

//////////////////////////////////////////////////////////////////////
bool ReadRX(CStorage &File, RX_MODEL &rx) {
	const int N = 48;
	char *names[N] = { "<R>",    "<a_drv>","<n_drv>","<a_drx>","<n_drx>",
		"<A_rv>", "<B_rv>", "<Q_rv>", "<m_rv>", "<m_rv>",
		"<m_rv>", "<m_rv>", "<m_rv>", "<A_rx>", "<B_rx>",
		"<Q_rx>", "<m_rx>",	"<A_e>",  "<m_ed>", "<Q_e>",
		"<m_ee>", "<n_srx>","<A_s>",  "<p1>",   "<p2>",
		"<p3>",   "<Q_srx>","<C0>",   "<C1>",   "<q1>",
		"<q2>",   "<q3>",   "<Q_ds>", "<n_mrx>","<A_m>",
		"<Q_mz>", "<p_m>",  "<Q_mrx>","<C_m>",  "<m_mrx>",
		"<A_drx>","<m_drx>","<Q_dd>", "<t_gg1>","<A1>",
		"<B1>",   "<A2>",	"<B2>" };

	double *ptrs[N] = { &rx.R,    &rx.a_drv,&rx.n_drv,&rx.a_drx,&rx.n_drx,
		&rx.A_rv, &rx.B_rv, &rx.Q_rv, &rx.m_rv, &rx.m_rv,
		&rx.m_rv, &rx.m_rv, &rx.m_rv, &rx.A_rx, &rx.B_rx,
		&rx.Q_rx, &rx.m_rx,	&rx.A_e,  &rx.m_ed, &rx.Q_e,
		&rx.m_ee, &rx.n_srx,&rx.A_s,  &rx.p1,   &rx.p2,
		&rx.p3,   &rx.Q_srx,&rx.C0,   &rx.C1,   &rx.q1,
		&rx.q2,   &rx.q3,   &rx.Q_ds, &rx.n_mrx,&rx.A_m,
		&rx.Q_mz, &rx.p_m,  &rx.Q_mrx,&rx.C_m,  &rx.m_mrx,
		&rx.A_drx,&rx.m_drx,&rx.Q_dd,&rx.t_gg1,&rx.A1,
		&rx.B1,   &rx.A2,	&rx.B2 };
	bool Valid[N];
	int nSuccess = 0;
	//Считываем и локально записываем успех\неудачу
	for (int i = 0; i<N; i++) {
		Valid[i] = ReadValue(File, names[i], ptrs[i]);
		if (Valid[i]) nSuccess++;
	}

	if (nSuccess > 0) {
		return true;	//Что-то считалось
	}
	return false;	//В файле данного типа нет нужных значений
}
