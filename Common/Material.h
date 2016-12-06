// Material.h: interface for the Material class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include "2DPoint.h"
#include "./Logger/Logger.h"
#include "Strings.h"

typedef struct RX_MODEL{
	double R;
	double a_drv;
	double n_drv;

	double a_drx;
	double n_drx;

	double A_rv;
	double B_rv;
	double Q_rv;
	double m_rv;

	double A_rx;
	double B_rx;
	double Q_rx;
	double m_rx;

	double A_e;
	double m_ed;
	double Q_e;
	double m_ee;

	double A_drx;
	double m_drx;
	double Q_dd;	//20 doubles

//SRX
	double n_srx;

	double A_s;
	double p1;
	double p2;
	double p3;
	double Q_srx;

	double C0;
	double C1;
	double q1;
	double q2;
	double q3;
	double Q_ds;	//32 doubles

//MRX
	double n_mrx;

	double A_m;
	double Q_mz;
	double p_m;
	double Q_mrx;

	double C_m;
	double m_mrx;	//39 doubles

//GRAIN GROWTH
	double t_gg1;
	double A1;
	double B1;
	double A2;
	double B2;	//44 doubles

}RX_MODEL;

//Загрузка значений RX параметров
bool ReadRX(CStorage &file, RX_MODEL &rx);

class CMaterial 
{
public:
	
	CMaterial();
	virtual ~CMaterial();

	CMaterial(const CMaterial& src);
	const CMaterial& operator = (const CMaterial &Src);

	//Выбор способа (тип) загрузки материала из файла
	int Load(const CString &strFile);

	//Загружаем файл в строку для сохранения в IO::CMaterial
	CString LoadToString(const CString &strFile);
	//Загружаем материал из сохранённой строки (IO::CMaterial) в файл для чтения имеющимся методом ReadValue
	bool LoadFromString(const CString &strContent);

	//Выбор формулы расчёта напряжений
	double Si(double e, double v, double T, double Xd=-1, double X0 = 0) const;
	
	//Вычисление параметра mu = sigma/(3*vel_epsilon)
	double Mu(double e, double v, double T, double Xd=-1, double X0 = 0) const;
	
	//Вычисление связи K, dt, mu, параметра А
	double CalcK(DBL dt, DBL mju) const;

	//Получение коэффициента объёмного сжатия (К)
	double GetK(double dt, double mju) const;

	double Lambda(double T);
	double MD(double T);
	double SH(double T);
	double Emiss(double T);

	CString GetName(void);
	bool LoadMaterial(const CString &strDir, const CString &strFile);
	void UnLoad();
	bool IsLoaded()const { return m_nType != 0; }
	void GetDrxParameters(double e, double t, double T, double ee, double dt, double D0, double &Xd, double &Z, double &D);
	void GetSrxParameters(double e, double ee, double T, double Tdef, double Z, double D0, double Xd, double t, double td, double &Xm, double &dXs, double &D, double &e0);

protected:
	CString m_strFName;	//путь до файла
	int m_nType;
	int n_e, n_v, n_t;
	double *M_t, *M_e, **M_v, ***Mat;

	double m_A;
	double m_m1, m_m2, m_m3, m_m4;

	DBL m_dK;
	DBL m_dMaxMju;

	std::vector<Math::C2DPoint> m_Lambda;	//Conductivity	(теплопроводность)
	std::vector<Math::C2DPoint> m_MD;	//Mass Dencity	(плотность)
	std::vector<Math::C2DPoint> m_C;	//Specific Heat	(удельная теплоемкость)
	std::vector<Math::C2DPoint> m_Emiss;	//Emission	(Коэффициент черноты)

	RX_MODEL m_RX;

protected:
	double Si_t1(double e, double v, double T) const;
	double Si_t2(double e, double v, double T) const;
	//сейчас не используется, разделил на Si_t3, Si_t4
	double Si_t34(double e, double v, double T, double Xd, double X0) const;

	double Si_t3(double e, double v, double T, double Xd, double X0) const;
	double Si_t4(double e, double v, double T, double Xd, double X0) const;
	
	double Si_t5(double e, double v, double T) const;
	
	int Load_t1(CStorage &File);
	int Load_t2(CStorage &File);

	int Load_t3(CStorage &File);
	int Load_t4(CStorage &File);
	
	int Load_t5(CStorage &File);
	int LoadThermo(CStorage &File);
    void DeleteArrays();

};

