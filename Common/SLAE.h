#pragma once
#include "AllTypes.h"
#include "BandMatrix.h"
#include "2DBCAtom.h"
#include "3DBCAtom.h"
#include "SimpleMath.h"
#include "./Logger/Logger.h"
#include "Strings.h"
#include "ResourceCommon.h"

class CSLAE
{
public:
	BANDMATRIX			m_matr; // матрица CBandMatrix DBL
	std::vector<DBL>	m_rp;   // правая часть Ax=b
	std::vector<DBL>	m_sol;  // решение Ax=b
	std::vector<bool>	m_flg;  // флажки для контроля граничных условий	//НЕ ИСПОЛЬЗУЮТСЯ
	bool			m_is_null;	// true, если только создана и не заполнена

	CSLAE();

	bool Init(size_t eqns, size_t band);
	void ZeroAll();

	// вращает локальную систему координат, заданную переменными k (x_k) и k+1 (y_k), на угол alpha
	// для матрицы
	void RotateMatrixLCS(size_t k, DBL alpha);
	// для правой части
	void RotateRPLCS(size_t k, DBL alpha);

	// внесение граничных условий для одного атомарного ГУ
	void SetBC(size_t k, const C2DBCAtom& bc);

	//для 3d случая
	void Set3DBC(size_t k, const C3DBCAtom& bc);

	void Gauss();
	void Gauss2();
	void Gauss3();
	//запись всего (матрицы, правой части и решения)
	void WriteToLog();
	//записываем в лог саму матрицу
	void WriteToLogMatrix();
	//записываем в лог правую часть
	void WriteToLogRightPart();
	//записываем в лог решение Ax = b
	void WriteToLogSolution();

	//запись всего (матрицы, правой части и решения) ПРИ ОТЛАДКЕ
	void WriteToLogOnDebug();
	//записываем в лог саму матрицу ПРИ ОТЛАДКЕ
	void WriteToLogMatrixOnDebug();
	//записываем в лог правую часть ПРИ ОТЛАДКЕ
	void WriteToLogRightPartOnDebug();
	//записываем в лог решение Ax = b ПРИ ОТЛАДКЕ
	void WriteToLogSolutionOnDebug();

	~CSLAE();
};

DBL Norm2(const std::vector<DBL>& arr);
DBL Length2(const std::vector<DBL>& arr1, const std::vector<DBL>& arr2);
DBL MaxErrorP(const std::vector<DBL>& arr1, const std::vector<DBL>& arr2);