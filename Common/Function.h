#pragma once
#include "IOInterface.h"
#include "2DPoint.h"
//-------------------------------------------------
/*
class CFunctionHeader : public CIOHeader
{
	void ciofunctionheader_func() {}
public:
	~CFunctionHeader() {}
};
*/
//-------------------------------------------------

template <class CParam, class CRet>
class CFunction : public IO::CClass
{
public:
	CFunction()
	{}
	
	virtual TYPEID GetType() const = 0;

	virtual CRet GetValue(const CParam& p) = 0;
	virtual CRet GetIntegral(const CParam& from, const CParam& to) = 0;

	virtual ~CFunction()
	{}
};
//-------------------------------------------------
typedef CFunction<DBL, DBL> C2DFunction;
//-------------------------------------------------

class C2DConstantFunction : public C2DFunction
{
	GETTYPE (C2DCONSTANTFUNC)

public:
	IO::CDouble m_value;

	C2DConstantFunction();

	virtual DBL GetValue(const DBL& p);

	virtual DBL GetIntegral(const DBL& from, const DBL& to);

	virtual ~C2DConstantFunction()
	{}
};

IODEF (C2DConstantFunction)
//-------------------------------------------------

class C2DPiecewiseConstantFunction : public C2DFunction
{
	GETTYPE (C2DPIECECONSTANTFUNC)
public:

	IO::CDoubleArray m_args;
	IO::CDoubleArray m_values;

	C2DPiecewiseConstantFunction();

	/*
	void SortAll()
	{
		//аргументы и соответстующие им значения должны быть отсортированы
	}
	*/

	virtual DBL GetValue(const DBL& p);

	virtual DBL GetIntegral(const DBL& from, const DBL& to);

	virtual ~C2DPiecewiseConstantFunction()
	{}
};

IODEF (C2DPiecewiseConstantFunction)
//-------------------------------------------------

class C2DPieceLinearFunction : public C2DFunction
{
	GETTYPE (C2DPIECELINEARFUNC)
public:

	IO::CDoubleArray m_args;
	IO::CDoubleArray m_values;

	C2DPieceLinearFunction();

	DBL GetValueOn(ptrdiff_t i, DBL p); // возвращает значение функции в точке p, которая должна находиться между аргументами i и i + 1

	DBL GetValue(const DBL& p);

	DBL GetIntegralBetween(ptrdiff_t from, DBL to); // интеграл от аргумента номер from до значения аргумента to (to должно находиться между аргументами from и from + 1)
	DBL GetIntegralBetween(DBL from, ptrdiff_t to); // интеграл от значения аргумента from до аргумента номер to (from должно находиться между аргументами to - 1 и to)
	DBL GetIntegralOn(ptrdiff_t i); // интеграл на отрезке i, i + 1
	DBL GetIntegralOn(ptrdiff_t i, DBL from, DBL to); // интеграл от значения аргумента from до значения аргумента to (from и to должны находиться между аргументами i и i + 1)

	virtual DBL GetIntegral(const DBL& from, const DBL& to);

	virtual ~C2DPieceLinearFunction()
	{}
};

IODEF (C2DPieceLinearFunction)
//-------------------------------------------------
