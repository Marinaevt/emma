#include "StdAfx.h"
#include "Function.h"
//-------------------------------------------------

IOIMPL (C2DConstantFunction, C2DCONSTANTFUNC)

C2DConstantFunction::C2DConstantFunction()
{
	RegisterMember(&m_value);
}

DBL C2DConstantFunction::GetValue(const DBL& p)
{
	return m_value();
}

DBL C2DConstantFunction::GetIntegral(const DBL& from, const DBL& to)
{
	return m_value() * (to - from);
}

//-------------------------------------------------

IOIMPL (C2DPiecewiseConstantFunction, C2DPIECECONSTANTFUNC)

C2DPiecewiseConstantFunction::C2DPiecewiseConstantFunction()
{
	RegisterMember(&m_args);
	RegisterMember(&m_values);
}

DBL C2DPiecewiseConstantFunction::GetValue(const DBL& p)
{
	std::vector<DBL> &args = m_args();
	std::vector<DBL> &vals = m_values();

	size_t argsN = args.size();

	if (p <= args[0])
		return vals[0];

	for (size_t i = 1; i < argsN; i++)
		if (p < args[i])
			return vals[i - 1];

	return vals[argsN - 2];
}

DBL C2DPiecewiseConstantFunction::GetIntegral(const DBL& from, const DBL& to)
{
	std::vector<DBL> &args = m_args();
	std::vector<DBL> &vals = m_values();

	size_t argsN = args.size();

	DBL f, t;
	
	f = min(from, to);
	f = max(f, args[0]);
	
	t = max(from, to);
	t = min(t, args[argsN - 1]);
	
	for (size_t i = 1; i < argsN; i++)
	{
		if (f <= args[i])
		{
			DBL tmp = 0;
			DBL prev = f;

			for (size_t j = i; j < argsN; j++)
			{
				if (t <= args[j])
				{
					if (j == i)
						return (t - f) * vals[j - 1];
					else
						return (tmp + (t - args[j - 1]) * vals[j - 1]);
				}
				else
				{
					tmp += (args[j] - prev) * vals[j - 1];
					prev = args[j];
				}
			}
		}
	}

	return 0.0;
}
//-------------------------------------------------

IOIMPL (C2DPieceLinearFunction, C2DPIECELINEARFUNC)

C2DPieceLinearFunction::C2DPieceLinearFunction()
{
	RegisterMember(&m_args);
	RegisterMember(&m_values);
}

// возвращает значение функции в точке p, которая должна находиться между аргументами i и i + 1
DBL C2DPieceLinearFunction::GetValueOn(ptrdiff_t i, DBL p) 
{
	std::vector<DBL> &args = m_args();	
	std::vector<DBL> &vals = m_values();

	return Math::i01toAB(vals[i], vals[i + 1], Math::iABto01(args[i], args[i + 1], p));
}

DBL C2DPieceLinearFunction::GetValue(const DBL& p)
{
	std::vector<DBL> &args = m_args();
	std::vector<DBL> &vals = m_values();

	size_t argsN = args.size();

	if (p <= args[0])
		return vals[0];

	for (size_t i = 1; i < argsN; i++){
		if (p < args[i])
			return GetValueOn(i - 1, p);
	}

	return vals[argsN - 1];
}

// интеграл от аргумента номер from до значения аргумента to (to должно находиться между аргументами from и from + 1)
DBL C2DPieceLinearFunction::GetIntegralBetween(ptrdiff_t from, DBL to) 
{
	return (m_values()[from] + GetValueOn(from, to)) * 0.5 * (to - m_args()[from]);
}

DBL C2DPieceLinearFunction::GetIntegralBetween(DBL from, ptrdiff_t to) // интеграл от значения аргумента from до аргумента номер to (from должно находиться между аргументами to - 1 и to)
{
	return (GetValueOn(to - 1, from) + m_values()[to]) * 0.5 * (m_args()[to] - from);
}

// интеграл на отрезке i, i + 1 (метод трапеций)
DBL C2DPieceLinearFunction::GetIntegralOn(ptrdiff_t i) 
{
	return (m_values()[i] + m_values()[i + 1]) * 0.5 * (m_args()[i + 1] - m_args()[i]);
}

// интеграл от значения аргумента from до значения аргумента to (from и to должны находиться между аргументами i и i + 1)
DBL C2DPieceLinearFunction::GetIntegralOn(ptrdiff_t i, DBL from, DBL to) 
{
	return (GetValueOn(i, from) + GetValueOn(i, to)) * 0.5 * (to - from);
}

DBL C2DPieceLinearFunction::GetIntegral(const DBL& from, const DBL& to)
{
	std::vector<DBL> &args = m_args();

	size_t argsN = args.size();	// количество аргументов (x)

	DBL f, t;
	
	f = min(from, to);
	f = max(f, args[0]);	//не выходит за рамки интервала аргументов (с начала)
	
	t = max(from, to);
	t = min(t, args[argsN - 1]);	//не выходит за рамки интервала аргументов (с конца)

	for (size_t i = 1; i < argsN; i++)
	{
		if (f <= args[i])	//нижняя граница
		{
			DBL tmp = 0.0;

			for (size_t j = i; j < argsN; j++)	//верхняя граница
			{
				if (t <= args[j])
				{
					if (j == i)
						return GetIntegralOn(j - 1, f, t);
					else
						return tmp + GetIntegralBetween(j - 1, t);
				}
				else
				{
					if (j == i)
						tmp += GetIntegralBetween(f, j);
					else
						tmp += GetIntegralOn(j - 1);
				}
			}
		}
	}

	return 0.0;
}