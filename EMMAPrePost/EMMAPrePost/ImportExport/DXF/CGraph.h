#pragma once

#include <iostream>
#include <fstream>
#include <vector>

#include "utility.h"

#include "CContour.h"
#include "CCurve.h"
#include "CPoint.h"

#include "CReader.h"

namespace DXF {

	class CContour;
	typedef std::vector<CContour> CContourArray;

//! Плоский чертёж.
/*!
	Содержит контуры, линии и точки.
	Умеет загружаться из файлов и сохраняться в файлы формата DXF.
*/
class CGraph {
public:
	CGraph () {}
	~CGraph ();
	
	bool LoadFromDXF (const char * path);
	bool LoadFromDXF (CReader & reader);
	
	bool SaveToDXF (const char * path);
	bool SaveToDXF (CReader & reader);

	//! Добавление контура в чертёж.
	/*!
		Записывает контур, предварительно регистрируя все его кривые и точки в чертеже. Если такой
		контур уже существует, то заново он не добавляется.
	*/
	void AddContour (CContour & contour);
	
	//! Добавление кривой в чертёж.
	/*!
		Если объект, расположенный по указателю "curve", уже существует в чертеже, то указатель
		"curve" освобождается, и возвращается указатель на существующий объект. В противном случае
		указатель "curve" записывается в массив кривых ("m_curves") и возвращается функцией.
		\return Указатель на кривую, расположенную в массиве кривых чертежа ("m_curves").
	*/
	CCurve * AddCurve (CCurve * curve);
	
	//! Добавление точки в чертёж.
	/*!
		Если объект, расположенный по указателю "point", уже существует в чертеже, то указатель
		"point" освобождается, и возвращается указатель на существующий объект. В противном случае
		указатель "point" записывается в массив точек ("m_points") и возвращается функцией.
		\return Указатель на точку, расположенную в массиве точек чертежа ("m_points").
	*/
	CPoint * AddPoint (CPoint * point);
	
	CContourArray::iterator ContoursBegin () { return m_contours.begin(); }
	CContourArray::iterator ContoursEnd () { return m_contours.end(); }
	
	CContourArray::const_iterator ContoursBegin () const { return m_contours.begin(); }
	CContourArray::const_iterator ContoursEnd () const { return m_contours.end(); }

	C2DOutline * MakeOutline () const;
	
private:
	CPointPtrArray m_points;
	CCurvePtrArray m_curves;
	CContourArray m_contours;
};

} // namespace DXF