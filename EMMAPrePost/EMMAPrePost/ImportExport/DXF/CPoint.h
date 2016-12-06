#pragma once

#include <vector>
#include <cmath>
#include "../../../../Common/2DOutline.h"
#include "CReader.h"

namespace DXF {

struct CPoint;
typedef std::vector<CPoint *> CPointPtrArray;

}

namespace DXF {

//! Точка.
/*!
	Просто точка, набор координат.
*/
struct CPoint {
	CPoint (double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
	
	bool LoadFromDXF (CReader & reader);
	bool SaveToDXF (CReader & reader);
	
	void set (double new_x, double new_y, double new_z) { this->x = new_x; this->y = new_y; this->z = new_z; }
	
	bool operator == (const CPoint & point);
	bool operator != (const CPoint & point) { return !(*this == point); }

	int ExportToOutline (C2DOutline * outline) const {
		C2DNode * node = dynamic_cast<C2DNode *>(IO::CreateObject(C2DNODE));
		
		node->SetPoint(Math::C2DPoint(x, y));
		
		return outline->AddNode(node);
	} //раскомментировать

	double x;
	double y;
	double z;
};

CPoint makeCPoint (double x = 0, double y = 0, double z = 0);

}