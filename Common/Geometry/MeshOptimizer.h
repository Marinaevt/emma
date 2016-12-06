#pragma once
#include "GeometryPrimitives.h"

namespace FemLibrary
{
	FEMLIBRARY_API bool MovePointMSW_Test(const std::vector<Vector2D>& vertexes, Vector2D& new_point);
	FEMLIBRARY_API bool MovePointNewton_Test(const std::vector<Vector2D>& vertexes, Vector2D& new_point);
	FEMLIBRARY_API std::vector<double> Equation2_Test(const LineSeg& s1, const LineSeg& s2);
	FEMLIBRARY_API std::vector<double> Equation3_Test(const LineSeg& s1, const LineSeg& s2, const LineSeg& s3);
}
