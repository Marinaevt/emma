#include "StdAfx.h"
#include "Startup.h"


void Startup(){
	//LOGGER.Init(CString(_T("..\\..\\Logs\\startup.cpp_Startup.txt")));	

// IO part

	C2DTaskCreatorRegister();
	C2DRigidCreatorRegister();
	C2DPlaneFEMCreatorRegister();
	C2DSymmetryCreatorRegister();
	C2DPressureCreatorRegister();
	C2DMesh3CreatorRegister();

	C2DConstantFunctionCreatorRegister();
	C2DPiecewiseConstantFunctionCreatorRegister();
	C2DPieceLinearFunctionCreatorRegister();

	C2DNodeCreatorRegister();
	C2DCurveCreatorRegister();
	C2DCircleArcCreatorRegister();
	C2DContourCreatorRegister();
	C2DOutlineCreatorRegister();
}