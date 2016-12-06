#include "stdafx.h"
#include "PrePostStartup.h"

#include "../../Common/Startup.h"
#include "EMMADoc.h"

#include "Gen2DDoc.h"
#include "2DSpecimen.h"
#include "2DSketch.h"
#include "2DMesher.h"
#include "MaterialsList.h"
#include "MaterialDoc.h"

#include "Gen2DOperation.h"
#include "2DTool.h"
#include "2DTrajectory.h"
#include "2DBuilder.h"
#include "Results.h"
#include "2DResults.h"

#include "CalcStage.h"
#include "2DCalcStage.h"

#include "Gen3DDoc.h"

#include "./Rolling/Grooves.h"
#include "./Rolling/RollingDoc.h"
#include "./Rolling/RRollingOperation.h"
#include "./Rolling/TRollingOperation.h"

#include "./Forming/2DFormingDoc.h"
#include "./Forming/2DFormingOperation.h"

#include "./Balloon/C2DBalloonDoc.h"
#include "./Balloon/C2DBalloonOperation.h"

#include "./Dragging/2DDraggingDoc.h"
#include "./Dragging/2DDraggingOperation.h"

void PrePostStartup()
{
	Startup();

	CEMMADocCreatorRegister();

	CGen2DDocCreatorRegister();
	C2DFormingDocCreatorRegister();
	C2DFormingOperationCreatorRegister();

	C2DBalloonDocCreatorRegister();
	C2DBalloonOperationCreatorRegister();

	C2DDraggingDocCreatorRegister();
	C2DDraggingOperationCreatorRegister();

	C2DSpecimenCreatorRegister();	//2D задача: Чертёж, КЭ сетка, Материал
	C2DSketchCreatorRegister();
	C2DMesherCreatorRegister();
	CMaterialsListCreatorRegister();
	CMaterialDocCreatorRegister();

	CGen2DOperationCreatorRegister(); //2D операция: Инструмент, Траектория, Сборка
	C2DToolCreatorRegister();
	C2DTrajectoryCreatorRegister();
	C2DBuilderCreatorRegister();

	CCalcStageCreatorRegister();	//2D расчёты и результаты
	C2DCalcStageCreatorRegister();
	CResultsCreatorRegister();
	C2DResultsCreatorRegister();

	CGen3DDocCreatorRegister();

	CGroovesCreatorRegister();
	
	CRollingDocCreatorRegister();
	CRRollingOperationCreatorRegister();
	CTRollingOperationCreatorRegister();
	
}