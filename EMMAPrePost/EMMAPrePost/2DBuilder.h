#pragma once
#include "emmadoc.h"
#include "Resource.h"
#include "Gen2DDoc.h"
#include "Gen2DOperation.h"
#include "CalcStage2DView.h"
#include "RegisteredMessages.h"
#include "../../Common/2DTask.h"
#include "../../Common/2DPlaneFEM.h"
#include "../../Common/2DRigid.h"
#include "../../Common/2DSymmetry.h"
#include "../../Common/2DPressure.h"
#include "../../Common/DlgShowError.h"

//! Класс Сборки
class C2DBuilder : public CEMMADoc{
	GETTYPE (C2DBUILDER)
public:
	C2DBuilder(void);
	virtual ~C2DBuilder(void);

public:

	virtual bool IsPaneVisible(size_t nID) {return nID == IDR_2DBUILDER_PANE;}

	bool Build();	// Сборка всего проекта для расчёта в ядре
	bool PrepareTools(C2DTask* task, CGen2DOperation* pTools);	//Сохраняем инструменты и траектории
	bool SaveTask(C2DTask* task);	//Сохраняем задание в файл

	void DrawMesh(CGen2DDoc *mesh, GLParam &parameter);						// Отрисовка КЭ сетки
	void DrawTools(CGen2DOperation *operation, GLParam &parameter);			// Отрисовка Инструментов
	void DrawTrajectories(CGen2DOperation *operation, GLParam &parameter);	// Отрисовка Траекторий с Инструментами

	bool GetBoundingBox(CRect2D &rect);	// Получение границ всех объектов в Сборке

	virtual void DrawGL(GLParam &parameter);
};

IODEF(C2DBuilder)