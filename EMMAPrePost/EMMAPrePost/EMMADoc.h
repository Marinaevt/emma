#pragma once
#include "ViewDescriptor.h"
#include "../../Common/IOInterface.h"
#include "RegisteredMessages.h"
#include "../../Common/GLParam.h"
#include "../../Common/3DPoint.h"

class CEMMADoc : public IO::CClass
{
	GETTYPE (CEMMADOC)
public:
	CEMMADoc(CString strName = _T(""));
	CEMMADoc(CRuntimeClass *pViewClass, size_t nID, CString strName = _T(""));
	virtual ~CEMMADoc(void);

public:
	//virtual bool IsParamNeeded() { return false; } //Не используется

	void RegisterMembers();	//Регистрирует переменные для сохранения
	virtual bool CanAdd()	 /*inl*/ {return false;}//Запрещает (false) или разрешает (true) добавление подэлементов
	virtual CEMMADoc *Add()	 /*inl*/ {return nullptr;} //Добавляет подэлемент
	bool IsRemovable();	//Удаляемы ли подэлементы
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc){return false;}
	
	virtual void DrawGL(){}
	virtual void DrawGL(GLParam &parameter){}
	virtual bool IsPaneVisible(size_t nID) /*inl*/ {return false;}//Возвращать true если панель с идентификатором ID должна отображаться для данного документа
////////////////ROMA
	virtual void FillPropList(CMFCPropertyGridCtrl *pGrid){}
	virtual void UpdateProp(double *pEval, UNLONG *pId){}	//обновляет данные документа (Doc)
	virtual void UpdatePropList(CMFCPropertyGridCtrl *pGrid) {}	//обновляет все значения таблицы
	
	virtual bool LoadBody(CStorage& file);

	virtual bool GetBoundingBox(CRect2D &rect){return false;}	// Получаем обрамляющий прямоугольник на 2D плоскости
	virtual bool GetBoundingBox(CRect3D &rect) { return false; } // Получаем обрамляющий прямоугольник в 3D пространстве
	virtual void DrawGL3D() {}
	virtual void DrawGL3D(GLParam &parameter) {}

public:
	void Remove() /*inl*/ {if(m_pParentDoc) m_pParentDoc->RemoveBranch(m_hTreeItem);} // удаляет сам себя
	void RemoveBranch(HTREEITEM hItem); // удаляет ветку, содержащую hItem
	size_t FindBranch(HTREEITEM hItem); // находит ветку, содержащую hItem
	CEMMADoc *FindDoc(HTREEITEM hItem); // находит документ, соответствующий hItem
	void	RemoveAt(size_t nItem); // удаляет ветку(подраздел)
	size_t InsertSubDoc(CEMMADoc *pDoc, ptrdiff_t nIndex = -1);
	size_t GetSubDocNum()				/*inl*/ {return m_arSubDocs().size();}	//кол-во поддерикторий
	void	SetParent(CEMMADoc *pDoc)	/*inl*/ {m_pParentDoc = pDoc;}
	CEMMADoc *SubDoc(size_t nIndex)	/*inl*/ { return dynamic_cast<CEMMADoc *>(m_arSubDocs()[nIndex]); }
	CEMMADoc *GetRoot();
	CEMMADoc *GetParent();

	void	SetName(UNINT nID)			/*inl*/ {int nLoadStr = m_strName().LoadStringW(nID);}
	void	SetName(const CString &Str) /*inl*/ {m_strName() = Str;}
	CString	GetName()					/*inl*/ { return m_strName(); }

	//Получаем путь до файла сохранения всего проекта
	CString GetFilePath() { return m_strFilePath(); }
	//Устанавливаем путь до файла сохранения всего проекта
	void SetFilePath(const CString &newPath) { m_strFilePath() = newPath; }

	HTREEITEM GetTreeItem()				/*inl*/ {return m_hTreeItem;}
	HTREEITEM AddToTree(CTreeCtrl &TreeCtrl, HTREEITEM hParent, HTREEITEM hInsertAfter, int Icon_name);
	HTREEITEM ModifyTree(CTreeCtrl &TreeCtrl, int Icon_name);

//CViewDescriptor functions
	virtual void Deactivate(){};// Переопределяется, чтобы выгрузить данные документа, когда он неактивен
	virtual void *GetInterface(size_t nInterfaceID){return nullptr;}
	virtual void SetView(CWnd *pView);
	CWnd *GetView();
	CRuntimeClass *GetViewClass();
	HMENU	GetHMenu();
	HACCEL	GetHAccel();
	UNINT	GetResID();

public:

	int ChooseIconByID(int Icons_name, int Icon_name); //выбираем номер иконки из длинной картинки

	int ShowIconNumber2D(int IDS_IconName); // Возвращает номер нужной иконки в списке для 2D Задачи
	int ShowIconNumberForming(int IDS_IconName); // Возвращает номер нужной иконки в списке для Формовки
	int ShowIconNumberRolling(int IDS_IconName); // Возвращает номер нужной иконки в списке для Прокатки
	int ShowIconNumberDragging(int IDS_IconName); // Возвращает номер нужной иконки в списке для Волочения
	int ShowIconNumber3D(int IDS_IconName); // Возвращает номер нужной иконки в списке для 3D Задачи

	void SetIconName(int name); // Вводит RESOURCE_ID иконки
	int GetIconName(); // Возвращает RESOURCE_ID иконки


protected:
	// Хендел соответствующей ветки дерева
	HTREEITEM m_hTreeItem;

	// Имя раздела
	IO::CString m_strName;

	// Путь до файла сохранения
	IO::CString m_strFilePath;

	// Родительский раздел
	CEMMADoc *m_pParentDoc;

	// Подразделы
	IO::CPtrArray<CIOEMMADOCPTRARRAY> m_arSubDocs;

protected:
	//Обновление всех представлений
	void InitialUpdateView();
	//Обновление всех панелей (с тулбарами, кнопками и прочими таблицами)
	void InitialUpdatePane();

	virtual CViewDescriptor *GetViewDescriptor(){return m_pViewDescriptor;}

private:
	int m_icon; // переменная для распознавания иконки при создании дерева
	CViewDescriptor *m_pViewDescriptor; //Информация о представлении и фрейме 
private:
	void InitViewDescriptor(CRuntimeClass *pViewClass, size_t nID);
	
};

IODEF (CEMMADoc)

enum eInterfaceIDs
{
	IID_I2DSKETCH = 1
};