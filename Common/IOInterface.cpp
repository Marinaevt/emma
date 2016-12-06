#include "StdAfx.h"
#include "IOInterface.h"

//----------------------------------------------------------------
namespace IO {

//----------------------------------------------------------------
bool LoadType(CStorage& file, TYPEID& type)
{
	//LOGGER.Init(::CString(_T("..\\..\\Logs\\IOInterface_LoadType")));
	bool ret = (sizeof(TYPEID) == file.Read((void*)&type, sizeof(TYPEID)));

	if (!ret)
	{
		DLOG(::CString(_T("IO::Loadtype() can not load type")), log_error);
		return false;
	}

	DLOG(::CString(_T("IO::Loadtype() | type read: ")) + AllToString(type), log_info);
	
	return true;
}
//----------------------------------------------------------------

bool SaveType(CStorage& file, TYPEID type)
{
	try { 
		file.Write((void*)&type, sizeof(TYPEID)); 
	}
	catch(CException* pEx) 
	{
		DLOG(::CString(_T("IO::SaveType(")) + AllToString(type) + _T(") can not write type to the file"), log_error);
		pEx->Delete();
		return false;
	}

	DLOG(::CString(_T("IO::SaveType(") + AllToString(type) + _T(") | type written")), log_info);
	return true;
}
//----------------------------------------------------------------

CInterface* LoadObject(CStorage& file)
{
	return IOMANAGER.Load(file);
}
//----------------------------------------------------------------
CInterface* CreateObject(TYPEID type)
{
	if (!IOMANAGER.checkData(type))
		return nullptr;

	return IOMANAGER.getData(type)();
}
//----------------------------------------------------------------

bool CInterface::Load(CStorage& file)
{
	TYPEID type;

	if (!LoadType(file, type)) {
		DLOG(::CString(_T("IO::CInterface | can not LoadType()")), log_error);
		return false;
	}

	if (GetType() != type)
	{
		DLOG(::CString(_T("IO::CInterface | loaded type rejected")), log_error);
		return false;
	}

	DLOG(::CString(_T("IO::CInterface | loaded type accepted")), log_info);

	return true;
}
//----------------------------------------------------------------

bool CInterface::Save(CStorage& file)
{
	return SaveType(file, GetType());
}
//----------------------------------------------------------------

bool CInterface::Copy(CInterface* ptr)
{
	if (!ptr) return false;	
	return (GetType() == ptr->GetType());
}
//---------------- CClass -------------------------------------

bool CClass::Load(CStorage& file)
{
	DLOG(::CString(_T("IO::CClass<")) + AllToString(GetType()) + ::CString(_T("> | start reading ->")), log_info);

	bool ret = CInterface::Load(file);
	if (!ret)
		return false;

	ret = LoadBody(file);
	DLOG(::CString(_T("IO::CClass<")) + AllToString(GetType()) + ::CString(_T("> | start reading ->")), log_info);

	return ret;
}

bool CClass::LoadBody(CStorage& file)
{
	DLOG(::CString(_T("IO::CClass<")) + AllToString(GetType()) + ::CString(_T("> | start reading body ->")), log_info);
	bool ret = CBase::Load(file);
	DLOG(::CString(_T("IO::CClass<")) + AllToString(GetType()) + ::CString(_T("> | end reading body   <-")), log_info);
	return ret;
}

bool CClass::Save(CStorage& file)
{
	bool ret;

	DLOG(::CString(_T("IO::CClass<")) + AllToString(GetType()) + ::CString(_T("> | writing header")), log_info);

	ret = IO::CInterface::Save(file);
	if (!ret) return ret;

	DLOG(::CString(_T("IO::CClass<")) + AllToString(GetType()) + ::CString(_T("> | writing base")), log_info);
	ret = CBase::Save(file);
	DLOG(::CString(_T("IO::CClass<")) + AllToString(GetType()) + ::CString(_T("> | end writing")), log_info);
	return ret;
}

//---------------- CMaterial -----------------------------------

//IOIMPL(IO::CMaterial,CIOMATERIAL)

CMaterial::CMaterial(void)
{
	RegisterMember(&m_filename);
	RegisterMember(&m_filecontent);
	//m_material = nullptr;
	//RegisterMember(&m_material);
}

CMaterial::CMaterial(const IO::CMaterial &pMat) {
	m_filename = pMat.m_filename;
	m_filecontent = pMat.m_filecontent;
	m_material = pMat.m_material;
}

bool CMaterial::Load(CStorage& file)
{
	bool bInterface = CInterface::Load(file);
	bool bBase = CBase::Load(file);

	bool bMatLoadFromString = m_material.LoadFromString(m_filecontent());	//должна быть не пуста
	if (bMatLoadFromString) {
		m_filename() = m_material.GetName();
	}
	else {
		return false; //сразу отказ, чтобы уменьшить проверки
	}
	bool bMatLoad = (m_material.Load(m_filename()) == 0);
	// после загрузки m_filename можно загружать и сам материал
	bool ret = bInterface && bBase && (bMatLoadFromString && bMatLoad);

	return ret;
}

bool CMaterial::Save(CStorage& file)
{
	bool ret;

	ret = CInterface::Save(file);
	if (!ret) {
		return ret;
	}

	ret = CBase::Save(file);

	return ret;
}

bool CMaterial::ReLoadMaterial(const ::CString& filename)
{
	if (m_material.Load(filename) != 0) {	//0 - успешно загрузили
		return false;
	}

	//запоминаем путь до файла и его содержимое
	m_filename() = filename;
	m_filecontent() = m_material.LoadToString(filename);
	
	if (m_filecontent() == _T("")) {	//содержимое файла не должно быть пусто
		return false;
	}

	return true;
}

//---------------- CManager -----------------------------------
CInterface* CManager::Load(CStorage& file)
{
	CCreator p_creator;
	CClass* p_object;
	TYPEID type;
	bool ret;

	DLOG(::CString(_T("IO::CManager | start reading object")), log_info);

	ret = LoadType(file, type);	
	if (!ret)
	{
		DLOG(::CString(_T("IO::CManager | can't read type")), log_error);
		return nullptr;
	}

	//DLOG(::CString(_T("CIOManager | header read, m_type == ")) + AllToString(p_header->m_type), log_info);

	if (checkData(type))
	{
		p_creator = getData(type); // получаем создателя

		if (!p_creator)
			return nullptr;

		DLOG(::CString(_T("IO::CManager | creator function acquired")), log_info);

		p_object = p_creator(); // создаем объект

		if (!p_object)
			return nullptr;

		DLOG(::CString(_T("IO::CManager | object created")), log_info);

		LOGGER.IncShift();
		ret = p_object->LoadBody(file); // загружаем данные в объект
		LOGGER.DecShift();

		if (!ret)
		{
			delete p_object;
			DLOG(::CString(_T("IO::CManager | can't load object of type:")) + AllToString(type), log_error);
			return nullptr;
		}

		DLOG(::CString(_T("IO::CManager | object loaded")), log_info);
		DLOG(::CString(_T("IO::CManager | end reading object")), log_info);		

		return p_object;
	}

	DLOG(::CString(_T("IO::CManager | can't find object creator of type:")) + AllToString(type), log_error);
	return nullptr;
}
//----------------------------------------------------------------

bool CInterfaceAdapter::Load(CStorage& file)
{
	CInterface *tmp = LoadObject(file);

	if (tmp)
	{
		delete m_obj;

		m_obj = tmp;
		return true;
	}

	return false;
}
//----------------------------------------------------------------

bool CInterfaceAdapter::Save(CStorage& file)
{
	if (m_obj)
		return m_obj->Save(file);

	return false;
}
//----------------------------------------------------------------

bool CInterfaceAdapter::Copy(CInterface* ptr)
{
	if (!ptr)
		return false;

	if (!m_obj)
		m_obj = CreateObject(ptr->GetType());

	return m_obj->Copy(ptr);
}
//----------------------------------------------------------------

void CBase::RegisterMember(CInterface* member)
{
	if (!member)
		return;

	for (size_t i = 0; i < m_members.size(); i++)
	{
		if (m_members[i] == member)
		{
			DLOG(::CString(_T("IO::CBase: member #")) + AllToString(i+1) + _T("/") + AllToString(m_members.size()) + _T(" already registered"), log_error);
			return;
		}
	}
	
	m_members.push_back(member);
}
//----------------------------------------------------------------

bool CBase::Load(CStorage& file)
{
	bool ret = true;

	DLOG(::CString(_T("IO::CBase | members: ")) + AllToString(m_members.size()), log_info);
	LOGGER.IncShift();

	for (size_t i = 0; i < m_members.size(); i++)
	{
		DLOG(::CString(_T("IO::CBase | start reading member: ")) + AllToString(i), log_info);
		ret = m_members[i]->Load(file);

		if (!ret){
			return false;
		}

		DLOG(::CString(_T("IO::CBase | end reading member: ")) + AllToString(i), log_info);
	}
	LOGGER.DecShift();

	DLOG(::CString(_T("IO::CBase | end reading")), log_info);

	return ret;
}
//----------------------------------------------------------------

bool CBase::Save(CStorage& file)
{
	bool ret = true;

	DLOG(::CString(_T("IO::CBase | members: ")) + AllToString(m_members.size()), log_info);
	LOGGER.IncShift();

	for (size_t i = 0; i < m_members.size(); i++)
	{
		DLOG(::CString(_T("IO::CBase | start writing member: ")) + AllToString(i), log_info);
		ret = m_members[i]->Save(file);

		if (!ret)
			return false;

		DLOG(::CString(_T("IO::CBase | end writing member: ")) + AllToString(i), log_info);
	}
	LOGGER.DecShift();

	DLOG(::CString(_T("IO::CBase | end writing")), log_info);

	return ret;
}
//----------------------------------------------------------------

bool CBase::Copy(CInterface* ptr)
{
	
	//if (!IO::CInterface::Copy(ptr))
	//	return false;

	CBase *right_ptr = dynamic_cast<CBase*>(ptr);

	size_t N = min(m_members.size(), right_ptr->m_members.size());
	
	for (size_t i = 0; i < N; i++)
		if (!m_members[i]->Copy(right_ptr->m_members[i]))
			return false;
	
	return true;
}

//---------------- CString -------------------------------------

bool CString::Load(CStorage& file)
{
	DLOG(::CString(_T("IO::CString: start reading")), log_info);

	bool ret = IO::CInterface::Load(file);

	if (!ret)
		return false;

	IO::CInt count;
	
	ret = count.Load(file);

	if (!ret)
	{
		DLOG(::CString(_T("IO::CString: count of elements can't be loaded completely")), log_error);
		return false;
	}

	 DLOG(::CString(_T("IO::CString | count read:")) + AllToString(count()), log_info);

	if (count() < 0)
	{
		DLOG(::CString(_T("IO::CString: loaded count of elements < 0")), log_error);
		return false;
	}

	wchar_t* buf = new wchar_t[count() + 2];
	buf[count()] = buf[count() + 1] = 0;

	if (!count())
	{
		m_data = buf;
		return true;
	}
		
	ret = ((count() * sizeof(wchar_t)) == file.Read((void*)buf, count() * sizeof(wchar_t)));

	if (!ret)
		DLOG(::CString(_T("IO::CString: chars can't be loaded completely")), log_error);
	else
		DLOG(::CString(_T("IO::CString | bytes read: ")) + AllToString(count()), log_info);

	 DLOG(::CString(_T("IO::CString: end reading")), log_info);

	m_data = buf;

	delete [] (wchar_t*) buf;

	return ret;	
}
//----------------------------------------------------------------

bool CString::Save(CStorage& file)
{
	 DLOG(::CString(_T("IO::CString: start writing")), log_info);

	bool ret = IO::CInterface::Save(file);

	if (!ret)
		return false;

	//m_data.FreeExtra();

	IO::CInt count;
	count() = m_data.GetAllocLength();
	DLOG(::CString(_T("IO::CString | GetLength = ")) + AllToString(count()), log_info);

	ret = count.Save(file);

	if (!ret)
		return false;

	DLOG(::CString(_T("IO::CString | count written: ")) + AllToString(count()), log_info);

	if (!count())
		return true;

	try { 
		file.Write((void*)(LPCWSTR)m_data, count() * sizeof(wchar_t));
	}
	catch(CException* pEx) {
		DLOG(::CString(_T("IO::CString: chars can't be saved completely")), log_error);
		pEx->Delete();
		return false; 
	}

	DLOG(::CString(_T("IO::CString | bytes written: ")) + AllToString(count()), log_info);
	DLOG(::CString(_T("IO::CString: end writing")), log_info);

	return true;
}
//----------------------------------------------------------------

::CStringW& CString::operator ()()
{
	return m_data;
}
//----------------------------------------------------------------

const ::CStringW& CString::operator ()() const
{
	return m_data;
}

} // namespace IO