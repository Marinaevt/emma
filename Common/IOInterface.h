#pragma once
#include "Singleton.h"
#include "SimpleManager.h"
#include "AllTypes.h"
#include "Tchar.h"
#include "./Logger/Logger.h"
#include "Strings.h"
#include "CustomArray.h"
#include "BandMatrix.h"
#include "Material.h"
#include "2DPoint.h"

// -----------------------------------------------------------------
// Макросы для работы с иерархиями, использующими систему IO
// -----------------------------------------------------------------
// Корневой класс наследовать от IO::CClass
//
#define IOMANAGER CSingleton<IO::CManager>::Instance()

#define GETTYPE(type_id) \
	public: virtual TYPEID GetType() const {return type_id;}

// после объявления класса
#define IODEF(name) \
	IO::CClass* name##Creator(); \
	void name##CreatorRegister();

// в CPP файле класса
#define IOIMPL(name, type_id) \
	IO::CClass* name##Creator() { return new name ; } \
	void name##CreatorRegister() \
	{ IOMANAGER.registerData( type_id , name##Creator ); }


namespace IO {
	//----------------------------------------------------------------

	bool LoadType(CStorage& file, TYPEID& type);
	bool SaveType(CStorage& file, TYPEID type);
	//----------------------------------------------------------------

	//!  Класс CInterface - интерфейс всех классов, которые автоматически записывают себя в файл и считывают себя оттуда.
	/*!
	IO::CInterface - класс, от которого должны быть (напрямую или косвенно) унаследованы любые классы, для которых
	в программе предусмотрено автоматическое сохранение/загрузка.
	Виртуальная функция GetHeader, создающая и возвращающая заголовок данного объекта, должна быть определена в наследниках.
	*/
	class CInterface
	{
	public:
		CInterface(void) {}

		virtual bool Load(CStorage& file);
		virtual bool Save(CStorage& file);

		virtual TYPEID GetType() const = 0;

		virtual bool Copy(CInterface* ptr);

		virtual ~CInterface(void) {}
	};
	//----------------------------------------------------------------

	CInterface* LoadObject(CStorage& file);
	CInterface* CreateObject(TYPEID type);
	//----------------------------------------------------------------

	//****************************************************************
	// Ввод-вывод иерархии объектов
	//****************************************************************

	//!  Класс IO::CBase - класс, реализующий автоматическое сохранение и загрузку своих мемберов.
	/*!
	Класс, унаследованный от IO::CBase должен для всех своих мемберов, требующих автоматического
	ввода/вывода, вызвать в своем конструкторе функцию RegisterMember, и передать ей указатели
	на указанных мемберов.
	*/
	class CBase : public CInterface
	{
	protected:
		std::vector<CInterface*> m_members;

		void RegisterMember(CInterface* member);
	public:
		CBase(void) {}

		virtual bool Load(CStorage& file);
		virtual bool Save(CStorage& file);

		virtual bool Copy(CInterface* ptr);

		virtual ~CBase(void) {}
	};
	//----------------------------------------------------------------

	class CClass : public CBase
	{
	public:
		virtual TYPEID GetType() const = 0;

		virtual bool Load(CStorage& file);
		virtual bool LoadBody(CStorage& file);
		virtual bool Save(CStorage& file);
	};

	// Функция-создатель - создает объект конкретного класса из иерархии,
	// должна быть описана для каждого класса иерархии
	typedef CClass* (*CCreator)();

	//!  Шаблон CIOManager
	/*!
	CIOManager - шаблон фабрики объектов из некоторой иерархии классов, которая однозначно
	определяется Заголовком (CHeader), который указывается в качестве параметра шаблона.
	CIOManager хранит указатели на функции, создающие конкретные объекты иерархии.
	Все действия с указателями наследуются от шаблона SimpleManager, где в качестве
	первого параметра шаблона указывается TYPEID - число, обозначающее тип класса
	в иерархии, а вторым параметром шаблона указывается указатель на функцию-создатель
	объектов иерархии.
	Функция Load загружает объект из файла и создает его в памяти.
	Этот шаблон реализован на основе техники Singleton, т.е. для каждого Заголовка, а,
	следовательно и иерархии, существует только один CIOManager, обратиться к которому
	можно из любого места программы так: CSingleton<CIOManager<CHeader> >::Instance().
	*/
	class CManager : public SimpleManager<TYPEID, CCreator>
	{
	public:
		/*!
		Функция Load служит для загрузки объекта иерархии из файла, возвращает указатель на IO::CInterface.
		Т.к. иерархия классов, объекты которой хранятся в файле, однозначно определяется CHeader, то
		функция Load:
		<ul>
		<li> 1. вызывает функцию LoadType, которая считывает тип хранящегося в файле объекта;
		<li> 2. по считанному типу объекта находит функцию-создатель объектов такого типа;
		<li> 3. вызывает эту функцию, которая возвращает указатель типа IO::CInterface, под которым находится необходимый объект;
		<li> 4. у необходимого объекта вызывается функция Load, которая загружает его данные.
		</ul>
		*/
		CInterface* Load(CStorage& file);

		friend CSingleton<CManager>;

	private:
		CManager(void) {}
		CManager(const CManager&);
		CManager& operator=(const CManager&);
		~CManager(void) {}
	};
	//----------------------------------------------------------------

	//!  Класс CIOInterfaceAdapter - так называемый adapter или wrapper для класса IO::CInterface
	/*!
	IO::CInterfaceAdapter позволяет встраивать в статическую систему IO (на основе IO::CBase) динамические объекты,
	унаследованные от класса IO::CInterface.
	Т.е. если мы заранее не знаем точный тип объекта, принадлежащего иерархии с корнем в IO::CInterface, но хотим,
	чтобы он автоматически сохранялся и загружался посредством методов Save и Load класса IO::CBase, то надо
	пользоваться классом IO::CInterfaceAdapter в качестве типа мембера вместо IO::CInterface*.
	*/
	class CInterfaceAdapter : public CInterface
	{
	protected:
		CInterface* m_obj;

	public:
		CInterfaceAdapter() : m_obj(nullptr) {}

		GETTYPE(0)

		virtual bool Load(CStorage& file);
		virtual bool Save(CStorage& file);
		virtual bool Copy(CInterface* ptr);

		CInterfaceAdapter(const CInterfaceAdapter& obj) {
			Copy(obj.m_obj);
		}

		CInterfaceAdapter& operator = (const CInterfaceAdapter& obj) {
			Copy(obj.m_obj);
			return *this;
		}

		void SetObject(CInterface* obj)
		{
			if (obj != m_obj)
			{
				m_obj = obj;
			}
		}

		operator CInterface*() { return m_obj; }

		operator CInterface&() { return *m_obj; }

		CInterface* operator () () { return m_obj; }

		virtual ~CInterfaceAdapter() {
			delete m_obj;
		}
	};

	//----------------------------------------------------------------

	// макрос для конвертации CIOTemplate.m_data в m_data


#define DECLARE_REFERENCE( cls, obj, pref ) cls & pref##obj = obj()
#define REF_INT( obj ) DECLARE_REFERENCE (int, obj, r_)
#define REF_DBL( obj ) DECLARE_REFERENCE (DBL, obj, r_)

#define DECLARE_MEMBER(cls, obj) DECLARE_REFERENCE (cls, obj, m_)
#define MEM_INT( obj ) DECLARE_MEMBER (int, obj)
#define MEM_DBL( obj ) DECLARE_MEMBER (DBL, obj)

	//----------------------------------------------------------------

	//!  Шаблон CTemplate
	/*!
	CTemplate предназначен для быстрого создания класса, автоматически сохраняющего/загружающего
	единственного своего мембера, тип которого задается первым параметром шаблона.
	В качестве второго параметра шаблона нужно указать число, которое будет записано в файл в
	качестве типа конкретного класса.
	*/
	template <class CData, TYPEID iType>
	class CTemplate : public CInterface
	{
		CData m_data;

	public:

		CTemplate(void) {}

		virtual TYPEID GetType() const { return iType; }

		virtual bool Load(CStorage& file)
		{
			bool ret = CInterface::Load(file);

			if (ret) {
				ret = (sizeof(m_data) == file.Read((void*)&m_data, sizeof(CData)));
			}

			if (!ret) {
				DLOG(::CString(_T("IO::CTemplate<")) + AllToString(iType) + _T("> can't be loaded completely"), log_error);
			}
			else {
				DLOG(::CString(_T("IO::CTemplate<")) + AllToString(iType) + ::CString(_T("> | bytes read: ")) + AllToString(sizeof(CData)), log_info);
			}

			return ret;
		}

		virtual bool Save(CStorage& file)
		{
			bool ret = CInterface::Save(file);

			if (ret)
			{
				try{ 
					file.Write((void*)&m_data, sizeof(CData)); 
				}
				catch(CException* pEx) {
					DLOG(::CString(_T("IO::CTemplate<")) + AllToString(iType) + _T("> can't be saved completely"), log_error);
					ret = false;
					pEx->Delete();

				}
			}

			DLOG(::CString(_T("IO::CTemplate<")) + AllToString(iType) + ::CString(_T("> | 8-DBL, 9-INT, 17-BOOL | bytes written: ")) + AllToString(sizeof(CData)), log_info);

			return ret;
		}

		virtual bool Copy(CInterface* ptr)
		{
			if (CInterface::Copy(ptr))
			{
				CTemplate* right_ptr = dynamic_cast<CTemplate*>(ptr);
				m_data = right_ptr->m_data;
				return true;
			}

			return false;
		}

		CTemplate(const CTemplate& obj)
		{
			m_data = obj.m_data;
		}

		CTemplate& operator=(const CTemplate& obj)
		{
			m_data = obj.m_data;
			return *this;
		}

		//Aksenov{{
		CTemplate& operator=(const CData& data)
		{
			m_data = data;
			return *this;
		}
		//}}

		operator CData& ()
		{
			return m_data;
		}

		CData& operator () ()
		{
			return m_data;
		}

		const CData& operator () () const
		{
			return m_data;
		}

		virtual ~CTemplate() {}
	};



	//----------------------------------------------------------------

	typedef CTemplate<DBL, CIODOUBLE> CDouble;
	typedef CTemplate<int, CIOINT> CInt;
	typedef CTemplate<size_t, CIOSIZE_T> CSize_t;	//Сохраняющийся size_t
	typedef CTemplate<UNINT, CIOUNINT> CUnInt;
	typedef CTemplate<bool, CIOBOOL> CBool;
	typedef CTemplate<Math::C2DLine, CIO2DLINE> C2DLine;
	typedef CTemplate<Math::C2DPoint, CIO2DPOINT> C2DPoint;

	//----------------------------------------------------------------

	//!  Шаблон IO::CSimpleArray
	/*!
	IO::CSimpleArray предназначен для быстрого создания класса, хранящего, записывающего и читающего
	массив простых данных, тип которых указывается первым параметром шаблона.
	В качестве второго параметра шаблона нужно указать число, которое будет записано в файл в
	качестве типа конкретного класса.
	*/
	template <class CData, TYPEID iType>
	class CSimpleArray : public CInterface
	{
	protected:
		std::vector<CData> m_data;
	public:
		//std::vector<CData>::reference ref;

		/*
		Ибо в vector<bool> reference есть весчь не совсем простая. 
		http://www.cplusplus.com/reference/v...ool/reference/ И связанно это с организацией хранения информации. 
		Т.е. на самом деле возвращаемым значением для типа bool является не bool&, а объект класса vector<T>::reference. 
		Почитайте Джосатиса он об этом прекрасно рассказывает.
		*/

		CSimpleArray() {}

		GETTYPE(iType)

			virtual bool Copy(CInterface* ptr)
		{
			if (CInterface::Copy(ptr))
			{
				CSimpleArray* right_ptr = dynamic_cast<CSimpleArray*>(ptr);
				m_data = right_ptr->m_data;
				return true;
			}

			return false;
		}

		CSimpleArray<CData, iType>& operator=(const CSimpleArray<CData, iType>& arr)
		{
			m_data = arr.m_data;
			return *this;
		}

		virtual bool Load(CStorage& file)
		{
			bool ret = CInterface::Load(file);

			if (!ret)
				return false;

			CSize_t count;

			ret = count.Load(file);

			if (!ret)
			{
				DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T(">: count of elements can't be loaded completely"), log_error);
				return false;
			}

			if (count() < 0)
			{
				DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T(">: loaded count of elements < 0"), log_error);
				return false;
			}

			DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T("> | elements count read: ") + AllToString(count()), log_info);

			m_data.resize(count());

			if (!count())
				return true;

			ret = (sizeof(CData) * count() == file.Read((void*)&m_data[0], sizeof(CData) * count()));

			if (!ret)
				DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T("> can't be loaded completely"), log_error);
			else
				DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T("> | bytes read: ") + AllToString(sizeof(CData) * count()), log_info);


			return ret;
		}

		virtual bool Save(CStorage& file)
		{
			bool ret = CInterface::Save(file);

			if (!ret)
				return false;

			CSize_t count;
			count() = m_data.size();

			ret = count.Save(file);

			if (!ret) {
				DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T("> | count() not saved ") + AllToString(count()), log_info);
				return false;
			}

			DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T("> | elements count saved: ") + AllToString(count()), log_info);

			if (!count())
				return true;

			try{ 
				file.Write((void*)&m_data[0], sizeof(CData) * count()); //analog GetData()
			}
			catch(CException* pEx) {
				DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T("> can't be saved completely"), log_error);
				pEx->Delete();
				return false;
			}
				DLOG(::CString(_T("IO::CSimpleArray<")) + AllToString(iType) + _T("> | bytes written: ") + AllToString(sizeof(CData) * count()), log_info);

			return true;
		}

		operator std::vector<CData>& ()
		{
			return m_data;
		}

		std::vector<CData>& operator () ()
		{
			return m_data;
		}

		const std::vector<CData>& operator () () const
		{
			return m_data;
		}

		
		// TODO: проверить ref
		CData& operator[](size_t index)
		{
			return m_data.at(index);
		}

		const CData& operator[](size_t index) const
		{
			return m_data.at(index);
		}
		//*/

		/*
		ref operator[](size_t index)
		{
			return m_data.at(index);
		}

		const ref operator[](size_t index) const
		{
			return m_data.at(index);
		}
		//*/

		//ref get(size_t index) { return m_data.at(index); }

		size_t GetSize() const {
			return m_data.size();
		}

		void SetSize(size_t nNewSize, size_t nGrowBy = -1) {
			m_data.resize(nNewSize);
		}
		virtual ~CSimpleArray() {}
	};
	//----------------------------------------------------------------

	typedef CSimpleArray<DBL, CIODOUBLEARRAY> CDoubleArray;
	typedef CSimpleArray<int, CIOINTARRAY> CIntArray;
	//typedef CSimpleArray<bool, CIOBOOLARRAY> CBoolArray;		//! запрещено, т.к. не возвращает ссылку на bool, bool пишется побитово!
	typedef CSimpleArray<size_t, CIOSIZE_TARRAY> CSize_tArray;
	//----------------------------------------------------------------

	//!  Шаблон CIOCustomArray
	/*!
	CIOCustomArray - аналог IO::CSimpleArray, основанный на CCustomArray
	*/
	template <class CData, TYPEID iType>
	class CCustomArray : public IO::CInterface
	{
	protected:
		::CCustomArray<size_t, CData> m_data;
	public:

		CCustomArray() {}

		virtual TYPEID GetType() const { return iType; }

		virtual bool Load(CStorage& file)
		{
			bool ret = CInterface::Load(file);

			if (!ret)
				return false;

			IO::CSize_t count;

			ret = count.Load(file);

			if (!ret)
			{
				DLOG(::CString(_T("IO::CCustomArray<")) + AllToString(iType) + _T("> : count of elements can't be loaded completely"), log_error);
				return false;
			}

			if (count() < 0)
			{
				DLOG(::CString(_T("IO::CCustomArray<")) + AllToString(iType) + _T("> : loaded count of elements < 0"), log_error);
				return false;
			}

			DLOG(::CString(_T("IO::CCustomArray<")) + AllToString(iType) + _T("> | elements count read: ") + AllToString(count()), log_info);

			m_data.resize(count());

			if (!count())
				return true;

			ret = (sizeof(CData) * count() == file.Read((void*)m_data.data(), sizeof(CData) * count()));

			if (!ret)
				DLOG(::CString(_T("IO::CCustomArray<")) + AllToString(iType) + _T("> can't be loaded completely"), log_error);
			else
				DLOG(::CString(_T("IO::CCustomArray<")) + AllToString(iType) + _T("> | bytes read: ") + AllToString(sizeof(CData) * count()), log_info);


			return ret;
		}

		virtual bool Save(CStorage& file)
		{
			bool ret = CInterface::Save(file);

			if (!ret)
				return false;

			IO::CSize_t count;
			count() = m_data.size();


			ret = count.Save(file);

			if (!ret)
				return false;

			DLOG(::CString(_T("IO::CCustomArray<")) + AllToString(iType) + _T("> | elements count saved: ") + AllToString(count()), log_info);

			if (!count())
				return true;

			try{ 
				file.Write((void*)m_data.data(), sizeof(CData) * count()); 
			}
			catch(CException* pEx) {
				DLOG(::CString(_T("IO::CCustomArray<")) + AllToString(iType) + _T("> can't be saved completely"), log_error);
				pEx->Delete();
				return false;
			}
			
			DLOG(::CString(_T("IO::CCustomArray<")) + AllToString(iType) + _T("> | bytes written: ") + AllToString(sizeof(CData) * count()), log_info);

			return true;
		}

		::CCustomArray<size_t, CData>& operator () ()
		{
			return m_data;
		}

		const ::CCustomArray<size_t, CData>& operator () () const
		{
			return m_data;
		}

		CData& operator[](ptrdiff_t index)
		{
			return m_data[index];
		}

		const CData& operator[](ptrdiff_t index) const
		{
			return m_data[index];
		}

		ptrdiff_t GetSize() const {
			return m_data.GetSize();
		}

		virtual ~CCustomArray() {}
	};
	//----------------------------------------------------------------
	typedef CCustomArray<DBL, CIOCUSTOMDOUBLEARRAY> CCustomDoubleArray;
	//----------------------------------------------------------------

	//!  Шаблон CIOBandMatrix
	/*!
	CIOBandMatrix - ленточная матрица
	*/
	template <class CData, TYPEID iType>
	class CBandMatrix : public CInterface
	{
	protected:
		::CBandMatrix<CData> m_data;
	public:

		CBandMatrix() {}

		virtual TYPEID GetType() const { return iType; }

		virtual bool Load(CStorage& file)
		{
			bool ret = CInterface::Load(file);

			if (!ret)
				return false;

			CInt band, rows;

			ret = band.Load(file);

			if (!ret)
			{
				DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T(">: band width can't be loaded completely"), log_error);
				return false;
			}

			ret = rows.Load(file);

			if (!ret)
			{
				DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T(">: rows count can't be loaded completely"), log_error);
				return false;
			}


			if ((band() < 0) || (rows() < 0))
			{
				DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T(">: loaded band width or rows count < 0"), log_error);
				return false;
			}

			DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T("> | band width read: ") + AllToString(band()), log_info);
			DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T("> | rows count read: ") + AllToString(rows()), log_info);

			m_data.resize(rows(), band());

			if (!rows || !band)
				return true;

			ret = (sizeof(CData) * m_data.size() == file.Read((void*)m_data.data(), sizeof(CData) * m_data.size()));

			if (!ret)
				DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T("> can't be loaded completely"), log_error);
			else
				DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T("> | bytes read: ") + AllToString(sizeof(CData) * m_data.size()), log_info);

			return ret;
		}

		virtual bool Save(CStorage& file)
		{
			bool ret = CInterface::Save(file);

			if (!ret)
				return false;

			CInt band, rows;
			band() = m_data.band();
			rows() = m_data.rows();

			ret = band.Save(file);

			if (!ret)
				return false;

			DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T("> | band width saved: ") + AllToString(band()), log_info);

			ret = rows.Save(file);

			if (!ret)
				return false;

			DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T("> | rows count saved: ") + AllToString(band()), log_info);


			if (!band() || !rows())
				return true;

			try{ 
				file.Write((void*)m_data.data(), sizeof(CData) * m_data.size()); 
			}
			catch(CException* pEx) {
				DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T("> can't be saved completely"), log_error);
				return false;
			}

			DLOG(::CString(_T("IO::CBandMatrix<")) + AllToString(iType) + _T("> | bytes written: ") + AllToString(sizeof(CData) * m_data.size()), log_info);

			return true;
		}

		::CBandMatrix<CData>& operator () ()
		{
			return m_data;
		}

		const ::CBandMatrix<CData>& operator () () const
		{
			return m_data;
		}

		virtual ~CBandMatrix() {}
	};
	//----------------------------------------------------------------
	typedef CBandMatrix<DBL, CIOBANDDBLMATRIX> CBandDoubleMatrix;
	//----------------------------------------------------------------

	//!  Шаблон IO::CIOArray
	/*!
	CIOArray предназначен для быстрого создания класса, хранящего, записывающего и читающего
	массив сложных данных, тип которых указывается первым параметром шаблона. По-скольку данные
	сложные, то функции Save и Load будут вызываться для каждого элемента массива.
	Имеет смысл, чтобы элементы такого массива были унаследованы (напрямую или косвенно) от
	IO::CInterface.
	В качестве второго параметра шаблона нужно указать число, которое будет записано в файл в
	качестве типа конкретного класса.
	*/
	template <class CData, TYPEID iType>
	class CIOArray : public CInterface
	{
	protected:
		std::vector<CData> m_data;
	public:

		CIOArray() {}

		GETTYPE(iType)

			virtual bool Load(CStorage& file)
		{
			bool ret = CInterface::Load(file);

			if (!ret)
				return false;

			CSize_t count;

			ret = count.Load(file);

			if (!ret)
			{
				DLOG(::CString(_T("IO::CIOArray<")) + AllToString(iType) + _T(">: count of elements can't be loaded completely"), log_error);
				return false;
			}

			if (count() < 0)
			{
				DLOG(::CString(_T("IO::CIOArray<")) + AllToString(iType) + _T(">: loaded count of elements < 0"), log_error);
				return false;
			}

			DLOG(::CString(_T("IO::CIOArray<")) + AllToString(iType) + _T("> | elements count read: ") + AllToString(count()), log_info);

			m_data.resize(count());

			if (!count())
				return true;

			for (size_t i = 0; i < count(); i++)
			{
				ret = m_data[i].Load(file);

				if (!ret)
				{
					DLOG(::CString(_T("IO::CIOArray<")) + AllToString(iType) + _T(">: element #") + AllToString(i) + _T("/") + AllToString(count()) + _T(" can't be loaded completely"), log_error);
					break;
				}
			}

			DLOG(::CString(_T("IO::CIOArray<")) + AllToString(iType) + _T("> | elements read: ") + AllToString(count()), log_info);

			return ret;
		}

		virtual bool Save(CStorage& file)
		{
			bool ret = CInterface::Save(file);

			if (!ret)
				return false;

			CSize_t count;
			count() = m_data.GetCount();

			ret = count.Save(file);

			if (!ret)
				return false;

			DLOG(::CString(_T("IO::CIOArray<")) + AllToString(iType) + _T("> | elements count written: ") + AllToString(count()), log_info);

			if (!count())
				return true;

			for (size_t i = 0; i < count(); i++)
			{
				ret = m_data[i].Save(file);

				if (!ret)
				{
					DLOG(::CString(_T("IO::CIOArray<")) + AllToString(iType) + _T(">: element #") + AllToString(i) + _T("/") + AllToString(count()) + _T(" can't be saved completely"), log_error);
					return false;
				}
			}

			DLOG(::CString(_T("IO::CIOArray<")) + AllToString(iType) + _T("> | elements written: ") + AllToString(count()), log_info);

			return true;
		}

		virtual bool Copy(CInterface* ptr)
		{

			if (!CInterface::Copy(ptr))
				return false;

			IO::CIOArray<CData, iType> *right_ptr = dynamic_cast<IO::CIOArray<CData, iType> *>(ptr);

			m_data = *right_ptr;

			return true;
		}

		operator std::vector<CData>& ()
		{
			return m_data;
		}

		std::vector<CData>& operator () ()
		{
			return m_data;
		}

		const std::vector<CData>& operator () () const
		{
			return m_data;
		}

		CData& operator[](size_t index)
		{
			return m_data[index];
		}

		const CData& operator[](size_t index) const
		{
			return m_data[index];
		}

		ptrdiff_t GetSize() const {
			return m_data.GetSize();
		}


		virtual ~CIOArray() {}
	};

	//----------------------------------------------------------------

	//!  Шаблон IO::CPtrArray
	/*!
	IO::CPtrArray - основной шаблон, если нужно автоматически сохранять/загружать целую иерархию
	классов. Предназначен для быстрого создания класса, хранящего, записывающего и читающего
	массив из указателей на различные объекты из заданной иерархии классов. Иерархия однозначно
	определяется Заголовком, т.е. классом, унаследованным от CIOHeader, почему он и указывается
	в качестве первого параметра шаблона. Элементы массива должны быть наследниками IO::CInterface.
	Этот шаблон связан с шаблоном CIOManager<CHeader>, который является фабрикой для указанной
	иерархии классов, а кроме того умеет считывать объекты иерархии из файла и создавать их
	в памяти.

	В качестве второго параметра шаблона нужно указать число, которое будет записано в файл в
	качестве типа конкретного класса-экземпляра этого шаблона.
	*/
	template <TYPEID iType>
	class CPtrArray : public CInterface
	{
	protected:
		std::vector<CInterface*> m_data;

	public:
		void DeleteObjects()
		{
			for (size_t i = 0; i < m_data.size(); i++)
			{
				delete m_data[i];
			}
		}

		CPtrArray() {}

		GETTYPE(iType)
			/*!
			Функция Load для загрузки объектов иерархии из файла:
			<ul>
			<li> 1. считывает количество сохраненных объектов;
			<li> 2. обращается к менеджеру IO::CManager и вызывает у него функцию Load для каждого сохраненного объекта.
			</ul>
			*/
			virtual bool Load(CStorage& file){
			DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | start reading"), log_info);

			bool ret = CInterface::Load(file);

			if (!ret)
				return false;

			CSize_t count;

			ret = count.Load(file);

			if (!ret)
			{
				DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> : count of elements can't be loaded completely"), log_error);
				return false;
			}

			if (count() < 0)
			{
				DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> : loaded count of elements < 0"), log_error);
				return false;
			}

			DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | elements count read: ") + AllToString(count()), log_info);

			DeleteObjects();
			m_data.clear();
			m_data.resize(count());

			if (!count())
				return true;

			LOGGER.IncShift();
			for (size_t i = 0; i < count(); i++)
			{
				DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | element [") + AllToString(i) + _T("] start reading"), log_info);

				m_data[i] = nullptr;
				m_data[i] = IOMANAGER.Load(file);

				if (!m_data[i])
				{
					DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + ::CString(_T("> | can't load object #")) + AllToString(i), log_error);
					return false;
				}
				DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | element [") + AllToString(i) + _T("] end reading"), log_info);
			}
			LOGGER.DecShift();

			DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | elements read: ") + AllToString(count()), log_info);
			DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | end reading"), log_info);

			return true;
		}

		//!  Функция Save вызывает функцию Save для каждого объекта, указатель на который хранится в массиве.
		virtual bool Save(CStorage& file)
		{
			DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | start writing"), log_info);
			bool ret = CInterface::Save(file);

			if (!ret)
				return false;

			CSize_t count;
			count() = m_data.size();

			if (!count.Save(file))
				return false;

			DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | elements count written: ") + AllToString(count()), log_info);

			if (!count())
				return true;

			LOGGER.IncShift();
			for (size_t i = 0; i < count(); i++)
			{
				DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | element [") + AllToString(i) + _T("] start writing"), log_info);
				ret = m_data[i]->Save(file);

				if (!ret)
					return false;
				DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | element [") + AllToString(i) + _T("] end writing"), log_info);
			}
			LOGGER.DecShift();

			DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | elements written: ") + AllToString(count()), log_info);
			DLOG(::CString(_T("IO::CPtrArray<")) + AllToString(iType) + _T("> | end writing"), log_info);

			return true;
		}

		virtual bool Copy(CInterface* ptr)
		{
			if (!CInterface::Copy(ptr))
				return false;

			CPtrArray<iType> *right_ptr = dynamic_cast<CPtrArray<iType> *>(ptr);

			DeleteObjects();
			m_data.clear();

			for (size_t i = 0; i < right_ptr->m_data.size(); i++)
			{
				CInterface* obj = CreateObject(right_ptr->m_data[i]->GetType());

				if (!obj)
					return false;

				m_data.push_back(obj);

				if (!obj->Copy(right_ptr->m_data[i]))
					return false;
			}
			return true;
		}

		operator std::vector<CInterface*>& ()
		{
			return m_data;
		}

		std::vector<CInterface*>& operator () ()
		{
			return m_data;
		}

		const std::vector<CInterface*>& operator () () const
		{
			return m_data;
		}

		CInterface* operator[](ptrdiff_t index)
		{
			return m_data[index];
		}

		const CInterface* operator[](ptrdiff_t index) const
		{
			return m_data[index];
		}

		size_t GetSize() const {
			return m_data.size();
		}

		virtual ~CPtrArray()
		{
			// почему я не написал удаление объектов из массива?
			// DeleteObjects();

			// при раскомментировании кладёт программу при просмотре превью-моделей
		}


	};
	//----------------------------------------------------------------

	//!  Класс IO::CString
	/*!
	IO::CString - строка (CString), которая умеет сохранять себя в файл и считывать себя из него.
	*/
	class CString : public CInterface
	{
	protected:
		CStringW m_data;

	public:

		GETTYPE(CIOSTRING)

		virtual bool Load(CStorage& file);
		virtual bool Save(CStorage& file);

		::CStringW& operator ()();
		const ::CStringW& operator ()() const;

		virtual ~CString(void) {}

	};
	//----------------------------------------------------------------

	// Оболочка для класса CMaterial, добавляющая ему свойства auto_io,
	// при этом данные, загружаемые классом CMaterial лежат в отдельном файле (текстовом),
	// а сам CIOMaterial хранится в том файле, где и объект, его использующий
	class CMaterial : public CBase
	{
		IO::CString m_filename; // имя файла, где хранятся данные для объекта CMaterial
		IO::CString m_filecontent;	//содержимое файла с материалом CMaterial

		::CMaterial m_material;
	public:
		CMaterial(void);
		CMaterial(const IO::CMaterial &pMat);

		GETTYPE(CIOMATERIAL)

		CMaterial& operator = (const CMaterial& m) {
			m_filename = m.m_filename;
			m_filecontent = m.m_filecontent;
			m_material = m.m_material;
			return *this;
		}


		virtual bool Load(CStorage& file);

		virtual bool Save(CStorage& file);

		bool ReLoadMaterial(const ::CString& filename);

		IO::CString& FileName()
		{
			return m_filename;
		}

		const IO::CString& FileName() const
		{
			return m_filename;
		}

		IO::CString& FileContent()
		{
			return m_filecontent;
		}

		const IO::CString& FileContent() const
		{
			return m_filecontent;
		}

		void SetFileContent(::CString sStr)
		{
			m_filecontent() = sStr;
		}

		::CMaterial& Material()
		{
			return m_material;
		}

		const ::CMaterial& Material() const
		{
			return m_material;
		}

		virtual ~CMaterial(void)
		{
			//if (m_material) delete m_material;
		}
	};
	//----------------------------------------------------------------

} // namespace IO