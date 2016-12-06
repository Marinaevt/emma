#pragma once
#include "AllTypes.h"
#include "IOInterface.h"
#include "2DPoint.h"
#include "CustomArray.h"
#include "BandMatrix.h"
#include "./Mesher/Mesh.h"

#include <vector>
#include <list>
#include <unordered_set>

#include <gl/GL.h>
#include "GLFont.h"
#include "Geometry/Utils.h"

int BitCount(int src, UNINT stop = 32); // возвращает количество установленных битов в src начиная с нулевого бита и заканчивая битом под номером stop
int Index2BitIndex(int src, int index); // возвращает номер index-го установленного бита (-1, если его нет)
int BitIndex2Number(int index); // возвращает число с установленным битом index
bool HasBit(int src, int bits); // возвращает true, если (src & bits)
int GetBitIndex(int src, int bits); // возвращает порядковый номер первого установленного бита из bits в src

void WriteBitFuncsToLog();	// --- отладка битовых функций ---

// ------- Points -------

typedef std::vector<Math::C2DPoint> C2DPoints;
typedef IO::CSimpleArray<Math::C2DPoint, CIO2DPOINTS> CIO2DPoints;
typedef IO::CSimpleArray<Math::C2DPoint, CIO2DPOINTS> CIO2DNodes;

// ------- Elements -------
class C2DElement3
{
public:
	union // индексы узлов
	{
		ptrdiff_t n[3];
		struct { ptrdiff_t n0, n1, n2; };
		struct { ptrdiff_t _0, _1, _2; };

		//struct { int first, second, third; };
	};

	union // индексы соседних (прилегающих) элементов
	{
		ptrdiff_t adj[3];
		struct { ptrdiff_t a0, a1, a2; };
	};

public:
	//
	int HaveEdge(ptrdiff_t node1, ptrdiff_t node2) const;

	//! Содержит ли элемент этот узел?
	bool IsHaveNode(ptrdiff_t nNode);
};

class C2DElement4
{
public:
	union
	{
		int n[4];
		struct { int _0, _1, _2, _3; };
		struct { int n0, n1, n2, n3; };
		struct { int first, second, third, fourth; };
	};

	union // индексы соседних (прилегающих) элементов
	{
		int adj[4];
		struct { int a0, a1, a2, a3; };
	};
};


class CEdge
{
public:
	union // индексы узлов
	{
		ptrdiff_t n[2];
		struct { ptrdiff_t n0, n1; };
		struct { ptrdiff_t _0, _1; };
	};

public:
	CEdge (ptrdiff_t node0, ptrdiff_t node1)
	{
		n0 = node0;
		n1 = node1;
	}

	// == , !=
	bool operator == (const CEdge& right)
	{
		return ((n0 == right.n0 && n1 == right.n1) || (n0 == right.n1 && n1 == right.n0));
	}
	bool operator != (const CEdge& right)
	{
		return !this->operator==(right);
	}

	// == , != CONST
	bool operator == (const CEdge& right) const { 
		return ((n0 == right.n0 && n1 == right.n1) || (n0 == right.n1 && n1 == right.n0));
	}
	bool operator != (const CEdge& right) const { 
		return !this->operator==(right); 
	}

	//Hash for unordered_set
	static size_t GetHash(const CEdge& p) { 
		//побитово сравниваем XOR т.к. отрезок может быть направлен в другую сторону (например не (3,4), а (4,3))
		return FemLibrary::Utils::hash_combine( p.n0, p.n1) ^ FemLibrary::Utils::hash_combine( p.n1, p.n0);
	}
};

enum eFields // физические сущности
{
	cm_x =  static_cast<size_t>(0x00001), // координата x центра
	cm_y =  static_cast<size_t>(0x00002), // координата y центра

	def   = static_cast<size_t>(0x00004), // деформация
	def_s = static_cast<size_t>(0x00008), // скорость деформаций

	str   = static_cast<size_t>(0x00010), // напряжение
	//str_s = (int)0x00020, // скорость напряжений
	avg_d = static_cast<size_t>(0x00020), // средние деформации

	eps =   static_cast<size_t>(0x00040), // модуль Юнга
	nju =   static_cast<size_t>(0x00080), //

	mju =   static_cast<size_t>(0x00100), // коэффициент вязкости

	sf_a =  static_cast<size_t>(0x00200), // коэффициент A функции формы
	sf_b =  static_cast<size_t>(0x00400), // коэффициент B функции формы
	sf_c =  static_cast<size_t>(0x00800), // коэффициент C функции формы
	sqr  =  static_cast<size_t>(0x01000), // площадь элемента

	int_d  = static_cast<size_t>(0x02000), // интенсивность деформаций
	int_ds = static_cast<size_t>(0x04000), // интенсивность скорости деформаций
	int_s  = static_cast<size_t>(0x08000), // интенсивность напряжений
	int_ss = static_cast<size_t>(0x10000), // резерв

	hyd_p  = static_cast<size_t>(0x20000), // накопленное гидростатическое давление
	smo_d  = static_cast<size_t>(0x40000), // сглаженные средние деформации	

	avg_ds  = static_cast<size_t>(0x80000), // средние скорости деформации
	smo_ds  = static_cast<size_t>(0x100000), // сглаженные скорости деформации

	sigma_y = static_cast<size_t>(0x200000), // напряжение по y
	sigma_x = static_cast<size_t>(0x400000), // напряжение по x
	sigma_fi = static_cast<size_t>(0x800000), // напряжение при осессимметрии (добавочное)

	//cm_x | cm_y | def | def_s | str | avg_d | eps | nju | mju | sf_a | sf_b | sf_c | sqr | int_d | int_ds | int_s | int_ss | hyd_p | smo_d | avg_ds // все поля
	all_f = static_cast<size_t>(-1)	//для x32 достаточно 0xffffffff, для x64 - нет.

	// и др.
};

enum eElementTypes // типы элементов
{
	tri =  static_cast<int>(0x01), // трехвершинные
	quad = static_cast<int>(0x02)  // четырехвершинные
};

enum e2DMeshTypes
{
	mesh3  = 1,
	mesh4  = 2,
	mesh34 = 3,
	border3 = 4,
	border4 = 5,
	border34 = 6,
	border43 = 6
};

template <int N>
class CFields
{
public:
	DBL f[N];

	DBL& operator() (int i)
	{
		return f[i];
	}
};

class C2DMeshInterface : public IO::CClass
{
protected:
	IO::CUnInt m_nfields;
	IO::CUnInt m_efields;

public:
	CIO2DPoints m_nodes;
	
	IO::CCustomDoubleArray m_nodefields;
	
	IO::CIntArray m_bordernodes,
				m_borderelems; //,
				//m_bordernodes_marks;	//номера узлов с засечками

public: 
	C2DMeshInterface(){
		RegisterMembers();
	}

	//Регистрируем члены класса для сохранения и загрузки
	void RegisterMembers();

	/// выделяет память под узлы и поля узлов
	bool InitNodes(size_t nodes, int nfields);

	/// выделяет память под поля узлов, исходя из количества узлов
	bool InitNodeFields(int nfields);

	int GetNodeFields()
	{
		return m_nfields();
	}

	virtual bool InitElements(size_t elements, int efields) = 0;

	int GetElementFields()
	{
		return m_efields();
	}

	virtual int GetMeshType() = 0;
	virtual int GetElementTypes() = 0;

	/* IO::CBase part */

	virtual TYPEID GetType() const = 0;	

	/* Нахождение граничного узла по номеру */
	Math::C2DPoint GetBorderNode(size_t nBorderNode) const;

	//Функция возвращает вектор, направленный по нормале к поверхности в граничной точке nBoundaryNode
	Math::C2DPoint GetNormalS(ptrdiff_t nBoundaryNode, int nSMode) const;

	//Функция, возвращает площадь кольца вокруг граничного узла nBoundaryNode
	DBL GetCircleSquare(ptrdiff_t nBoundaryNode) const;

	size_t GetBorderSize() const {return m_bordernodes.GetSize();}

	void AddBorderNode(int node);
	void AddBorderElem(int elem);

	DBL& GetNField(size_t node, int field);
	DBL GetNField(size_t node, int field) const;

	virtual void FillSTM(BANDMATRIX& m, std::vector<DBL>& rp, const CMaterial& mat, DBL dt, size_t startE, size_t endE) = 0;

	virtual Math::C2DRect GetBoundingBox() const;
	virtual bool GetBoundingBox(Math::C2DRect &rect);

	virtual void DrawGL(GLParam &parameter) = 0;
	virtual void DrawGL3D(GLParam &parameter) = 0;

	virtual std::vector< std::vector<size_t> > GetNodeLinks() const = 0;
	virtual void ReorderNodes(const std::vector<size_t>& neworder);

	virtual ~C2DMeshInterface() {}
};


class C2DMesh3 : public C2DMeshInterface
{
public:
	IO::CSimpleArray<C2DElement3, C2DEL3ARRAY> m_elements;
	IO::CCustomDoubleArray m_elemfields;

	GETTYPE (CIO2DMESH3)

	C2DMesh3(){
		RegisterMembers();
	}

	C2DMesh3(const CMesh &mesh);

	void RegisterMembers();

	/// выделяет память под элементы и поля элементов
	virtual bool InitElements(size_t elements, int efields);

	/// выделяет память под поля элементов, исходя из количества элементов
	bool InitElemFields(int efields);

	virtual int GetMeshType()
	{
		return mesh3;
	}

	virtual int GetElementTypes()
	{
		return tri;
	}

	void SortNodesInElements(bool ccw = true);
	void FillBorderNodes();	//Заполнение граничных узлов
	void GenerateMarks();	//Создание горизонтальных рисок (использовалось)

	virtual void FillSTM(BANDMATRIX& m, std::vector<DBL>& rp, const CMaterial& mat, DBL dt, size_t startE, size_t endE);
	void CalcFormCoeffs();

	DBL& GetEField(ptrdiff_t elem, int field);

	virtual void DrawGL(GLParam &parameter);
	virtual void DrawGL3D(GLParam &parameter);

	virtual std::vector< std::vector<size_t> > GetNodeLinks() const;

	virtual void ReorderNodes(const std::vector<size_t>& neworder);
	void WriteToLog() const;
};

IODEF (C2DMesh3)
