#pragma once
#include "CustomArray.h"
#include "AllTypes.h"
#include "./Logger/Logger.h"

//---------------------------------------------------------------------------

#define MATRIX CMatrix<size_t, DBL>

template <typename SizeT, typename Cell>
class CMatrix : public CCustomArray<SizeT, Cell>
{
protected:

	size_t rows_; // rows number
	size_t cols_; // cols number

public:

	typedef CCustomArray<SizeT, Cell> Array;
	typedef Cell Cell;
	typedef SizeT SizeT;

	CMatrix()
	{
		rows_ = cols_ = 0;
	}

	CMatrix(SizeT rows, SizeT cols) : Array(rows * cols),
		rows_(rows), cols_(cols)
	{}

	CMatrix(const CMatrix& matrix) : CCustomArray(matrix),
		rows_(matrix.rows_), cols_(matrix.cols_)
	{}

	bool resize(SizeT rows = 1, SizeT cols = 1, Cell value = Cell(0))
	{
		if (rows && cols) {
			if (Array::resize(rows * cols))
			{
				rows_ = rows;
				cols_ = cols;
				for (size_t i = 0; i < size_; i++)
					data_[i] = value;
				return true;
			}
			return false;

		}
		else {
			free();
		}
		return true;
	}

	void reset(const Cell& c = Cell(0))
	{
		for (SizeT i = 0; i < size_; i++)
			data_[i] = c;
	}

	void free()
	{
		Array::free();
		rows_ = cols_ = 0;
	}

	void swap(CMatrix& matrix)
	{
		Array::swap(matrix);

		rows_ = matrix.rows_;
		cols_ = matrix.cols_;
	}

	size_t rows() const
	{
		return rows_;
	}

	size_t cols() const
	{
		return cols_;
	}

	void symmetrize()
	{
		for (SizeT i = 0; i < rows_; i++)
			for (SizeT j = 0; j < cols_; j++)
				cell(j, i) = cell(i, j);
	}

	Cell& cell(SizeT i, SizeT j)
	{
		return cell2d(i, j, cols_);
	}

	const Cell& cell(SizeT i, SizeT j) const
	{
		return cell2d(i, j, cols_);
	}

	Cell& bandCell(SizeT i, SizeT j)
	{
		return cell(i, (i + j) % cols_);
	}

	CMatrix& operator= (const CMatrix& matrix)
	{
		if (data_ != matrix.data_) {
			Array::operator=(matrix);
			rows_ = rows;
			cols_ = cols;
		}
		return *this;
	}

	void WriteToLog(bool val) // val = true - выводить значение ячейки, val = false - выводить *, если не 0, иначе - 0
	{
		DLOG(CString(_T("Matrix log start ->")), log_info);
		for (SizeT i = 0; i < rows_; i++)
		{
			CString tmp = _T("");
			for (SizeT j = 0; j < cols_; j++)
			{
				if (val)
					tmp += AllToString(cell(i, j)) + CString(_T(","));
				else
					if (cell(i, j)) tmp += _T("*");
					else tmp += _T("0");
			}

			DLOG(CString(_T("[")) + AllToString(i) + CString(_T("]:")) + tmp, log_info);
		}	
		DLOG(CString(_T("Matrix log end   <-")), log_info);
	}

};


// стандартная ленточная матрица
#define BANDMATRIX CBandMatrix<DBL>

// symmetric band matrix class
// stores only upper-half of band (with diagonal, of course)
template <typename Cell>
class CBandMatrix : public CCustomArray<size_t, Cell>
{
	static Cell zero_; // остальные ячейки

	protected:

	size_t rows_; // кол-во строк
	size_t band_; // размер полосы

	public:

	typedef CCustomArray<size_t, Cell> Array;
	typedef Cell Cell;

	CBandMatrix()
	{
		rows_ = band_ = 0;
	}

	CBandMatrix(size_t rows, size_t band) : Array(rows * band),
													   rows_(rows),
													   band_(band)
	{}

	CBandMatrix(const CBandMatrix& matrix) : CCustomArray(matrix),
										     rows_(matrix.rows_),
											 band_(matrix.band_)
	{}

	bool resize(size_t rows = 0, size_t band = 0, Cell value = Cell(0))
	{
		if (rows && band) {
			if (Array::resize(rows * band))
			{
				rows_ = rows;
				band_ = band;
				for (size_t i = 0; i < size_; i++)
					data_[i] = value;
				return true;
			}
			return false;

		}
		else {
			free();
		}
		return true;
	}

	void reset(const Cell& c = Cell(0))
	{
		for (SizeT i = 0; i < size_; i++)
			data_[i] = c;
	}

	//обнуление строк и полосы
	void free()	
	{
		Array::free();
		rows_ = band_ = 0;
	}

	void swap(CBandMatrix& matrix)
	{
		Array::swap(matrix);

		rows_ = matrix.rows_;
		band_ = matrix.band_;
	}

	//возвращает кол-во строк
	size_t rows() const
	{
		return rows_;
	}

	//!возвращает размер полосы (или размер матрицы?)
	size_t band() const
	{
		return band_;
	}

	//возвращает кол-во строк (???)
	size_t cols() const
	{
		return rows();
	}

	///! прямое индексирование
	Cell& direct_cell(size_t i, size_t j)
	{
		return cell2d(i, j, band_);
	}
		
	///! прямое индексирование
	const Cell& direct_cell(size_t i, size_t j) const
	{
		return cell2d(i, j, band_);
	}

	///! универсальное ("безопасное") индексирование
	Cell& cell(size_t i, size_t j)
	{
		if ((j >= i) && (j - i < band_))
			return cell2d(i, j - i, band_);

		if ((j < i) && (i - j < band_))
			return cell2d(j, i - j, band_);

		return zero_;
	}

	CBandMatrix& operator= (const CBandMatrix& matrix)
	{
		if (data_ != matrix.data_) {
			Array::operator=(matrix);
			rows_ = rows;
			band_ = band;
		}
		return *this;
	}

	//записывает ленточную (НЕ полную) матрицу
	// true - выводить значение ячейки, false - выводить *, если не 0, иначе - 0
	void WriteToLog(bool val){

		DLOG(CString(_T("Matrix log start ->")), log_info);
		for (SizeT i = 0; i < rows(); i++){
			CString tmp = _T("");
			//for (SizeT j = 0; j < cols(); j++)
			for (SizeT j = i; j < min(cols(), i+band_); j++){

				if (val) {
					tmp += AllToString(cell(i, j)) + CString(_T(","));
				}else{
					if (cell(i, j)) {
						tmp += AllToString((int)(min(1, fabs(cell(i, j)) + 1))); //_T("*");
					}else{
						tmp += _T("0");
					}
				}
			}

			CString count = _T("0");
			if (i>9){
				count.Empty();	//если больше 9 строк, то перестаём вставлять 0: "01, 02,...09, _10"
			}

			DLOG(CString(_T("[")) + count + AllToString(i) + CString(_T("]:")) + tmp, log_info);
		}	
		DLOG(CString(_T("Matrix log end   <-")), log_info);
	}

	//записывает квадратную (полную) матрицу ПРИ ОТЛАДКЕ
	// true - выводить значение ячейки, false - выводить *, если не 0, иначе - 0
	void WriteToLogFullMatrix(bool val) {

		LOG(CString(_T("Matrix log start ->")), log_info);
		size_t row_size = rows();
		size_t col_size = cols();

		for (SizeT i = 0; i < row_size; i++) {
			CString tmp = _T("");
			for (SizeT j = 0; j < col_size; j++) {

				if (val) {
					tmp += AllToString(cell(i, j)) + CString(_T(","));
				}
				else {
					if (cell(i, j)) {
						tmp += AllToString((int)(min(1, fabs(cell(i, j)) + 1))); //_T("*");
					}
					else {
						tmp += _T("0");
					}
				}
			}

			CString count = _T("0");
			if (i>9) {
				count.Empty();	//если больше 9 строк, то перестаём вставлять 0: "01, 02,...09, _10"
			}

			LOG(CString(_T("[")) + count + AllToString(i) + CString(_T("]:")) + tmp, log_info);
		}
		LOG(CString(_T("Matrix log end   <-")), log_info);
	}

	//записывает квадратную (полную) матрицу ПРИ ОТЛАДКЕ
	// true - выводить значение ячейки, false - выводить *, если не 0, иначе - 0
	void WriteToLogFullMatrixOnDebug(bool val) {

		DLOG(CString(_T("Matrix log start ->")), log_info);
		size_t row_size = rows();
		size_t col_size = cols();
		
		for (SizeT i = 0; i < row_size; i++) {
			CString tmp = _T("");
			for (SizeT j = 0; j < col_size; j++) {

				if (val) {
					tmp += AllToString(cell(i, j)) + CString(_T(","));
				}
				else {
					if (cell(i, j)) {
						tmp += AllToString((int)(min(1, fabs(cell(i, j)) + 1))); //_T("*");
					}
					else {
						tmp += _T("0");
					}
				}
			}

			CString count = _T("0");
			if (i>9) {
				count.Empty();	//если больше 9 строк, то перестаём вставлять 0: "01, 02,...09, _10"
			}

			DLOG(CString(_T("[")) + count + AllToString(i) + CString(_T("]:")) + tmp, log_info);
		}
		DLOG(CString(_T("Matrix log end   <-")), log_info);
	}
/*
	void WriteToLog()
	{
		DLOG(CString(_T("BandMatrix log start ->")), log_info);
		for (int i = 0; i < rows(); i++)
		{
			CString tmp = _T("");
			for (int j = 0; j < rows(); j++)
				tmp += AllToString(cell(i, j)) + _T(",");

			DLOG(CString(_T("[")) + AllToString(i) + CString(_T("]:")) + tmp, log_info);
		}	
		DLOG(CString(_T("BandMatrix log end   <-")), log_info);
	}
*/
};

template <typename Cell>
Cell CBandMatrix<Cell>::zero_ = Cell(0);