#pragma once

#include "memory.h"
//---------------------------------------------------------------------------

// Шаблон класса "Массив элементов"
// Cell - тип элементов массива (как правило, double),
// SizeT - тип переменной, отвечающей за размер массива (как правило, size_t).

template <typename SizeT, typename Cell>
class CCustomArray
{
	protected:

	Cell* data_;
	SizeT size_;

	public:

	typedef Cell Cell;
	typedef SizeT SizeT;
	
	
	CCustomArray()
	{
		data_ = nullptr;
		size_ = 0;
	}
	
	CCustomArray(SizeT size, Cell* data = nullptr)
	{
		if (data_ && size) {
			data_ = nullptr;
			resize(size);
		}
	}
	

	CCustomArray(const CCustomArray& array)
	{
		if (data_ != array.data_) {
			data_ = nullptr;
			resize(array.size_);
			memcpy(data_, array.data_, size_ * sizeof(Cell));
		}
	}

	bool resize(SizeT newSize)
	{
		if (newSize == size_) {
			return true;
		}

		free();

		if (newSize) {
			data_ = new Cell[newSize];
			//if (data_) {
			size_ = newSize;
			return true;
			//}
			//return false;
		}

		return false;
	}

	void free()
	{
		//if (data_) 
		delete [] data_;
		//data_ = nullptr;
		size_ = 0;
	}

	void swap(CCustomArray& array)
	{
		Cell* tempData = data_;
		SizeT tempSize = size_;

		data_ = array.data_;
		size_ = array.size_;

		array.data_ = tempData;
		array.size_ = tempSize;
	}

	SizeT size() const
	{
		return size_;
	}

	Cell* data()
	{
		return data_;
	}

	const Cell* data() const
	{
		return data_;
	}

	Cell& operator[] (SizeT i)
	{
		return data_[i];
	}

	const Cell& operator[] (SizeT i) const
	{
		return data_[i];
	}

	CCustomArray& operator= (const CCustomArray& array)
	{
		if (data_ != array.data_) {
			data_ = nullptr;
			resize(array.size_);
			memcpy(data_, array.data_, size_ * sizeof(Cell));
		}
		return *this;
	}

	// индексирование 2мерной матрицы, когда заранее известно
	// конкретное число столбцов
	template <SizeT Cols>
	Cell& cell2d(SizeT i, SizeT j)
	{
		return data_[i * Cols + j];
	}

	template <SizeT Cols>
	const Cell& cell2d(SizeT i, SizeT j) const
	{
		return data_[i * Cols + j];
	}

	// индексирование 2мерной матрицы
	Cell& cell2d(SizeT i, SizeT j, SizeT cols)
	{
		return data_[i * cols + j];
	}

	const Cell& cell2d(SizeT i, SizeT j, SizeT cols) const
	{
		return data_[i * cols + j];
	}

	~CCustomArray()
	{
		free();
	}
};