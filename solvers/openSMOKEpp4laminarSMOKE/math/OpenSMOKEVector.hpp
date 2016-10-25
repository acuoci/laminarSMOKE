/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

//#include "math/OpenSMOKEVector.h"
#include "math/OpenSMOKEMatrix.h"
#include <typeinfo>
#include <stdio.h>
#include <string.h>

#if OPENSMOKE_USE_MKL == 1
	#include "mkl.h"
#elif OPENSMOKE_USE_OPENBLAS == 1
	#include "lapacke.h"
#endif

namespace OpenSMOKE
{

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(void)
	{
		SetCounters();
		Initialize(0);
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(OpenSMOKEVector<T, IndexPolicy> const& rval)
	{
		SetCounters();
		Initialize(rval.dimensions_);
		Copy(rval);
	/*
		if(rval.shadow_ == false)
		{
			Initialize(rval.dimensions_);
			if(dimensions_ != 0)
				memcpy(vector_, rval.vector_, (dimensions_+this->index_)*sizeof(T));
		}
		else
		{
			Initialize(0);
			OpenSMOKEVector<T, IndexPolicy> aux; 
			aux = rval;
			Swap(this, &aux);
			count_--;
			whoAmI_ = aux.whoAmI_;
		}
		*/
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(const int n)
	{
		SetCounters();
		Initialize(n);
		if(dimensions_ != 0)
			memset(vector_+this->index_, 0, dimensions_*sizeof(T));
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(const int n, const T v1, ...)
	{
		SetCounters();
		Initialize(n);
		T *w = vector_ + this->index_;
		va_list pointerList;
		va_start(pointerList, v1);
		*w = v1;
		for(int i=1+this->index_;i<n+this->index_;i++)
			*++w = va_arg(pointerList, T);
		va_end(pointerList);
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(const int n, const T* values)
	{
		SetCounters();
		Initialize(n);
		if(dimensions_ != 0)
		{
			T* w = vector_;
			w+=this->index_;
			memcpy(w, values, n*sizeof(T));
		}
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(const int n, OpenSMOKEVector<T, IndexPolicy> const& rhs)
	{
		SetCounters();
		Initialize(n);
	
		if(n > rhs.dimensions_)
			ErrorMessage("OpenSMOKEVector(const int n, OpenSMOKEVector<T, IndexPolicy> const& rhs) requires n<=rhs.dimensions");

		if(dimensions_ != 0)
		{
			T* w = vector_;
			T* v = rhs.vector_;
			w+=this->index_;
			v+=rhs.index_;
			memcpy(w, v, n*sizeof(T));
		}
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(const int n, const int i, OpenSMOKEVector<T, IndexPolicy> const& rhs)
	{
		SetCounters();
		Initialize(n);

		if(n > rhs.dimensions_-i+1)
			ErrorMessage("OpenSMOKEVector(const int n, const int i, OpenSMOKEVector<T, IndexPolicy> const& rhs) requires n<=rhs.dimensions-i+1");
		
		if(dimensions_ != 0)
		{
			T* w = vector_;
			T* v = rhs.vector_;
			w+=this->index_;
			v += i;
			memcpy(w, v, n*sizeof(T));
		}
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
	{
		SetCounters();
		Load(fileName, fileFormat);
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(std::istream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		SetCounters();
		Load(fInput, fileFormat);
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::OpenSMOKEVector(const std::vector<T> values)
	{
		SetCounters();
		Initialize(values.size());
		if(dimensions_ != 0)
		{
			T* w = this->vector_;
			w+=this->index_;
			for (typename std::vector<T>::const_iterator it = values.begin(); it!=values.end(); ++it)
			{
				*w = *it;
				w++;
			}
		}
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>::~OpenSMOKEVector(void)
	{
		if( matrixAsVector_ == true || subVectorAsVector_ == true )
		{
			countInScope_--;
			vector_ = 0;
			dimensions_ = 0;
			matrixAsVector_ = false;
			subVectorAsVector_ = false;
			return;
		}

		if(dimensions_ != 0)
			delete[] vector_;

		vector_ = 0;
		dimensions_ = 0;
		countInScope_--;
	}

	template<typename T, typename IndexPolicy>
	inline T* OpenSMOKEVector<T, IndexPolicy>::GetHandle()
	{
		return vector_ + this->index_;
	}
	
	template<typename T, typename IndexPolicy>
	inline const T* OpenSMOKEVector<T, IndexPolicy>::GetHandle() const	
	{
		return vector_ + this->index_;
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Copy(OpenSMOKEVector<T, IndexPolicy> const& orig)
	{
		if (dimensions_!=orig.dimensions_)
			ChangeDimensions(orig.Size(), this, false);
			//ErrorMessage("OpenSMOKEVector<T, IndexPolicy>::Copy(OpenSMOKEVector<T, IndexPolicy> const& orig) Dimension check failure");
						
		if(dimensions_ != 0)
			memcpy(vector_+this->index_, orig.vector_+orig.index_, dimensions_*sizeof(T));
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Initialize(const int size)
	{
		matrixAsVector_		= false;
		subVectorAsVector_	= false;
		shadow_				= false;
		
		if (size < 0)
		{
			ErrorMessage("Size must be a positive integer");
		}
		else if(size == 0)
		{
			dimensions_ = 0;
			//vector_		= 0;
		}
		else
		{
			dimensions_ = size;
			vector_ = new T[size+this->index_]; 
			if(!vector_)
				ErrorMessage("Size must be a positive integer");
		}
	}

	template<typename T, typename IndexPolicy>
	inline const T& OpenSMOKEVector<T, IndexPolicy>::At(const int i) const
	{
		return vector_[i];
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::GetValue(const int i) const
	{
		if( (i<IndexPolicy::index_) || (i>dimensions_-1+this->index_) )
			ErrorMessage("Vector index outside the ranges");
		return vector_[i];
	}
	
	template<typename T, typename IndexPolicy>
	inline void OpenSMOKEVector<T, IndexPolicy>::SetValue(const int i, const T value)
	{
		if( (i<IndexPolicy::index_) || (i>dimensions_-1+this->index_) )
			ErrorMessage("Vector index outside the ranges");
		vector_[i] = value;
	}

	template<typename T, typename IndexPolicy>
	inline T& OpenSMOKEVector<T, IndexPolicy>::operator() (int i)
	{
		if( (i<this->index_) || (i>dimensions_-1+this->index_) )
			ErrorMessage("Vector index outside the ranges");
		return vector_[i];
	}

	template<typename T, typename IndexPolicy>
	inline T& OpenSMOKEVector<T, IndexPolicy>::operator[] (const int i)
	{
		return vector_[i];
	}

	template<typename T, typename IndexPolicy>
	inline const T& OpenSMOKEVector<T, IndexPolicy>::operator[] (const int i) const
	{
		return vector_[i];
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::operator() (const int n, OpenSMOKEVector<T, IndexPolicy> const& rhs)
	{
		Initialize(n);

		if(n>rhs.dimensions_)
			ErrorMessage("operator(const int n, OpenSMOKEVector<T, IndexPolicy> const& rhs) requires n<=rhs.dimensions");
		
		if(dimensions_ != 0)
		{
			T* w = vector_;
			T* v = rhs.vector_;
			w+=this->index_;
			v+=rhs.index_;
			memcpy(w, v, n*sizeof(T));
		}
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::operator() (std::vector<T> const& rhs)
    {
		Initialize(rhs.size());
		std::cout << rhs.size() << std::endl;
		if(dimensions_ != 0)
		{
			T* w = this->vector_;
			w+=this->index_;
			for (typename std::vector<T>::const_iterator it = rhs.begin(); it!=rhs.end(); ++it)
			{
				std::cout << *w << " " << *it << std::endl;
				*w= *it;
				std::cout << *w << " " << *it << std::endl;
				w++;
			}
		}
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::operator() (const int n, const int start, OpenSMOKEVector<T, IndexPolicy> const& rhs)
	{
		Initialize(n);

		int max_n = rhs.dimensions_ - start + this->index_;

		if( n>max_n )
			ErrorMessage("operator(const int n, const int start, OpenSMOKEVector<T, IndexPolicy> const& rhs) requires n>max_n");
		else if ( max_n<1 )
			ErrorMessage("operator(const int n, const int start, OpenSMOKEVector<T, IndexPolicy> const& rhs) requires max_n>0");
		else if(start<this->index_)
			ErrorMessage("operator(const int n, const int start, OpenSMOKEVector<T, IndexPolicy> const& rhs) requires start>=index");
		else if (start>rhs.dimensions_+this->index_-1)
			ErrorMessage("operator(const int n, const int start, OpenSMOKEVector<T, IndexPolicy> const& rhs) requires start<=rhs.dimensions");

		if(dimensions_ != 0)
		{
			T* w = vector_;
			T* v = rhs.vector_ + start;
			w+=this->index_;
			memcpy(w, v, n*sizeof(T));
		}
	}
        
        template<typename T, typename IndexPolicy>
        void OpenSMOKEVector<T, IndexPolicy>::operator = (T value)
	{
                if(value == 0 && dimensions_ > 0)
                        memset(vector_+this->index_,0,dimensions_*sizeof(T));
                else
		{
                        for(int i=this->index_;i<dimensions_+this->index_;i++)
                                vector_[i] = value;
		}
	}

	template<typename T, typename IndexPolicy>
	OpenSMOKEVector<T, IndexPolicy>& OpenSMOKEVector<T, IndexPolicy>::operator =(OpenSMOKEVector<T, IndexPolicy> const& orig)
	{
			if (&orig!=this)
				Copy(orig);

			return *this;
/*
		int who = whoAmI_;
		if(dimensions_ != orig.dimensions_)
		{
			delete vector_;
			Initialize(orig.dimensions_);
			count_--;
			countInScope_--;
		}

		whoAmI_ = who;
		if(dimensions_ > 0)
			memcpy(vector_+this->index_, orig.vector_+orig.index_, dimensions_*sizeof(T));
		
		return *this;
		*/
	}

	template<typename T, typename IndexPolicy>
	template<typename IndexPolicyRHS>
	OpenSMOKEVector<T, IndexPolicy>& OpenSMOKEVector<T, IndexPolicy>::operator += (OpenSMOKEVector<T, IndexPolicyRHS> const& rval)
	{
		if(dimensions_ != rval.Size())
			ErrorMessage("operator+=(const OpenSMOKEVector<T, IndexPolicy>& rval) dimension check failure");
		if(whoAmI_ == rval.WhoAmI())
			Add(this);
		else
			Sum(dimensions_, vector_+this->index_, rval.Vector()+rval.Index());
		return *this;
	}

	template<typename T, typename IndexPolicy>
	template<typename IndexPolicyRHS>
	OpenSMOKEVector<T, IndexPolicy>& OpenSMOKEVector<T, IndexPolicy>::operator-=(OpenSMOKEVector<T, IndexPolicyRHS> const& rval)
	{
		if(dimensions_ != rval.Size())
			ErrorMessage("operator-=(const OpenSMOKEVector<T, IndexPolicy>& rval) Dimension check failure");
		
		if(whoAmI_ == rval.WhoAmI())
			Sub(this);
		else
			Difference(dimensions_, vector_+this->index_, rval.Vector()+rval.Index());
		return *this;
	}

	
	template<typename T, typename IndexPolicy>
	inline OpenSMOKEVector<T, IndexPolicy>& OpenSMOKEVector<T, IndexPolicy>::operator +=(T const& rval)
	{
		Sum(dimensions_ , rval, vector_+this->index_); 
		return *this;
	}

	template<typename T, typename IndexPolicy>
	inline OpenSMOKEVector<T, IndexPolicy>& OpenSMOKEVector<T, IndexPolicy>::operator -=(T const& rval)
	{
		Sum(dimensions_ , -rval, vector_+this->index_); 
		return *this;
	}
	
	template<typename T, typename IndexPolicy>
	inline OpenSMOKEVector<T, IndexPolicy>& OpenSMOKEVector<T, IndexPolicy>::operator *=(T const& rval)
	{
		Prod(dimensions_ , rval, vector_+this->index_); 
		return *this;
	}

	template<typename T, typename IndexPolicy>
	inline OpenSMOKEVector<T, IndexPolicy>& OpenSMOKEVector<T, IndexPolicy>::operator /=(T const& rval)
	{
		if (rval==0)
			ErrorMessage("operator/=(T const& rval) Division by zero");

		Prod(dimensions_ , 1./rval, vector_+this->index_); 
		return *this;
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::SumElements() const
	{
		T sum = 0;
		T* x = vector_+this->index_;
		for(int i=this->index_;i<dimensions_+this->index_;i++)
			sum += *x++;
		return sum;
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::SumAbsElements() const
	{
		T sum = 0;
		T* x = vector_+this->index_;
		for(int i=this->index_;i<dimensions_+this->index_;i++)
			sum += abs(*x++);
		return sum;
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::DeleteLastElements(const int n) 
	{
		if (n<=0)
			return;
		
		if(dimensions_ <= n)
		{
			OpenSMOKEVector<T, IndexPolicy> v;
			Swap(&v,this);
			return;
		}

		OpenSMOKEVector<T, IndexPolicy> v(dimensions_ - n);
		memmove(v.vector_+this->index_, vector_+this->index_, (dimensions_-n)*sizeof(T));
		Swap(&v,this);		
	}

	template<typename T, typename IndexPolicy>
	inline void OpenSMOKEVector<T, IndexPolicy>::PrintOnVideo() const
	{
		std::cout << typeid(OpenSMOKEVector<T, IndexPolicy>).name() << std::endl;
		std::cout << "Counters:      " << whoAmI_ << "/" << countInScope_ << "/" << count_ << std::endl;
		std::cout << "Size:          " << dimensions_ << std::endl;
		for(int i=this->index_;i<dimensions_+this->index_;i++)
			std::cout << i << " " << vector_[i] << std::endl;
		std::cout << std::endl;
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Load(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			std::ifstream fInput(fileName.c_str(), std::ios::in);
			if (fInput.fail())
				ErrorMessage("The " + fileName + " file does not exist");

			Load(fInput, fileFormat);
			
			fInput.close();
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			std::ifstream fInput(fileName.c_str(), std::ios::in | std::ios::binary);
			if (fInput.fail())
				ErrorMessage("The " + fileName + " binary file does not exist");

			Load(fInput, fileFormat);

			fInput.close();
		}
	}
		
	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Load(std::istream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			int n;
			fInput >> n;
                        Initialize(n);
			for(int i=this->index_;i<dimensions_+this->index_;i++)
				fInput >> vector_[i];
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			int n;

			// Reading vector size
			if(!fInput.read(reinterpret_cast<char *>(&n), sizeof(int)))
				ErrorMessage("I was unable to read from binary file (1)");
			
			// Initializing
			Initialize(n);

			// Reading vector elements
			for(int i=this->index_;i<dimensions_+this->index_;i++)
			if(!fInput.read(reinterpret_cast<char *>(&vector_[i]), sizeof(T)))
				ErrorMessage("I was unable to read from binary file (2)");
			else
				std::cout << n << " " << vector_[i] << std::endl;
		}
	}

		
	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Save(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			std::ofstream fOutput(fileName.c_str(), std::ios::out);
			if (fOutput.fail())
				ErrorMessage("The " + fileName + " file cannot be open");

			Save(fOutput, fileFormat);
			
			fOutput.close();
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			std::ofstream fOutput(fileName.c_str(), std::ios::out | std::ios::binary);
			if (fOutput.fail())
				ErrorMessage("The " + fileName + " binary file cannot be open");

			Save(fOutput, fileFormat);

			fOutput.close();
		}
	}
		
	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Save(std::ostream& fOutput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fOutput << dimensions_ << std::endl;
			for(int i=this->index_;i<dimensions_+this->index_;i++)
				fOutput << std::setprecision(16) << std::scientific << vector_[i] << std::endl;
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			// Writing vector size
			if(!fOutput.write( reinterpret_cast<char *>(&dimensions_), sizeof(int)))
				ErrorMessage("I was unable to write on binary file");

			// Writing vector elements
			for(int i=this->index_;i<dimensions_+this->index_;i++)
			if(!fOutput.write( reinterpret_cast<char *>(&vector_[i]), sizeof(T)))
				ErrorMessage("I was unable to write on binary file");
		}
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::Max(int *imax) const
	{
		if(dimensions_ == 0)
			return 0;
		T xmax = OpenSMOKE::Max(dimensions_, vector_ + this->index_, imax);
		if(imax != 0)
			(*imax)+=this->index_;
		return xmax;
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::MaxAbs(int *imax) const
	{
		if(dimensions_ == 0)
			return 0;
		T xmax = OpenSMOKE::MaxAbs(dimensions_, vector_ + this->index_, imax);
		if(imax != 0) 
			(*imax)+=this->index_;
		return xmax;
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::Min(int *imin) const
	{
		if(dimensions_ == 0)
			return 0;
		T xmin = OpenSMOKE::Min(dimensions_, vector_ + this->index_, imin);
		if(imin != 0) 
			(*imin)+=this->index_;
		return xmin;
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::MinAbs(int *imin) const
	{
		if(dimensions_ == 0)
			return 0;
		T xmin = OpenSMOKE::MinAbs(dimensions_, vector_ + this->index_, imin);
		if(imin != 0) 
			(*imin)+=this->index_;
		return xmin;
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::MinMax(int* iMin, T* min, int* iMax, T* max) const
	{
		if(dimensions_ == 0)
		{
			*iMin = 0;
			*iMax = 0;
			*min = 0;
			*max = 0;
			return;
		}

		*iMin = *iMax = this->index_ ;
		*min = *max = vector_[this->index_];
	
		for(int i=1+this->index_;i<dimensions_+this->index_;i++)
		{
			if(vector_[i] < *min)
			{
				*min = vector_[i];
				*iMin = i;
			}
			else if(vector_[i] > *max)
			{
				*max = vector_[i];
				*iMax = i;
			}
		}
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::Norm1() const
	{
		T norm = 0.;
		T* x = vector_+this->index_;
		for(int i=this->index_;i<dimensions_+this->index_;i++)
			norm += abs(*x++);
		return norm;
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::Norm2() const
	{
		T* x = vector_+this->index_;
		return SqrtSumSqr(dimensions_, x);
	}

	template<typename T, typename IndexPolicy>
	inline T OpenSMOKEVector<T, IndexPolicy>::NormInf() const
	{
		return MaxAbs();
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Append(const T f)
	{
		OpenSMOKEVector<T, IndexPolicy> v(dimensions_+1);
		if(dimensions_ != 0)
			memcpy(v.vector_, vector_,(dimensions_+this->index_)*sizeof(T));
		v[dimensions_+this->index_] = f;
		Swap(&v,this);
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Append(OpenSMOKEVector<T, IndexPolicy> const& w)
	{
		int sz = w.dimensions_;
		OpenSMOKEVector<T, IndexPolicy> v(dimensions_ + sz);
		if(dimensions_ != 0)
			memcpy(v.vector_, vector_, (dimensions_+this->index_)*sizeof(T));
	
		memcpy(v.vector_ + dimensions_ + this->index_, w.vector_ + w.index_, sz*sizeof(T));
		Swap(&v,this);
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Insert(const int i, const T f)
	{
		if(i < this->index_) 
			ErrorMessage("Insert(const int i, const T f) requires i>=index");
		
		if(i < dimensions_+this->index_)
		{
			OpenSMOKEVector<T, IndexPolicy> v(dimensions_+1);
			memmove(v.vector_, vector_, i*sizeof(T));
			v[i] = f;
			memmove(v.vector_+i+1, vector_+i, (dimensions_-i+1)*sizeof(T));
			Swap(&v,this);		
		}
		else
		{
			OpenSMOKEVector<T, IndexPolicy> v(i+1-this->index_);
			memmove(v.vector_ + this->index_, vector_ + this->index_, dimensions_ * sizeof(T));
			v[i] = f;
			Swap(&v,this);		
		}
	}

	template<typename T, typename IndexPolicy>
	int OpenSMOKEVector<T, IndexPolicy>::InsertElementInSortedVector(const T f)
	{
		const int j = LocateInSortedVector(f);
		Insert(j+this->index_, f);
		return j+this->index_;
	}

	template<typename T, typename IndexPolicy>
	int OpenSMOKEVector<T, IndexPolicy>::InsertElementInFirstNSortedElements(const int n, const T f)
	{
		if(n >= this->dimensions_)
		{
			OpenSMOKEVector<T, IndexPolicy> v(this->dimensions_+10);
			if(this->dimensions_ != 0)
				memcpy(v.vector_, this->vector_, (this->dimensions_+1)*sizeof(T));
			Swap(&v,this);
		}
	
		const int j = LocateInFirstNSortedElements(n,f);
		if(j != n)
		{
			memmove(this->vector_+j+1+this->index_, this->vector_+j+this->index_, (this->dimensions_-j-1)*sizeof(T));
		}
		(*this)[j+this->index_] = f;
		return j+this->index_;
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::Insert(const int i, OpenSMOKEVector<T, IndexPolicy> const& w)
	{
		int sz = w.dimensions_;
		
		OpenSMOKEVector<T, IndexPolicy> v(dimensions_ + sz);
		if(i < this->index_) 
			ErrorMessage("Insert(const int i, OpenSMOKEVector<T, IndexPolicy> const& w) requires i>=index");

		if(i < dimensions_+this->index_)
		{
			memmove(v.vector_, vector_,i*sizeof(T));
			memmove(v.vector_ + i, w.vector_ + w.index_, sz*sizeof(T));
			memmove(v.vector_ + sz + i, vector_ + i, (dimensions_ - i + this->index_)*sizeof(T));
			Swap(&v,this);		
		}
		else
			Append(w);
	}

	template<typename T, typename IndexPolicy>
	int OpenSMOKEVector<T, IndexPolicy>::LocateInSortedVector(const T f)
	{
		int n = this->dimensions_;

		if(n < this->index_)
			return this->index_-1;
		if(n == this->index_)
		{
			if(this->vector_[this->index_] < f)	return this->index_;
			else                                    return this->index_-1;
		}
		
		int upper,lower,jm;
		lower = 0;
		upper = n + 1;
		if(this->vector_[n-1+this->index_] > this->vector_[this->index_])
		{
			while(upper - lower > 1)
			{
				jm = (upper + lower) >> 1;
				if(f > this->vector_[jm-1+this->index_])
					lower = jm;
				else
					upper = jm;
			}
			return lower-1+this->index_;
		}
		else
		{
			while(upper - lower > 1)
			{
				jm = (upper + lower) >> 1;
				if(f < this->vector_[jm-1+this->index_])
					lower = jm;
				else
					upper = jm;
			}
			return lower-1+this->index_;
		}
	}

	template<typename T, typename IndexPolicy>
	int OpenSMOKEVector<T, IndexPolicy>::LocateInFirstNSortedElements(const int position, const T f)
	{
		int n = position;
		
		if(n > this->dimensions_)
			n = this->dimensions_;
	
		if(n < 1)
			return 0;
	
		if(n == 1)
		{
			if(this->vector_[1-1+this->index_] < f)	return 1;
			else									return 0;
		}
	
		int upper,lower;
		lower = 0;
		upper = n+1;
	
		if(this->vector_[n-1+this->index_] > this->vector_[1-1+this->index_])
		{
			while(upper - lower > 1)
			{
				const int jm = (upper + lower) >> 1;
				if(f > this->vector_[jm-1+this->index_])
					lower = jm;
				else
					upper = jm;
			}
			return lower;
		}
		else
		{
			while(upper - lower > 1)
			{
				const int jm = (upper + lower) >> 1;
				if(f < this->vector_[jm-1+this->index_])
					lower = jm;
				else
					upper = jm;
			}
			return lower;
		}
	}

	// Friend Functions
	template<typename T, typename IndexPolicy>
	void Swap(OpenSMOKEVector<T, IndexPolicy>* lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		Swap(&lval->vector_,     &rval->vector_);
		Swap(&lval->dimensions_, &rval->dimensions_);
	}

	template<typename T, typename IndexPolicy>
	void ChangeDimensions(const int size, OpenSMOKEVector<T, IndexPolicy>* result, bool reset)
	{
		if(size != result->dimensions_)
		{
			if(result->matrixAsVector_ == true)
				ErrorMessage("ChangeDimensions", "You can not change the dimensions when you are using a matrix as a vector");
	
			if(result->subVectorAsVector_ == true)
				ErrorMessage("ChangeDimensions", "You can not change the dimensions when you are using a subvector as a vector");
		
                        if (result->dimensions_ != 0)
                                delete[] result->vector_;

			result->Initialize(size);
		}
	
		// Resetting
		if(reset == true && result->dimensions_ != 0)
			memset(result->vector_+result->index_, 0, (result->dimensions_)*sizeof(T) );
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB, typename IndexPolicyC>
	void Add(const OpenSMOKEVector<T, IndexPolicyA>& lval, const OpenSMOKEVector<T, IndexPolicyB>& rval, OpenSMOKEVector<T, IndexPolicyC>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Sum(const OpenSMOKEVector<T, IndexPolicy>& lval, const OpenSMOKEVector<T, IndexPolicy>& rval, OpenSMOKEVector<T, IndexPolicy>* result)",  
							"Vector dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"Sum(const OpenSMOKEVector<T, IndexPolicy>& lval, const OpenSMOKEVector<T, IndexPolicy>& rval, OpenSMOKEVector<T, IndexPolicy>* result)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
            const T* ptrval = rval.GetHandle();
			      T* ptresult = result->GetHandle();

            vdAdd( n, ptlval, ptrval, ptresult );
		}
		#else
		{
			if(result->whoAmI_ == lval.whoAmI_)
				(*result) += rval;
			else if(result->whoAmI_ == rval.whoAmI_)
				Add(lval,result);
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				Sum(result->dimensions_, lval.vector_ + lval.index_, rval.vector_ + rval.index_, result->vector_ + result->index_);
			}
		}
		#endif
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyC>
	void Add(const OpenSMOKEVector<T, IndexPolicyA>& lval, const T rval, OpenSMOKEVector<T, IndexPolicyC>* result)
	{	
		if(result->whoAmI_ == lval.whoAmI_)
			(*result) += rval;
		else
		{
			ChangeDimensions(lval.dimensions_, result, false);
			Sum(result->dimensions_, lval.vector_ + lval.index_, rval, result->vector_ + result->index_);
		}
	}


	template<typename T, typename IndexPolicyA, typename IndexPolicyB>
	void Add(OpenSMOKEVector<T, IndexPolicyA>* lvalAndResult, const OpenSMOKEVector<T, IndexPolicyB>& rval)
	{
		(*lvalAndResult) += rval;
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB>
	void Add(const OpenSMOKEVector<T, IndexPolicyA>& lval, OpenSMOKEVector<T, IndexPolicyB>* rvalAndResult)
	{
		(*rvalAndResult) += lval;
	}

	template<typename T, typename IndexPolicy>
	void Add(OpenSMOKEVector<T, IndexPolicy>* lvalRvalAndResult)
	{
		Sum(lvalRvalAndResult->dimensions_, lvalRvalAndResult->vector_+lvalRvalAndResult->index_);
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB, typename IndexPolicyC>
	void Sub(const OpenSMOKEVector<T, IndexPolicyA>& lval, const OpenSMOKEVector<T, IndexPolicyB>& rval, OpenSMOKEVector<T, IndexPolicyC>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Difference(const OpenSMOKEVector<T, IndexPolicy>& lval, const OpenSMOKEVector<T, IndexPolicy>& rval, OpenSMOKEVector<T, IndexPolicy>* result)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"Difference(const OpenSMOKEVector<T, IndexPolicy>& lval, const OpenSMOKEVector<T, IndexPolicy>& rval, OpenSMOKEVector<T, IndexPolicy>* result)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
            const T* ptrval = rval.GetHandle();
			      T* ptresult = result->GetHandle();

            vdSub( n, ptlval, ptrval, ptresult );
		}
		#else
		{
			if(result->whoAmI_ == lval.whoAmI_)
				(*result) -= rval;
			else if(result->whoAmI_ == rval.whoAmI_)
				Sub(lval,result);
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				Difference(result->dimensions_, lval.vector_+lval.index_, rval.vector_ + rval.index_, result->vector_ + result->index_); 
			}
		}
		#endif

	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB>
	void Sub(OpenSMOKEVector<T, IndexPolicyA>* lvalAndResult, const OpenSMOKEVector<T, IndexPolicyB>& rval)
	{
		(*lvalAndResult) -= rval;
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB>
	void Sub(const OpenSMOKEVector<T, IndexPolicyA>& lval, OpenSMOKEVector<T, IndexPolicyB>* rvalAndResult)
	{
		if(lval.dimensions_ != rvalAndResult->dimensions_)
			ErrorMessage("Difference(const OpenSMOKEVector<T, IndexPolicy>& lval, OpenSMOKEVector<T, IndexPolicy>* rvalAndResult)", "Dimension check failure");
		DifferenceBis(lval.dimensions_, lval.vector_+lval.index_, rvalAndResult->vector_+rvalAndResult->index_);
	}

	template<typename T, typename IndexPolicy>
	void Sub(OpenSMOKEVector<T, IndexPolicy>* lvalRvalAndResult)
	{
		T *w = lvalRvalAndResult->vector_+lvalRvalAndResult->index_;
		for(int i=1;i<=lvalRvalAndResult->dimensions_;i++)
			*w++ = 0;
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB>
	void Product(const T lval, const OpenSMOKEVector<T, IndexPolicyA>& rval, OpenSMOKEVector<T, IndexPolicyB>* result)
	{
    
	#if OPENSMOKE_USE_MKL == 1
    	{
        	if(rval.dimensions_ != result->dimensions_)
            	ChangeDimensions(rval.dimensions_, result, false);

        	result->Copy(rval); 
        	int one = 1;
        	int n = rval.dimensions_;
        	T* ptresult = result->GetHandle();
        	dscal(&n, &lval, ptresult, &one);
    	}
	#else
    	if (rval.WhoAmI() == result->WhoAmI())
        	Prod(result->dimensions_, lval, result->vector_ + result->Index());
    	else
    	{
        	ChangeDimensions(rval.dimensions_, result, false);
        	Prod(rval.dimensions_, lval, rval.vector_ + rval.Index(), result->vector_ + result->Index());
    	}
	#endif
	}

	template<typename T, typename IndexPolicy>
	void Product(const T lval, OpenSMOKEVector<T, IndexPolicy>* rvalAndResult)
	{
	#if OPENSMOKE_USE_MKL == 1
    	{
        	int one = 1;
        	int n = rvalAndResult->dimensions_;
        	T* ptr_rvalAndResult = rvalAndResult->GetHandle();

        	dscal(&n, &lval, ptr_rvalAndResult, &one);
    	}
	#else
    	Prod(rvalAndResult->dimensions_, lval, rvalAndResult->vector_+rvalAndResult->Index()); 
	#endif
	}

	template<typename T, typename IndexPolicy>
	void DotProduct(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, T* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage("DotProduct(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, T* result)", "Dimension check failure");
		
		*result = Dot(lval.dimensions_, lval.vector_+lval.index_, rval.vector_+rval.index_); 
	}

	template<typename T, typename IndexPolicy>
	void UDotProduct(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, T* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage("UDotProduct(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, T* result)", "Dimension check failure");
		
		*result = UDot(lval.dimensions_, lval.vector_+lval.index_, rval.vector_+rval.index_); 
	}

	template<typename T, typename IndexPolicy>
	T Dot(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval)
	{
		T result = 0;
		DotProduct(lval, rval, &result);
		return result;
	}

	template<typename T, typename IndexPolicy>
	T UDot(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval)
	{
		T result = 0;
		UDotProduct(lval, rval, &result);
		return result;
	}

	template<typename T, typename IndexPolicy>
	void DotProduct(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, const int start, const int elements, T* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage("Dot(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, const int start, const int elements, T* result)", "Dimension check failure");

		*result = Dot(elements, lval.vector_+start, rval.vector_+start);
	}

	template<typename T, typename IndexPolicy>
	void UDotProduct(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, const int start, const int elements, T* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage("UDot(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, const int start, const int elements, T* result)", "Dimension check failure");

		*result = UDot(elements, lval.vector_+start, rval.vector_+start);
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB, typename IndexPolicyC>
	void ElementByElementProduct(OpenSMOKEVector<T, IndexPolicyA> const& lval, OpenSMOKEVector<T, IndexPolicyB> const& rval, OpenSMOKEVector<T, IndexPolicyC>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"ElementByElementProduct(OpenSMOKEVector<T, IndexPolicy> const& l, OpenSMOKEVector<T, IndexPolicy> const& r, OpenSMOKEVector<T, IndexPolicy>* s)",
							"Dimension check failure");
		
        #if OPENSMOKE_USE_MKL == 1
        {
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"ElementByElementProduct(OpenSMOKEVector<T, IndexPolicy> const& l, OpenSMOKEVector<T, IndexPolicy> const& r, OpenSMOKEVector<T, IndexPolicy>* s)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
			const T* ptrval = rval.GetHandle();
				  T* ptresult = result->GetHandle();

			vdMul( n, ptlval, ptrval, ptresult );      
        }
        #else
        {
			if(rval.whoAmI_ == result->whoAmI_ || lval.whoAmI_ == result->whoAmI_)
			{
				OpenSMOKEVector<T, IndexPolicyC> aux;
				ElementByElementProduct(lval, rval, &aux);	
				Swap(&aux,result);						
			}
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				T* w = lval.vector_ + lval.index_;
				T* v = rval.vector_ + rval.index_;
				T* x = result->vector_ + result->index_;
				for(int i=1;i<=lval.dimensions_;i++)
					(*x++) = (*w++) * (*v++);
			}        
        }
        #endif
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB>
        void Division(const OpenSMOKEVector<T, IndexPolicyA>& lval, const T rval, OpenSMOKEVector<T, IndexPolicyB>* result)
	{
	#if OPENSMOKE_USE_MKL == 1
    	{
        	if(lval.dimensions_ != result->dimensions_)
            	ChangeDimensions(lval.dimensions_, result, false);
        
        	result->Copy(lval);
        	double rval_inv = 1. / double(rval);
        	int one = 1;
        	int n = lval.dimensions_;
        	T* ptresult = result->GetHandle();
        	dscal(&n, &rval_inv, ptresult, &one);        
    	}
	#else
    	if (lval.WhoAmI() == result->WhoAmI())
        	Div(result->dimensions_, result->vector_ + result->Index(), rval);
    	else
    	{
        	ChangeDimensions(lval.dimensions_, result, false);
        	Div(lval.dimensions_, lval.vector_ + lval.Index(), rval, result->vector_ + result->Index());
    	}
	#endif
	}

	template<typename T, typename IndexPolicy>
	void Division(OpenSMOKEVector<T, IndexPolicy>* lvalAndResult, const T rval)
	{
	#if OPENSMOKE_USE_MKL == 1
    	{
        	double rval_inv = 1. / double(rval);
        	int one = 1;
        	int n = lvalAndResult->dimensions_;
        	T* ptr_lvalAndResult = lvalAndResult->GetHandle();
        	dscal(&n, &rval_inv, ptr_lvalAndResult, &one);        
    	}
	#else
		Div(lvalAndResult->dimensions_, lvalAndResult->vector_ + lvalAndResult->Index(), rval);
	#endif
	}

	template<typename T, typename IndexPolicyA, typename IndexPolicyB, typename IndexPolicyC>
	void ElementByElementDivision(OpenSMOKEVector<T, IndexPolicyA> const& lval, OpenSMOKEVector<T, IndexPolicyB> const& rval, OpenSMOKEVector<T, IndexPolicyC>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"ElementByElementDivision(OpenSMOKEVector<T, IndexPolicy> const& l, OpenSMOKEVector<T, IndexPolicy> const& r, OpenSMOKEVector<T, IndexPolicy>* s)",
							"Dimension check failure");
		
        #if OPENSMOKE_USE_MKL == 1
        {
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
			const T* ptrval = rval.GetHandle();
				  T* ptresult = result->GetHandle();

			vdDiv( n, ptlval, ptrval, ptresult );      
        }
        #else
        {
			if(rval.whoAmI_ == result->whoAmI_ || lval.whoAmI_ == result->whoAmI_)
			{
				OpenSMOKEVector<T, IndexPolicyC> aux;
				ElementByElementDivision(lval, rval, &aux);	
				Swap(&aux,result);						
			}
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				T* w = lval.vector_ + lval.index_;
				T* v = rval.vector_ + rval.index_;
				T* x = result->vector_ + result->index_;
				for(int i=1;i<=lval.dimensions_;i++)
					(*x++) = (*w++) / (*v++);
			}     
        }
        #endif
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::CopyTo(T* rhs) const
	{		
		if(dimensions_ != 0)
			memcpy(rhs, vector_+this->index_, dimensions_*sizeof(T));
	}

	template<typename T, typename IndexPolicy>
	void OpenSMOKEVector<T, IndexPolicy>::CopyFrom(const T* rhs)
	{		
		if(dimensions_ != 0)
			memcpy(vector_+this->index_, rhs, dimensions_*sizeof(T));
	}
        
	template<typename T, typename IndexPolicy>
	void Exp(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Exp(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdExp( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = std::exp(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T, typename IndexPolicy>
	void Ln(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Ln(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdLn( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = log(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T, typename IndexPolicy>
	void Log10(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Log10(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdLog10( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = log10(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T, typename IndexPolicy>
	void Sin(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Sin(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdSin( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = std::sin(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T, typename IndexPolicy>
	void Cos(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Cos(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdCos( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = std::cos(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T, typename IndexPolicy>
	void Sqrt(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Sqrt(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdSqrt( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = sqrt(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T, typename IndexPolicy>
	void Sqr(const OpenSMOKEVector<T, IndexPolicy>& lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage(	"Sqr(const OpenSMOKEVector<T, IndexPolicy>& lval, OpenSMOKEVector<T, IndexPolicy>* rval)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
                  T* ptrval = rval->GetHandle();

            vdSqr( n, ptlval, ptrval);
		}
		#else
		{
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = lval.vector_[i]*lval.vector_[i]; 
		}
		#endif
	}

	template<typename T, typename IndexPolicy>
	void Inversion(const OpenSMOKEVector<T, IndexPolicy>& lval, OpenSMOKEVector<T, IndexPolicy>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage(	"Inversion(const OpenSMOKEVector<T, IndexPolicy>& lval, OpenSMOKEVector<T, IndexPolicy>* rval)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
                  T* ptrval = rval->GetHandle();

            vdInv( n, ptlval, ptrval);
		}
		#else
		{
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = 1./lval.vector_[i]; 
		}
		#endif
	}

	template<typename T, typename IndexPolicy>
	void Pow(const OpenSMOKEVector<T, IndexPolicy>& lval, OpenSMOKEVector<T, IndexPolicy>* rval, const double power)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage(	"Pow(const OpenSMOKEVector<T, IndexPolicy>& lval, OpenSMOKEVector<T, IndexPolicy>* rval, const double power)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
                  T* ptrval = rval->GetHandle();

            vdPowx( n, ptlval, power, ptrval);
		}
		#else
		{
                for(int i=rval->index_;i<rval->dimensions_+rval->index_;i++)
                    rval->vector_[i] = std::pow(lval.vector_[i], power);
		}
		#endif
	}

	template<typename T, typename IndexPolicy>
	void Pow(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, OpenSMOKEVector<T, IndexPolicy>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Pow(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, OpenSMOKEVector<T, IndexPolicy>* result)",
							"Dimension check failure");
		
		if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"Pow(OpenSMOKEVector<T, IndexPolicy> const& lval, OpenSMOKEVector<T, IndexPolicy> const& rval, OpenSMOKEVector<T, IndexPolicy>* result)",
								"Vector dimension check failure");

        #if OPENSMOKE_USE_MKL == 1
        {
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
			const T* ptrval = rval.GetHandle();
				  T* ptresult = result->GetHandle();

			vdPow( n, ptlval, ptrval, ptresult );      
        }
        #else
        {
                for(int i=rval.index_;i<rval.dimensions_+rval.index_;i++)
                    result->vector_[i] = std::pow(lval.vector_[i],rval.vector_[i]);     
        }
        #endif
	}

	template<typename T, typename IndexPolicy>
	bool operator == (const OpenSMOKEVector<T, IndexPolicy> &lval, const OpenSMOKEVector<T, IndexPolicy> &rval)
	{
		if(lval.whoAmI_ == rval.whoAmI_) 
			return true;
		
		bool flag = true;
		if(lval.dimensions_ != rval.dimensions_) 
			flag = false;
		else
		{
			if(memcmp(lval.vector_ + lval.index_,rval.vector_ + rval.index_, rval.dimensions_*sizeof(T)) == 0)
				flag = true;
			else flag = false;
		}
		return flag;
	}

	template<typename T, typename IndexPolicy>
	void Sort(OpenSMOKEVector<T, IndexPolicy> *result)
	{
		Sort(result->dimensions_, result->vector_+result->index_); 
	}

	template<typename T, typename IndexPolicy>
	void Sort(OpenSMOKEVector<T, IndexPolicy> *result, OpenSMOKEVector<int, OpenSMOKE::OneIndexPolicy >  *iS)
	{
		ChangeDimensions(result->dimensions_, iS, true);
		for(int i=1;i<=result->dimensions_;i++)
			(*iS)[i] = i;
		Sort(result->dimensions_, result->vector_ + result->index_,iS->vector_ + iS->index_);
	}

	template<typename T, typename IndexPolicy>
	void Reorder(OpenSMOKEVector<T, IndexPolicy> *v, OpenSMOKEVector<int, OpenSMOKE::OneIndexPolicy>  &iS)
	{
		int size =v->Size();

		if(size != iS.Size())
			ErrorMessage("void Reorder(OpenSMOKEVector<T, IndexPolicy> *v, OpenSMOKEVector<int, OpenSMOKE::OneIndexPolicy>  &iS)", "iS has a wrong size");

		T *w = new T[size+1];
		if(!w)
			ErrorMessage("void Reorder(OpenSMOKEVector<T, IndexPolicy> *v, OpenSMOKEVector<int, OpenSMOKE::OneIndexPolicy>  &iS)", "no enough memory");

		for(int i=1;i<=size;i++)
			w[i] = (*v)[iS[i]];
		
		Swap(&w, &v->vector_);

		delete[] w;
	}
        
        template<typename T, typename IndexPolicy>
	void Reverse(OpenSMOKEVector<T, IndexPolicy> *result)
	{
            const int n = result->Size();
            const int m = n/2;	
            if(n < 2)
                return;
	
            T *v = result->vector_ + result->index_;
            T *w = result->vector_ + n + result->index_ -1;
	
            for(unsigned int i=1;i<=m;i++,v++,w--)
                Swap(v,w);
	}
}
