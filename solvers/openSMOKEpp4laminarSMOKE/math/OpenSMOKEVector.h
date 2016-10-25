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

#ifndef OpenSMOKE_OpenSMOKEVector_Hpp
#define OpenSMOKE_OpenSMOKEVector_Hpp

#include "OpenSMOKEStdInclude.h"
#include "OpenSMOKEBaseClass.h"
#include "OpenSMOKEUtilities.h"
#include "OpenSMOKEFunctions.h"


namespace OpenSMOKE
{

	//!  A class for vectors
	/*!
		 This is a user-friendly class to manage basic operations on vectors
	*/

	template<typename T, typename IndexPolicy = OneIndexPolicy>
	class OpenSMOKEVector : public OpenSMOKEBaseClass, IndexPolicy
	{
		template<typename T_, typename IndexPolicy_>
		friend class OpenSMOKEMatrix;
		
		public:
		
		/**
		* Default constructor (Type 1)
		*/
		OpenSMOKEVector(void);

		/**
		* Copy-initializer constructor (Type 2)
		*/
		OpenSMOKEVector(OpenSMOKEVector<T, IndexPolicy> const& rval);

		/**
		* Constructor (Type 3): the vector is sized and all the elements are set equal to 0
		*/
		OpenSMOKEVector(const int n);

		/**
		* Constructor (Type 4): the vector is sized and all the elements are provided by the user
		*/
		OpenSMOKEVector(const int n, const T v1, ...);

		/**
		* Constructor (Type 5): the vector is sized and all the elements are provided by the user
		*/
		OpenSMOKEVector(const int n, const T* values);

		/**
		* Constructor (Type 6): the vector is sized and the elements are set equal to first n elements of rhs vector
		*/
		OpenSMOKEVector(const int n, OpenSMOKEVector<T, IndexPolicy> const& rhs);

		/**
		* Constructor (Type 7): the vector is sized and the elements are set equal to the n elements of rhs vector starting from i position
		*/
		OpenSMOKEVector(const int n, const int i, OpenSMOKEVector<T, IndexPolicy> const& rhs);

		/**
		*@brief Construct the vector from file (Type 8)
		*@param fileName name of file to be read
		*@param fileFormat format of file to be read
		*/
		OpenSMOKEVector(const std::string fileName, const OpenSMOKE_File_Format fileFormat);
		
		/**
		*@brief Construct the vector from file (Type 9)
		*@param fInput stream to be read
		*@param fileFormat format of file to be read
		*/
		OpenSMOKEVector(std::istream& fInput, const OpenSMOKE_File_Format fileFormat);

		/**
		* Constructor (Type 10): all the elements are provided by the user
		*/
		OpenSMOKEVector(const std::vector<T> values);

		/**
		* Default destructor
		*/
		~OpenSMOKEVector(void);


		/**
		* Returns a pointer to internal data
		*/
		inline T* GetHandle();
		
		/**
		* Returns a pointer to internal data
		*/
		inline const T* GetHandle() const;
		

		/**
		* Returns the vector dimension (number of elements)
		*/
		inline int Size() const { return dimensions_; }

		/**
		* Returns WhoAmI
		*/
		inline unsigned int WhoAmI() const { return whoAmI_; }

		/**
		* Returns Index
		*/
		inline int Index() const { return this->index_; }

		/**
		* Returns the vector 
		*/
		inline T* Vector() const { return vector_; }


		/**
		* Returns the i element of the vector (with range control)
		*/
		inline T GetValue(const int i) const;

		/**
		* Assigns the i element of the vector (with range control)
		*/
		inline void SetValue(const int i, const T val);

		/**
		*@brief Get a const reference of a given value
		*@param i the query position
		*@return the value
		*/
		inline const T& At(const int i) const;

		/**
		* Returns the i element of the vector (with range control)
		*/
		inline T& operator () (int i);

		/**
		* Returns the i element of the vector (without range control)
		*/
		inline T& operator [] (const int i);

		/**
		* Returns the i element of the vector (without range control)
		*/
		inline const T& operator [] (const int i) const;

		/**
		* The vector is reinitialized and set equal to the first n elements of rhs vector
		*/
		void operator() (const int n, OpenSMOKEVector<T, IndexPolicy> const& rhs);

		/**
		* The vector is reinitialized and set equal to the elements of rhs vector
		*/
		void operator() (std::vector<T> const& rhs);

		/**
		* The vector is reinitialized and set equal to the n elements of rhs vector starting from start position
		*/
		void operator() (const int n, const int start, OpenSMOKEVector<T, IndexPolicy> const& rhs);
                
		/**
		* Assignment operator
		* Example: a=3 means \f$ a_i=a_i+3 \f$;
		*/
		void operator = (const T c);

		/**
		* Assignment operator
		*/
		OpenSMOKEVector<T, IndexPolicy>& operator =(OpenSMOKEVector<T, IndexPolicy> const& orig);

		/**
		* Summing another vector
		* Example: a+=b means \f$ a=a+b \f$;
		*/
		template<typename IndexPolicyRHS>
		OpenSMOKEVector<T, IndexPolicy>& operator +=(OpenSMOKEVector<T, IndexPolicyRHS> const& rval);

		/**
		* Difference between 2 vectors
		* Example: a-=b means \f$ a=a-b \f$;
		*/	
		template<typename IndexPolicyRHS>
		OpenSMOKEVector<T, IndexPolicy>& operator -=(OpenSMOKEVector<T, IndexPolicyRHS> const& rval);

		/**
		* Adds the same constant to every element
		* Example: a+=2 means \f$ a_i=a_i+2 \f$;
		*/
		inline OpenSMOKEVector<T, IndexPolicy>& operator +=(T const& rval);

		/**
		* Subtracts the same constant from every element
		* Example: a-=2 means \f$ a_i=a_i-2 \f$;
		*/	
		inline OpenSMOKEVector<T, IndexPolicy>& operator -=(T const& rval);

		/**
		* Multiplies every element by the same constant 
		* Example: a*=2 means \f$ a_i=a_i*2 \f$;
		*/
		inline OpenSMOKEVector<T, IndexPolicy>& operator *=(T const& rval);

		/**
		* Divides every element by the same constant
		* Example: a/=2 means \f$ a_i=a_i/2 \f$;
		*/	
		inline OpenSMOKEVector<T, IndexPolicy>& operator /=(T const& rval);

		/**
		* Returns the sum of the vector elements
		*/
		inline T SumElements() const;

		/**
		* Returns the sum of the absolute values of vector elements
		*/
		inline T SumAbsElements() const;


		/**
		* Remove the last n elements
		*/	
		void DeleteLastElements(const int n);

		/**
		* Print on video (verbose)
		*/
		inline void PrintOnVideo() const;
		
		/**
		*@brief Load the vector from file
		*@param fileName name of file to be read
		*@param fileFormat format of file to be read
		*/
		void Load(const std::string fileName, const OpenSMOKE_File_Format fileFormat);
		
		/**
		*@brief Load the vector from file
		*@param fInput stream to be read
		*@param fileFormat format of file to be read
		*/
		void Load(std::istream& fInput, const OpenSMOKE_File_Format fileFormat);

		/**
		*@brief Save the vector on file
		*@param fileName name of file to be written
		*@param fileFormat format of file to be written
		*/
		void Save(const std::string fileName, const OpenSMOKE_File_Format fileFormat);
		
		/**
		*@brief Save the vector on file
		*@param fOutput stream to be written
		*@param fileFormat format of file to be written
		*/
		void Save(std::ostream& fOutput, const OpenSMOKE_File_Format fileFormat);

		/**
		*@brief Returns the maximum together with the corresponding index
		*/
		inline T Max(int *imax = NULL) const;

		/**
		*@brief Returns the maximum (absolute value) together with the corresponding index
		*/
		inline T MaxAbs(int *imax = NULL) const;

		/**
		*@brief Returns the minimum together with the corresponding index
		*/
		inline T Min(int *imin = NULL) const;
		
		/**
		*@brief Returns the minimum (absolute value) together with the corresponding index
		*/
		inline T MinAbs(int *imin = NULL) const;

		/**
		*@brief Returns the minimum and the maximum together with corresponding indices
		*/
		void MinMax(int* iMin, T* min, int* iMax, T* max) const;

		/**
		*@brief Calculates the norm1
		*@return Norm 1 of the vector \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$.
		*/
		inline T Norm1() const;

		/**
		*@brief Calculates the Euclidean Norm
		*@return Eucledian norm of the vector
		*/
		inline T Norm2() const;
		
		/**
		*@brief Calculates the norm inf
		*@return Norm Inf of the vector
		*/	
		inline T NormInf() const;

		/**
		*@brief Given a vector v, appends the f value in the last position 
		*@param f The value to be appended
		*/
		void Append(const T f);
		
		/**
		*@brief Given a vector v, appends the w vector in the last position 
		*@param w The vector to be appended
		*/
		void Append(OpenSMOKEVector<T, IndexPolicy> const& w);

		/**
		*@brief Given a vector v, inserts the f value in the i position
		*@param i The position where to insert
		*@param f The value to insert
		*/
		void Insert(const int i, const T f);

		/**
		*@brief Given a sorted vector v, inserts the f value in the correct position
		*@param f The value to insert
		*/
		int InsertElementInSortedVector(const T f);

		/**
		*@brief TODO
		*@param n TODO
		*@param f The value to insert
		*/
		int InsertElementInFirstNSortedElements(const int n, const T f);

		/**
		*@brief Given a vector v, inserts the w vector in the i position
		*@param i The position where to insert
		*@param w The vector to insert
		*/
		void Insert(const int i, OpenSMOKEVector<T, IndexPolicy> const& w);

        /**
		*@brief The elements of the given vector v are copied in the array rhs
		*@param rhs The array where to copy the vector v
		*/
		void CopyTo(T* rhs) const;

		/**
		*@brief The elements of the given the array rhs are copied in v
		*@param rhs The array to copy in vector v
		*/
		void CopyFrom(const T* rhs);

		/**
		*@brief Given a vector monotonically increasing, find the index j such that an assigned f value lies between v[j] and v[j+1]
		*@param f The value to compare
		*/		
		int LocateInSortedVector(const T f);

		/**
		*@brief Given a vector monotonically increasing, find the index j such that an assigned f value lies between v[j] and v[j+1]
		        but only in the first N elements
		*@param f The value to compare
		*/		
		int LocateInFirstNSortedElements(const int position, const T f);

	// Friend functions
	public:	
	
		/**
		* Changing the dimensions of an existing vector 
		*/
		template<typename T_, typename IndexPolicy_>
		friend void ChangeDimensions(const int size, OpenSMOKEVector<T_, IndexPolicy_>* result, bool reset);

		/**
		* Swapping the elements of any 2 OpenSMOKEVector<T, IndexPolicy>. 
		* Provides an efficient method of copying. When a BzzVector remains in the scope and one leaves, they swap	
		*/	
		template<typename T_, typename IndexPolicy_>
		friend void Swap(OpenSMOKEVector<T_, IndexPolicy_> *lval, OpenSMOKEVector<T_, IndexPolicy_> *rval);

		/**
		* Add a constant to a vector
		* Example: Add(a,b,&c) means \f$ c=a+b \f$;
		*/	
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyC_>
		friend void Add(const OpenSMOKEVector<T_, IndexPolicyA_>& A, const T_ B, OpenSMOKEVector<T_, IndexPolicyC_>* C);

		/**
		* Summing 2 vectors in a third vector
		* Example: Add(a,b,&c) means \f$ c=a+b \f$;
		*/	
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_, typename IndexPolicyC_>
		friend void Add(const OpenSMOKEVector<T_, IndexPolicyA_>& A, const OpenSMOKEVector<T_, IndexPolicyB_>& B, OpenSMOKEVector<T_, IndexPolicyC_>* C);

		/**
		* Summing 2 vectors in the first vector
		* Example: Sum(&a,b) means \f$ a=a+b \f$;
		*/	
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_>
		friend void Add(OpenSMOKEVector<T_, IndexPolicyA_>* A, const OpenSMOKEVector<T_, IndexPolicyB_>& B);

		/**
		* Summing 2 vectors in the second vector
		* Example: Add(a,&b) means \f$ b=a+b \f$;
		*/		
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_>
		friend void Add(const OpenSMOKEVector<T_, IndexPolicyA_>& A, OpenSMOKEVector<T_, IndexPolicyB_>* B);

		/**
		* Summing 2 identical vectors 
		* Example: Sum(&a) means \f$ a=a+a \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend void Add(OpenSMOKEVector<T_, IndexPolicy_>* lvalRvalAndResult);

		/**
		* Difference between 2 vectors
		* Example: Difference(a,b,&c) means \f$ c=a-b \f$;
		*/	
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_, typename IndexPolicyC_>
		friend void Sub(const OpenSMOKEVector<T_, IndexPolicyA_>& lval, const OpenSMOKEVector<T_, IndexPolicyB_>& rval, OpenSMOKEVector<T_, IndexPolicyC_>* result);

		/**
		* Difference between 2 vectors
		* Example: Difference(&a,b) means \f$ a=a-b \f$;
		*/	
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_>
		friend void Sub(OpenSMOKEVector<T_, IndexPolicyA_>* lvalAndResult, const OpenSMOKEVector<T_, IndexPolicyB_>& rval);

		/**
		* Difference between 2 vectors
		* Example: Difference(b,&a) means \f$ a=b-a \f$;
		*/
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_>
		friend void Sub(const OpenSMOKEVector<T_, IndexPolicyA_>& lval, OpenSMOKEVector<T_, IndexPolicyB_>* rvalAndResult);
	
		/**
		* Difference between the same vector
		* Example: Difference(&a) means \f$ a=a-a \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend void Sub(OpenSMOKEVector<T_, IndexPolicy_>* lvalRvalAndResult);

		/**
		* Element by element product
		* Example: ElementByElementProduct(a,b,&c) means \f$ c_i=a_i*b_i \f$;
		*/
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_, typename IndexPolicyC_>
		friend void ElementByElementProduct(OpenSMOKEVector<T_, IndexPolicyA_> const& l, OpenSMOKEVector<T_, IndexPolicyB_> const& r, OpenSMOKEVector<T_, IndexPolicyC_>* s);

                /**
                * The division between the elements of a vector and a constant value
                * Example: Division(x,3.,&y) means \f$ y=x/3 \f$;
                */
                template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_>
                friend void Division(const OpenSMOKEVector<T_, IndexPolicyA_>& lval, const T_ rval, OpenSMOKEVector<T_, IndexPolicyB_>* result);

                /**
                * The division between the elements of a vector and a constant value
                * Example: Division(&x,3.) means \f$ x=x/3 \f$;
                */
                template<typename T_, typename IndexPolicy_>
                friend void Division(OpenSMOKEVector<T_, IndexPolicy_>* lvalAndResult, const T_ rval);

		/**
		* Element by element division
		* Example: ElementByElementDivision(a,b,&c) means \f$ c_i=a_i/b_i \f$;
		*/
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_, typename IndexPolicyC_>
		friend void ElementByElementDivision(OpenSMOKEVector<T_, IndexPolicyA_> const& lval, OpenSMOKEVector<T_, IndexPolicyB_> const& rval, OpenSMOKEVector<T_, IndexPolicyC_>* result);

		/**
		* The product between a constant value and a vector	
		* Example: Product(3.,x,&y) means \f$ y=3x \f$;
		*/
		template<typename T_, typename IndexPolicyA_, typename IndexPolicyB_>
		friend void Product(const T_ lval, const OpenSMOKEVector<T_, IndexPolicyA_>& rval, OpenSMOKEVector<T_, IndexPolicyB_>* result);

		/**
		* The product between a constant value and a vector	
		* Example: Product(3.,&x) means \f$ x=3x \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend void Product(const T_ lval, OpenSMOKEVector<T_, IndexPolicy_>* rvalAndResult);

		/**
		* The dot product between 2 vectors (with size check)
		* Example: DotProduct(x,y,&z) means \f$ z=todo \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend void DotProduct(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_> const& rval, T_* result);

		/**
		* The dot product between 2 vectors (with size check)
		* Example: DotProduct(x,y,&z) means \f$ z=todo \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend void UDotProduct(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_> const& rval, T_* result);

		/**
		* The dot product between 2 vectors (with size check)
		* Example: z=Dot(x,y) means \f$ z=todo \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend T_ Dot(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_> const& rval);

		/**
		* The dot product between 2 vectors (with size check)
		* Example: z=Dot(x,y) means \f$ z=todo \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend T_ UDot(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_> const& rval);

		/**
		* The dot product between 2 vectors (with size check)
		* Example: z=Dot(x,y,start,elements) means \f$ z=todo \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend void DotProduct(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_> const& rval, const int start, const int elements, T_* result);

		/**
		* The dot product between 2 vectors (with size check)
		* Example: z=UDot(x,y,start,elements) means \f$ z=todo \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend void UDotProduct(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_> const& rval, const int start, const int elements, T_* result);

		/**
		* The exp of every element of vector a is calculated and put in b
		* Example: Exp(a,&b) means \f$ b_i=exp(a_i) \f$;
		*/                
        	template<typename T_, typename IndexPolicy_>
		friend void Exp(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval);

		/**
		* The natural log of every element of vector a is calculated and put in b
		* Example: Ln(a,&b) means \f$ b_i=ln(a_i) \f$;
		*/                
        	template<typename T_, typename IndexPolicy_>
		friend void Ln(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval);

		/**
		* The denary log of every element of vector a is calculated and put in b
		* Example: Log10(a,&b) means \f$ b_i=log10(a_i) \f$;
		*/                
        	template<typename T_, typename IndexPolicy_>
		friend void Log10(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval);

		/**
		* The sin of every element of vector a is calculated and put in b
		* Example: Sin(a,&b) means \f$ b_i=sin(a_i) \f$;
		*/                
        	template<typename T_, typename IndexPolicy_>
		friend void Sin(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval);

		/**
		* The cos of every element of vector a is calculated and put in b
		* Example: Cos(a,&b) means \f$ b_i=cos(a_i) \f$;
		*/                
        	template<typename T_, typename IndexPolicy_>
		friend void Cos(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval);
		
		/**
		* The root square of every element of vector a is calculated and put in b
		* Example: Sqrt(a,&b) means \f$ b_i=(a_i)^(1/2) \f$;
		*/                
        template<typename T_, typename IndexPolicy_>
		friend void Sqrt(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval);
 
		/**
		* The square of every element of vector a is calculated and put in b
		* Example: Sqr(a,&b) means \f$ b_i=(a_i)^2 \f$;
		*/                
        template<typename T_, typename IndexPolicy_>
		friend void Sqr(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval);

		/**
		* Every element of vector a is inverted and put in b
		* Example: Inversion(a,&b) means \f$ b_i=1/a_i \f$;
		*/                
       	 template<typename T_, typename IndexPolicy_>
		friend void Inversion(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval);

		/**
		* Every element of vector a is raised to the power of c and put in b
		* Example: Pow(a,&b,c) means \f$ b_i=(a_i)^c \f$;
		*/                
        	template<typename T_, typename IndexPolicy_>
		friend void Pow(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_>* rval, const double power);

		/**
		* Every element of vector a is raised to the power of b and put in c
		* Example: Pow(a,b,&c) means \f$ c_i=(a_i)^(b_i) \f$;
		*/
		template<typename T_, typename IndexPolicy_>
		friend void Pow(OpenSMOKEVector<T_, IndexPolicy_> const& lval, OpenSMOKEVector<T_, IndexPolicy_> const& rval, OpenSMOKEVector<T_, IndexPolicy_>* result);

		/**
		* Compares two vectors
		*/
		template<typename T_, typename IndexPolicy_>
		friend bool operator == (const OpenSMOKEVector<T_, IndexPolicy_> &lval, const OpenSMOKEVector<T_, IndexPolicy_> &rval);

		/**
		* Sort a vector
		*/
		template<typename T_, typename IndexPolicy_>
		friend void Sort(OpenSMOKEVector<T_, IndexPolicy_> *result);

		/**
		* Sort a vector and keeps the pattern of ordering in vector pattern
		*/
		template<typename T_, typename IndexPolicy_>
		friend void Sort(OpenSMOKEVector<T_, IndexPolicy_> *result, OpenSMOKEVector<int, OpenSMOKE::OneIndexPolicy > *pattern);

		/**
		* Reorder a vector
		*/
		template<typename T_, typename IndexPolicy_>
		friend void Reorder(OpenSMOKEVector<T_, IndexPolicy_> *v, OpenSMOKEVector<int, OpenSMOKE::OneIndexPolicy>  &iS);

                /**
		* Reverse a vector
		*/
		template<typename T_, typename IndexPolicy_>
		friend void Reverse(OpenSMOKEVector<T_, IndexPolicy_> *result);

	protected:

		T* vector_;					/**< vector elements  */  
		int dimensions_;			/**< vector dimension */  
				
		bool shadow_;						/**< TODO */  
		bool matrixAsVector_;				/**< TODO */  
		bool subVectorAsVector_;			/**< TODO */  

	private:

		/**
		*@brief Initialize a vector (memory allocation)
		*@param size The number of vector elements
		*/
		void Initialize(const int size);	

		/**
		*@brief Copy values of anatoher vector
		*@param orig The vector to be copied
		*/
		void Copy(OpenSMOKEVector<T, IndexPolicy> const& orig);
	};

	typedef OpenSMOKEVector<unsigned int, OpenSMOKE::OneIndexPolicy >		OpenSMOKEVectorUnsignedInt;
	typedef OpenSMOKEVector<int, OpenSMOKE::OneIndexPolicy >				OpenSMOKEVectorInt;
	typedef OpenSMOKEVector<float, OpenSMOKE::OneIndexPolicy >				OpenSMOKEVectorFloat;
	typedef OpenSMOKEVector<double, OpenSMOKE::OneIndexPolicy >				OpenSMOKEVectorDouble;
	typedef OpenSMOKEVector<std::string, OpenSMOKE::OneIndexPolicy >		OpenSMOKEVectorString;
	typedef OpenSMOKEVector<bool, OpenSMOKE::OneIndexPolicy >				OpenSMOKEVectorBool;
}

#include "OpenSMOKEVector.hpp"

#endif	// OpenSMOKE_OpenSMOKEVector_Hpp

