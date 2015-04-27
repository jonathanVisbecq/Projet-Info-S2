/*************************************************************************************
	*
	* Aliases and some overload to provide convenient manipulation of basic linear algebra
	* operations.
	*
	*************************************************************************************/

#ifndef MY_ARRAY_HPP
#define MY_ARRAY_HPP

#include "stl_headers.hpp"

/*
	* Typedefs to emulate linear algebra vector/matrix calculus
	*/
template<size_t dim>
using Array = std::array<double,dim>;

template<size_t dim1,size_t dim2>
using Matrix = std::array<std::array<double,dim2>,dim1>;


/*
	* Add two Arrays element-wise
	*/
template<size_t dim>
Array<dim> operator+(const Array<dim>& a1, const Array<dim>& a2)
{
				Array<dim> a_tp(a1);
				for(int i=0; i<dim; ++i)
								a_tp.at(i) += a2.at(i);

				return a_tp;
}

/*
	* Substract two Arrays element-wise
	*/
template<size_t dim>
Array<dim> operator-(const Array<dim>& a1, const Array<dim>& a2)
{
				Array<dim> a_tp(a1);
				for(int i=0; i<dim; ++i)
								a_tp.at(i) -= a2.at(i);

				return a_tp;
}

/*
	* Multiply two Arrays element-wise
	*/
template<size_t dim>
Array<dim> operator*(const Array<dim>& a1, const Array<dim>& a2)
{
				Array<dim> a_tp(a1);
				for(int i=0; i<dim; ++i)
								a_tp.at(i) *= a2.at(i);

				return a_tp;
}

/*
	* Multiply an array by a scalar
	*/
template<size_t dim>
Array<dim> operator*(const double& s, const Array<dim>& a)
{
				Array<dim> a_tp(a);
				for(int i=0; i<dim; ++i)
								a_tp.at(i) *= s;

				return a_tp;
}

/*
	* Divide an array by a scalar
	*/
template<size_t dim>
Array<dim> operator/(const Array<dim>& a, const double& s)
{
				Array<dim> a_tp(a);
				for(int i=0; i<dim; ++i)
								a_tp.at(i) /= s;

				return a_tp;
}

/*
	* Multiplication of a Matrix and an Array
	*/
template<size_t dim1, size_t dim2>
Array<dim1> operator*(const Matrix<dim1,dim2>& m, const Array<dim2>& a)
{
				Array<dim1> a_tp;
				for(int i=0; i<dim1; ++i)
								for(int j=0; j<dim2; ++j)
												a_tp.at(i) += m.at(i).at(j) * a.at(j);

				return a_tp;
}


/*
	* Scalar product of two Arrays
	*/
template<size_t dim>
inline
double scalar_prod(const Array<dim>& a1, const Array<dim>& a2)
{
				return std::inner_product(a1.begin(),a1.end(),a2.begin(),0.);
}


/*
	* Return an Array<dim> with all elements equals to 'val'
	*/
template<size_t dim>
Array<dim> makeFill(double val)
{
				Array<dim> a;
				for(int i=0; i<dim; ++i)
								a.at(i) = val;

				return a;
}

/*
	* Return a Matrix<dim1,dim2> with all elements equals to 'val'
	*/
template<size_t dim1, size_t dim2>
Matrix<dim1,dim2> makeFill(double val)
{
				Matrix<dim1,dim2> m;
				for(int i=0; i<dim1; ++i)
								for(int j=0; j<dim2; ++j)
												m.at(i).at(j) = val;

				return m;
}

/*
	* Overload to serialize Array<dim>
	*/
template<size_t dim>
std::ostream& operator<<(std::ostream& stream, const Array<dim>& a)
{
				for(auto it=a.begin(); it!=a.end(); ++it)
								stream << *it << std::endl;

				return stream;
}

/*
	* Return indices that sort an array in descending order
	*/
template<size_t dim>
std::array<unsigned,dim> sort_idx_desc(const Array<dim>& a)
{
				std::array<unsigned,dim> idx;
				std::iota(idx.begin(), idx.end(), 0);

				std::sort(
								idx.begin(), idx.end(),
								[&](unsigned i1, unsigned i2) { return a[i1] > a[i2]; }
				);

				return idx;
}


#endif // MY_ARRAY_HPP
