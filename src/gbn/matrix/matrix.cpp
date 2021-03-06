#include "matrix.h"
#include <iostream>
#include "../../helpers.hpp"

// abstract Matrix class
Matrix::Matrix(const Index n, const Index m, const MatrixType type, bool is_stochastic)
: n(n), m(m), type(type), is_stochastic(is_stochastic), one_mask_n(init_one_mask(n)), one_mask_m(init_one_mask(m))
{}

Matrix::~Matrix()
{}

BitVec Matrix::init_one_mask(Index n) {
	BitVec b;
	for(Index i = n; i < MAX_PLACES; i++)
		b[i] = 1;

	return b;
}


// DynamicMatrix
DynamicMatrix::DynamicMatrix(const Index n, const Index m)
	: Matrix(n,m,DYNAMIC,true)
{}

DynamicMatrix::~DynamicMatrix()
{}

double DynamicMatrix::get(const BitVec& to, const BitVec& from) const 
{
	auto it = data.find({ to, from });
	if(it == data.cend())
		return 0;
	else
		return it->second;
}

void DynamicMatrix::set(const BitVec& to, const BitVec& from, double val) 
{
	data[{ to, from }] = val;
}

void DynamicMatrix::add(const BitVec& to, const BitVec& from, double val) 
{
	data[{ to, from }] += val;
}


// FMatrix
FMatrix::FMatrix(const Index k, const bool b)
	: Matrix(k,k,F,false), k(k), b(b)
{
	if(n != m)
		throw std::logic_error("Tried to create F matrix where n != m");	
}

FMatrix::~FMatrix()
{}

double FMatrix::get(const BitVec& to, const BitVec& from) const 
{
	auto to_copy = to | one_mask_m;

	if(to != from)
		return 0;

	if((b && to_copy.all()) || (!b && to.none()))
		return 0;

	return 1;
}

void FMatrix::set(const BitVec& /*to*/, const BitVec& /*from*/, double /*val*/) 
{
	throw std::logic_error("Tried to set entry of F-matrix.");
}

void FMatrix::add(const BitVec& /*to*/, const BitVec& /*from*/, double /*val*/) 
{
	throw std::logic_error("Tried to set entry of F-matrix.");
}



// OneBMatrix
OneBMatrix::OneBMatrix(const bool b)
	: Matrix(0,1,ONE_B,true), b(b)
{}

OneBMatrix::~OneBMatrix()
{}

double OneBMatrix::get(const BitVec& to, const BitVec& /*from*/) const 
{
	return (to.test(0) == b) ? 1 : 0;
}

void OneBMatrix::set(const BitVec& /*to*/, const BitVec& /*from*/, double /*val*/) 
{
	throw std::logic_error("Tried to set entry of OneB-matrix.");
}

void OneBMatrix::add(const BitVec& /*to*/, const BitVec& /*from*/, double /*val*/) 
{
	throw std::logic_error("Tried to set entry of OneB-matrix.");
}




bool is_stochastic(const Matrix& matrix)
{
	unsigned long long i_max_row = 1;
	unsigned long long i_max_col = 1;
	i_max_col = i_max_col << matrix.n;
	i_max_row = i_max_row << matrix.m;

	for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
	{
		double sum = 0;
		for(unsigned long long i_row = 0; i_row < i_max_row; i_row++)
		{
			sum += matrix.get(i_row, i_col);
		}
		if(std::abs(1-sum) > 1e-5)
			return false;
	}

	return true;
}




// TerminatorMatrix
TerminatorMatrix::TerminatorMatrix()
	:Matrix(1,0,TERMINATOR,true)
{}

TerminatorMatrix::~TerminatorMatrix()
{}

double TerminatorMatrix::get(const BitVec& /*to*/, const BitVec& /*from*/) const
{
	return 1;
}

void TerminatorMatrix::set(const BitVec& /*to*/, const BitVec& /*from*/, double /*val*/)
{
	throw std::logic_error("Tried to set entry of TerminatorMatrix.");
}
void TerminatorMatrix::add(const BitVec& /*to*/, const BitVec& /*from*/, double /*val*/)
{
	throw std::logic_error("Tried to set entry of TerminatorMatrix.");
}




// ZeroMatrix
ZeroMatrix::ZeroMatrix(const Index n, const Index m)
	:Matrix(n,m,ZERO,false)
{}

ZeroMatrix::~ZeroMatrix()
{}

double ZeroMatrix::get(const BitVec& /*to*/, const BitVec& /*from*/) const
{
	return 0;
}

void ZeroMatrix::set(const BitVec& /*to*/, const BitVec& /*from*/, double /*val*/)
{
	throw std::logic_error("Tried to set entry of ZeroMatrix.");
}
void ZeroMatrix::add(const BitVec& /*to*/, const BitVec& /*from*/, double /*val*/)
{
	throw std::logic_error("Tried to set entry of ZeroMatrix.");
}

// DiagonalMatrix
DiagonalMatrix::DiagonalMatrix(const Index k)
		: Matrix(k,k,DIAGONAL,false), k(k)
{
	if(n != m)
		throw std::logic_error("Tried to create DiagonalMatrix where n != m");
}

DiagonalMatrix::~DiagonalMatrix()
{}

double DiagonalMatrix::get(const BitVec& to, const BitVec& from) const
{
	if(to != from)
		return 0;

	for(auto d : data) {
	    if(d.first == to) return d.second;
	}
	return 0;
}

void DiagonalMatrix::set(const BitVec& to, const BitVec& from, double val)
{
	if(to != from)
		throw std::logic_error("Tried to set entry of DiagonalMatrix on n != m");

    data[to] = val;
}

void DiagonalMatrix::add(const BitVec& to, const BitVec& from, double val)
{
    data[to] += val;
}