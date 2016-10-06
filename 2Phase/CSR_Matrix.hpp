// -------------------------------- *- C++ -* ------------------------------ //
// -------------------------------* TEMPLATE *------------------------------ //
// ------------------------------------------------------------------------- //
//! \file  CSR_Matrix.hpp
//! \brief Compressed Sparse Row Storage Matrix datatype
//
// ------------------------------------------------------------------------- //
/** *************************************************************************
 *  \author     : Rami M. Younis
 *  \created    : 07/22/2013
 *  \revised    : 
 *  \warning    : 
 *  Target MCU  : Generic
 *  Editor TABs : 3
 *  Editor      : emacs
 *  Auto Style  : ellemtel
 ** ************************************************************************ */
/** **************************************************************************
 * Copyright (c) 2013, all rights reserved
 * FUture Reservoir Simulation Systems & Technology
 * McDougall School of Petroleum Engineering
 * The University of Tulsa, Tulsa, Oklahoma 74104, USA
 ** ************************************************************************ */
#ifndef __CSRMATRIX_HPP_INCLUDED_
#define __CSRMATRIX_HPP_INCLUDED_

#include <vector>
#include <ostream>
#include "mkl.h"
namespace GENSOL { // -------------------------------------------- BEGIN NAMESPACE 

  class EllipticInfo
  {
  public: 
    std::size_t elliptic_var_id;
    std::size_t n_vars_per_block;
    std::size_t n_blocks;
    std::size_t n_unblocked_vars;
  };

  // ----------------------------  CSR_Matrix  ------------------------------ //
  /** \class CSR_Matrix
   *  
   *  Compressed Sparse Row Storage Matrix datatype
   **/
  // ------------------------------------------------------------------------- //
  template< typename __T, typename __I >
  class CSR_Matrix
  {
  public:
    typedef __T *                             T_pointer;
    typedef __I *                             I_pointer;

  private:
    typedef std::vector<__T>                  T_collection;
    typedef std::vector<__I>                  I_collection;

  public:
    //.............................  LIFECYCLE  ..........................//
    CSR_Matrix(  ) : m_val( ), m_colind( ), m_rowptr( ) {  }

    CSR_Matrix(  __I _N, __I _NNZ ) { resize(_N,_NNZ); }

      
    //...........................  ASSIGNMENT  ...........................//      
    
    //.............................  ACCESS  .............................//
    __I N() const { return mN; }
    __I NNZ() const { return mNNZ; }
    __I offset() const { return m_rowptr[0]; }

    const T_collection & value()  const { return m_val; }
    const I_collection & colind() const { return m_colind; }
    const I_collection & rowptr() const { return m_rowptr; }
    __T                  operator()( std::size_t i, std::size_t j ) const
    {
      const __I offset = m_rowptr[0];
      double    value  = 0;
      if ( (i<mN) && (j<mN) ) 
	{
	  const __I r1 = m_rowptr[i] - offset;
	  const __I r2 = m_rowptr[i+1] - offset - 1;
	  if ( ( j >= (m_colind[r1]-offset) ) && ( j <= ( m_colind[r2]-offset ) ) )
	    {
	      __I r = r1;
	      while ( j > (m_colind[r]-offset) ) ++r;
	      if ( j == ( m_colind[r]-offset) ) value = m_val[r];
	    }
	}
      return value;
    };

    T_collection & value()  { return m_val; }
    I_collection & colind() { return m_colind; }
    I_collection & rowptr() { return m_rowptr; }   

    //.............................  OPERATORS  .........................//
    void resize( __I _N, __I _NNZ )
    {
      mN = _N; mNNZ = _NNZ;
      m_val.resize( mNNZ );
      m_colind.resize( mNNZ );
      m_rowptr.resize( mN + 1 );
    }
   
    void check_size( )
    {
      mN   = m_rowptr.size()-1;
      mNNZ = m_val.size();
    }
	CSR_Matrix transpose()
	{
		int job[8];
		job[0] = 0;
		job[1] = 0;
		job[2] = 0;
		int inform = 0;
		int N = this->N();
		CSR_Matrix B(this->N(), this->NNZ());
		mkl_dcsrcsc(job, &N, this->value().data(), this->colind().data(), this->rowptr().data(),
			B.value().data(), B.colind().data(), B.rowptr().data(), &inform);
		return B;
	}

	std::vector<double> operator*(std::vector<double> & V)
	{
		//char matdes[]= "G  C  ";
		char tran = 'N';
		I_collection pointE;
		I_collection pointB;
		double alpha = 1;
		double beta = 0;
		int m = this->mN;
		std::vector<double> r(m);
		pointB.resize(m);
		pointE.resize(m);
		for (unsigned i = 0; i < m; i++)
		{
			pointB[i] = this->rowptr()[i];
			pointE[i] = this->rowptr()[i + 1];
		}
		mkl_dcsrmv(&tran, &m, &m, &alpha, "G**C",
			this->value().data(), this->colind().data(), pointB.data(), pointE.data(), V.data(), &beta,
			r.data());
		return r;
	}
	std::vector<double> tranpose_product(std::vector<double> & V)
	{
		//char matdes[]= "G  C  ";
		char tran = 'T';
		I_collection pointE;
		I_collection pointB;
		double alpha = 1;
		double beta = 0;
		int m = this->mN;
		std::vector<double> r(m);
		pointB.resize(m);
		pointE.resize(m);
		for (unsigned i = 0; i < m; i++)
		{
			pointB[i] = this->rowptr()[i];
			pointE[i] = this->rowptr()[i + 1];
		}
		mkl_dcsrmv(&tran, &m, &m, &alpha, "G**C",
			this->value().data(), this->colind().data(), pointB.data(), pointE.data(), V.data(), &beta,
			r.data());
		return r;
	}
	std::vector<double> minus_tranpose_product(std::vector<double> & V)
	{
		//char matdes[]= "G  C  ";
		char tran = 'T';
		I_collection pointE;
		I_collection pointB;
		double alpha = -1;
		double beta = 0;
		int m = this->mN;
		std::vector<double> r(m);
		pointB.resize(m);
		pointE.resize(m);
		for (unsigned i = 0; i < m; i++)
		{
			pointB[i] = this->rowptr()[i];
			pointE[i] = this->rowptr()[i + 1];
		}
		mkl_dcsrmv(&tran, &m, &m, &alpha, "G**C",
			this->value().data(), this->colind().data(), pointB.data(), pointE.data(), V.data(), &beta,
			r.data());
		return r;
	}

	void m_m_multiply(CSR_Matrix & M)
	{
		//char matdes[]= "G  C  ";
		char tran = 'N';
		I_collection pointE;
		I_collection pointB;
		int m = n = k = M.N();
		double alpha = -1;
		double beta = -1;
		std::vector<double> r(m);
		pointB.resize(m);
		pointE.resize(m);
		for (unsigned i = 0; i < m; i++)
		{
			pointB[i] = this->rowptr()[i];
			pointE[i] = this->rowptr()[i + 1];
		}
		mkl_dcsrmv(&tran, &m, &m, &alpha, "G**C",
			this->value().data(), this->colind().data(), pointB.data(), pointE.data(), V.data(), &beta,
			r.data());
		return r;
	}

	int
		CSR_to_dens(double *b)
	{
			int job[6];
			job[0] = 1;
			job[1] = 0;
			job[2] = 0;
			job[3] = 2;
			job[4] = this->NNZ();
			int m = this->mN, n = this->mN;

			int lda = m;
			int info;
			mkl_ddnscsr(job, &m, &n, b,
				&lda, this->value().data(), this->colind().data(), this->rowptr().data(), &info);
			return info;
		}
	int
		dens_to_CSR(double *b)
	{
			int job[6];
			job[0] = 0;
			job[1] = 0;
			job[2] = 0;
			job[3] = 2;
			job[4] = this->NNZ();
			int m = this->mN, n = this->mN;
			int lda = 0;
			int info;
			mkl_ddnscsr(&job, &m, &n, b,
				&lda, this->value().data(), this->colind().data(), this->rowptr().data(), &info);
			return info;
		}
  public:
    EllipticInfo ainfo;
  private:
    __I mN;
    __I mNNZ;
    T_collection m_val;
    I_collection m_colind;
    I_collection m_rowptr;
  }; // CSR_Matrix

  template< typename A, typename B >
  std::ostream & spy ( std::ostream & ostr, 
		       const GENSOL::CSR_Matrix<A,B> & _out )
  {
	  ostr << "[";
    for ( std::size_t r = 0; r < _out.N( ); ++r )
      {

	const std::size_t I1  = _out.rowptr()[r] - 1;
	const std::size_t I2  = _out.rowptr()[r+1] - 1;
	std::size_t cc  = 0;
	for ( std::size_t nz = I1; nz < I2; ++nz )
	  {
	    const std::size_t c = _out.colind()[nz];
	    for ( ; cc < c; ++cc ) ostr << 0;
	    if (_out.value()[nz]!=0.0) ostr <<1;
	    ++cc;
	  }
	ostr << std::endl;
      }
	ostr << "]" << std::endl;
    return ostr;
  };

}; // ---------------------------------------------------------- END NAMESPACE

template< typename A, typename B >
std::ostream & operator << ( std::ostream & ostr, 
			     const GENSOL::CSR_Matrix<A,B> & _out )
{
  ostr << _out.N() << " X " << _out.N() << "\t NNZ = " << _out.NNZ() << std::endl;
  ostr << "VAL = [";
  for (std::size_t i = 0; i<_out.NNZ(); ++i )
    ostr << _out.value()[i] << ", ";
  ostr << " ]" << std::endl;
  ostr << "COL = [";
  for (std::size_t i = 0; i<_out.NNZ(); ++i )
    ostr << _out.colind()[i] << ", ";
  ostr << " ]" << std::endl;
  ostr << "ROW = [";
  for (std::size_t i = 0; i<_out.N()+1; ++i )
    ostr << _out.rowptr()[i] << ", ";
  ostr << " ]";
  return ostr;
}



#endif // __CSR_Matrix_HPP_INCLUDED_

//---------------------------------------------------------------------------//
//                           EOF CSR_Matrix.hpp
//---------------------------------------------------------------------------//

