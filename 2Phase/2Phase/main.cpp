#include "DiscreteProblem.hpp"
//#include "Linear/Intel_ILU0.hpp"
//#include "Linear/Intel_ILUT.hpp"
//#include "Linear/Intel_Prec_GMRES.hpp"
#include "Intel_Pardiso.hpp"
#include "NewtonSolver.hpp"

typedef GENSOL::Intel_Pardiso                         LINEARSOLVER;
//typedef GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 >  LINEARSOLVER;
typedef GENSOL::NewtonSolver< DiscreteProblem, LINEARSOLVER > STDN;

#include "fastl/containers/pod_vector_unbounded.hpp"
#include <fstream>

void
read_from_file( const char * filename, std::vector<double> &v )
{
   double tmp;
   std::ifstream strm ( filename );
   if ( strm.good( ) )
   {
      for ( std::size_t i=0; i<v.size(); ++i )
	 if ( strm >> tmp ) v[i] = tmp;
   }
}

void 
dump_solution( const char * filename, const DiscreteProblem::StateVector & v )
{
   std::ofstream strm( filename );
   for ( std::size_t i = 0; i<v.size(); ++i )
      strm << v[i].Po.value() << "\t"
	   << v[i].Sw.value() << std::endl;
   strm.close();
}

void 
dump_field( const char * filename, 
	    const std::vector<double> &phi,
	    const std::vector<double> &kx,
	    const std::vector<double> &ky,
	    const std::vector<double> &kz  )
{
   std::ofstream strm( filename );
   for ( std::size_t i = 0; i<phi.size(); ++i )
      strm << phi[i] << "\t"
	   << kx[i] << "\t"
	   << ky[i] << "\t"
	   << kz[i] << std::endl;
   strm.close();
}

int main ( )
{
   
   const std::size_t MAX_NLNITER = 6;
   const double      DT_INIT     = 5;
   const double      DT_CUT      = 0.5;
   const double      DT_GROW     = 2.0;
   const double      DT_MAX      = 70.0;
   const double      T_FINAL     = 100;
   
   const std::size_t NX=20, NY=20, NZ=1;
   const double      LX=800, LY=800, LZ=10.0;
   std::vector<double> vPORO ( NX*NY*NZ, 0.2 );
   std::vector<double> vKX(NX*NY*NZ, 200.33);
   std::vector<double> vKY(NX*NY*NZ, 200.33);
   std::vector<double> vKZ(NX*NY*NZ, 200.33);
   
   const double OPT_NLNITER = ( 3 < MAX_NLNITER ? MAX_NLNITER : 3 );
   DiscreteProblem model( NX, NY, NZ, LX, LY, LZ, vKX,vKY,vKZ, vPORO );

   LINEARSOLVER lnsolver( model.max_num_eqns(), 
			  model.max_num_nnz() );

   STDN newton( model, MAX_NLNITER, 1);

   DiscreteProblem::StateVector uOld, uNew;
   model.initialize_state( uOld );

   double DT   = DT_INIT;
   double time = 0.0;
   while ( time < T_FINAL )
   {
      uNew = uOld;
      STDN::report_t stdsmry = newton.solve_timestep( uNew, uOld, DT, model, lnsolver );

      if ( stdsmry.is_converged )
      {
         uOld =  uNew;
         time += DT;
         if ( stdsmry.niter < OPT_NLNITER ) DT *= DT_GROW;
         DT = std::min( DT, DT_MAX );
         DT = std::min( DT, T_FINAL - time );
         std::cout << "CONVERGED t = " << time << " days" << std::endl;
      }
      else
      {
         DT *= DT_CUT;
         std::cout << "FAILED " << std::endl;
      }
   }
   dump_solution( "./output/results.out", uNew );
   return -1;
}
