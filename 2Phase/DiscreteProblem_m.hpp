#pragma once

#include <vector>
#include "adetl/systems/ADvector.hpp"
#include "PropertyCalculator.hpp"
#include "DiscreteProblem.hpp"
//#include "DiscreteDomain.hpp"
#include <vector>
#include "CSR_Matrix.hpp"


typedef GENSOL::CSR_Matrix<double, int> CSR;
///////////////////////////////////////////////////////////////////////////////
class DiscreteProblem_m
{

public:
	typedef struct
	{
		bool        is_producer;
		std::size_t loc;
		double      WI;
		double	rw;
		std::vector<unsigned>	f_index;
		std::vector<double>		f_trans;
		std::vector<double>     QINJ;
		std::vector<double>      Pbh;
	}                                            Well;

	typedef struct
	{
		std::vector<double> val;
		std::vector<int> col, row;

	} Sparse_matrix;
	typedef struct
	{
		adetl::ADscalar<> qo = 0;
		adetl::ADscalar<> qw = 0;
		adetl::ADscalar<> P = 0;
		adetl::ADscalar<> WCT = 0;
		adetl::ADscalar<> qt = 0;
	}											wellrate;

	typedef EDFM::Hor_Well						Hor_well;
private:
	typedef PropertyCalculator::Phase     PhaseID;
	typedef PropertyCalculator::StatusID  StatusID;
	typedef struct
	{
		adetl::ADscalar<> Po;
		adetl::ADscalar<> Sw;
		StatusID   status;
	}                                            StateElement;


	typedef struct
	{
		double	phi;
		double	vol;
	}											fcell;
	typedef struct
	{
		unsigned ia;
		unsigned ib;
		double trans;
	}		connection;
	typedef EDFM								  FracType;
public:
	typedef DiscreteProblem::StateVector        StateVector;
	typedef struct
	{
		double MatBal[2];
		double NormSat[2];
	}                                            ConvergenceInfo;
private:
	typedef PropertyCalculator::CellProps        CellProps;
	typedef std::vector< CellProps >             CellVector;

	typedef struct
	{
		double            T;
		adetl::ADscalar<> L[2];
		adetl::ADscalar<> Pot[2];
	}                                           FaceProps;
	typedef std::vector< FaceProps >            FaceVector;

	typedef DiscreteDomain                      MeshType;

public:

	void transfer_data(EDFM& model);

	void record_data(std::vector<std::vector<double>> &qo,
		std::vector<std::vector<double>> &qw);
	void get_inj_index(std::vector<unsigned> &iidx);
	DiscreteProblem_m(std::size_t NX, std::size_t NY, std::size_t NZ,
		double LX, double LY, double LZ,
		const std::vector<double> &vKX,
		const std::vector<double> &vKY,
		const std::vector<double> &vKZ,
		const std::vector<double> &Phi_ref);
	DiscreteProblem_m(std::size_t NX, std::size_t NY, std::size_t NZ,
		double LX, double LY, double LZ,
		const std::vector<double> &vKX,
		const std::vector<double> &vKY,
		const std::vector<double> &vKZ,
		const std::vector<double> &Phi_ref,
		int model);
	std::size_t  max_num_eqns() const { return 2 * N_c; }
	std::size_t  max_num_nnz() const { return 7 * 2 * N_c; }

	void initialize_state(DiscreteProblem_m::StateVector &oldstate, DiscreteProblem_m::StateVector &newstate,
		DiscreteProblem_m::StateVector &aoldstate, DiscreteProblem_m::StateVector &anewstate);

	void bind_to_old_state(const StateVector &_old_state);

	bool discretize(const StateVector & state, double DT);



	bool is_converged(ConvergenceInfo & nrm);

	void load_hschedule(std::vector<std::vector<unsigned>> &hwells,  dlib::matrix<unsigned> &s);

	//template< typename V, typename M >
	//void extract_R_Jm(V &residual, M &jacobian, std::size_t i);

	void extract_R_J(std::vector<double> &r, CSR &m, std::size_t offset);

	double compute_total_WI(Hor_well &h_well);

	void dump_residual(std::ofstream &out)
	{
		fastl::spy(out, mResidual);
	}
	std::size_t read_well();

	/////////////////////////////////////////////////
	//void DiscreteProblem_m::extract_F_der(DiscreteProblem::StateVector &oldstate, DiscreteProblem::StateVector &state, double DT,
	//	std::vector<CSR> &F_y, std::vector<CSR> &F_yold, std::vector<CSR> &F_v);

	void DiscreteProblem_m::extract_obj_der(unsigned nsche,
		std::vector<adetl::ADscalar<>> &objv);

	void loadwell_data(const std::vector<double> &time, const dlib::matrix<double> &well_sche);

	void extract_F_der(DiscreteProblem_m::StateVector &oldstate, DiscreteProblem_m::StateVector &state, double DT, DiscreteProblem_m::Sparse_matrix &M );
	////////////////////////////////////////////////
protected:
	void setup_wells(std::size_t NX, std::size_t NY, std::size_t NZ,
		double LX, double LY, double LZ,
		const std::vector<double> &vKX,
		const std::vector<double> &vKY,
		const std::vector<double> &vKZ);

	void check_schedule_size(unsigned nw, unsigned nc, unsigned nwc);

	std::size_t eqnID(std::size_t c, PhaseID ph) const
	{
		std::size_t incr;
		switch (ph)
		{
		case PhaseID::O:
			incr = 0;
			break;
		case PhaseID::W:
			incr = 1;
			break;
		}
		return 2 * c + incr;
	}

	std::size_t eqnID(std::size_t c, std::size_t ph) const
	{
		PhaseID incr;
		switch (ph)
		{
		case 0:
			incr = PhaseID::O;
			break;
		case 1:
			incr = PhaseID::W;
			break;
		}
		return eqnID(c, incr);
	}
	std::size_t ftoc(std::size_t fracID) const
	{
		return fracID + N_m;
	}
	std::size_t ctof(std::size_t cellID) const
	{
		return cellID - N_m;
	}
	void initialize_transmissibility(const std::vector<double> & KX,
		const std::vector<double> & KY,
		const std::vector<double> & KZ);
	double safeguard_MAC(double upd);
	void compute_accumulation();
	void compute_flow();
	void compute_wells();
	void compute_face_properties();
	void compute_cell_properties(const StateVector & state);
	void compute_hwells();

	void extract_inform(const std::vector<std::vector<unsigned>> &x, std::vector<std::vector<unsigned>> &inf);

	////////////////////////////////////////////////////////////////////

	//bool discretizem(const StateVector & state, double DT);

	//void DiscreteProblem::compute_Qo();

	//void DiscreteProblem::compute_Qw();

	//void DiscreteProblem::compute_Pwf();

	//void DiscreteProblem::compute_WCT();

	//void DiscreteProblem::compute_objv();

	//template< typename V, typename M >
	//void extract_acc(V &r, M &m, std::size_t offset);

	//void assign_wellrate(std::vector<wellrate> &temp_v);

	void CSR_to_Sparse(CSR & x, Sparse_matrix &M);

	/////////////////////////////////////////////////////////////////////////



public:
	std::size_t			nw,Nj, Np, N_hj, N_hp, nct,nwc;
	std::vector<Well>	 mWells;
	std::vector<double>  schedule;
	double  sum_t;
	std::vector<std::vector<adetl::ADscalar<>>>   H_chedule;
	std::vector<std::vector<unsigned>>	H_constrain;
	std::vector<std::vector<unsigned>>	H_control;
	std::vector<wellrate>   Hwell_q;
	std::vector<std::vector<unsigned>> wct_inform, wcs_inform;
	double		row;
	std::vector<double> btt;
private:
	MeshType            mMesh;
	PropertyCalculator  mPropCalc;
	FracType            fracsys;
	std::size_t        N_f, N_m, N_mm, N_fm, N_ff, N_var, N_c, face_n;
	CellVector          mCells;
	std::vector<fcell>	   fCells;
	FaceVector          mFaces;
	FaceVector          fFaces;
	std::vector<connection> mmc;
	std::vector<connection> mfc;
	std::vector<connection> ffc;
	std::vector<double> mPhi_ref;
	std::vector<double> mAccum_old;
	adetl::ADvector     mAccum;
	adetl::ADvector     mFlow;
	adetl::ADvector     mResidual;
	adetl::ADvector     mResidualm;
	double              mDT;
	adetl::ADscalar<>   mTmpVars[3];
	std::size_t			nsche;


	/////////////////extension for adjoint code//////////////////////////////////
	std::vector<wellrate> temp_v, temp_y;

	////////////////////////////////////////////////////////////////////////////
};