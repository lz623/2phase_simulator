#include "DiscreteProblem.hpp"
#include "MeshBuilder.hpp"
#include "Point.hpp"
#include <fstream>

DiscreteProblem::DiscreteProblem( std::size_t NX, std::size_t NY, std::size_t NZ,
				  double LX, double LY, double LZ,
				  const std::vector<double> &vKX, const std::vector<double> &vKY, const std::vector<double> &vKZ,
				  const std::vector<double> &Phi_ref) :
   mMesh( MeshBuilder<UniformCartesian> ( NX,NY,NZ,LX,LY,LZ ) ),
   mPropCalc(  ),
   fracsys(NX, NY, NZ, LX, LY, LZ, vKX, Phi_ref),
   mFaces( mMesh.size_faces() ),
   mWells( ),
   mPhi_ref( Phi_ref )
{
	transfer_data(fracsys);
	mCells.resize(N_c);
	fFaces.resize(N_fm+N_ff);
	mAccum_old.resize(2 * N_c);
	mAccum.resize(2 * N_c, 0);
	mFlow.resize(2 * N_c, 0);
	mResidual.resize(N_tc, 0);
   initialize_transmissibility( vKX, vKY, vKZ );
   setup_wells( NX,NY,NZ, LX,LY,LZ, vKX, vKY, vKZ );
};
DiscreteProblem::DiscreteProblem(std::size_t NX, std::size_t NY, std::size_t NZ,
	double LX, double LY, double LZ,
	const std::vector<double> &vKX, const std::vector<double> &vKY, const std::vector<double> &vKZ,
	const std::vector<double> &Phi_ref,int model) :
	mMesh(MeshBuilder<UniformCartesian>(NX, NY, NZ, LX, LY, LZ)),
	mPropCalc(),
	fracsys(NX, NY, NZ, LX, LY, LZ, vKX, Phi_ref,model),
	mFaces(mMesh.size_faces()),
	mWells(),
	mPhi_ref(Phi_ref)
{
	transfer_data(fracsys);
	mCells.resize(N_c);
	fFaces.resize(N_fm + N_ff);
	mAccum_old.resize(2 * N_c);
	mAccum.resize(2 * N_c, 0);
	mFlow.resize(2 * N_c, 0);
	mResidual.resize(2 * N_c, 0);
	initialize_transmissibility(vKX, vKY, vKZ);
	setup_wells(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ);
};

void DiscreteProblem::transfer_data(EDFM& model)
{
	unsigned n = model.inter_nn.coor.size();
	connection tmp;
	N_fm = model.n_frac.size();
	N_f = N_fm;
	N_ff = model.inter_nn.i_index.size();
	N_m = mMesh.size_cells();
	N_c = N_fm + N_m;
	N_tc = 2*N_c + fracsys.N_w;
	fCells.resize(N_fm);
	nw = fracsys.N_w;
	for (unsigned i = 0; i < N_fm; i++)
	{
		tmp.ia = model.n_frac[i].g_index;
		tmp.ib = model.n_frac[i].f_index;
		tmp.trans = 0.00112712*model.n_frac[i].trans;
		fCells[i].phi = model.n_frac[i].por;
		fCells[i].vol = model.n_frac[i].A*model.n_frac[i].aper;
		mfc.push_back(tmp);
	}

	for (unsigned i = 0; i < N_ff; i++)
	{
		tmp.ia = model.inter_nn.i_index[i];
		tmp.ib = model.inter_nn.j_index[i];
		tmp.trans = 0.00112712 *model.inter_nn.trans[i];
		ffc.push_back(tmp);
	}
	is_producer.resize(nw);
	well_rate.resize(nw);
	for (unsigned i = 0; i < fracsys.h_well.size();i++)
	{
		is_producer[i]=fracsys.h_well[i].isproducer;
	}

}

void DiscreteProblem::extract_inform(const std::vector<std::vector<unsigned>> &x, std::vector<std::vector<unsigned>> &inf)
{
	unsigned key;
	inf.resize(nw);
	for (unsigned i = 0; i < nw;i++)
	{
		inf[i].resize(5);
		std::fill(inf[i].begin(), inf[i].end(), 0);
		for (unsigned j = 0; j <x[i].size();j++)
		{
			key = x[i][j];
			inf[i][key] = j;
		}
	}
}

void
DiscreteProblem::read_from_file_unsigned(std::string filename, dlib::matrix<unsigned> &v)
{
	double tmp;
	size_t n = 0;
	std::ifstream strm(filename);
	while (strm >> tmp)
	{
		n++;
	}
	v.set_size(n, 1);
	strm.close();
	strm.open(filename);
	for (unsigned i = 0; i < n; i++)
	{
		strm >> v(i, 0);
	}
}

void DiscreteProblem::loadwell_independent(const std::vector<double> &time, const dlib::matrix<double> &well_sche)
{
	dlib::matrix<unsigned> ct, cs;  // input well constrain and well control
	read_from_file_unsigned("input/wellcontrol.dat", ct);
	read_from_file_unsigned("input/wellconstrain.dat", cs);
	check_schedule_size(nw, time.size(), well_sche.size());
	schedule = time;
	nct = time.size();
	nwc = well_sche.size();
	for (unsigned i = 0; i < nw; i++)
	{
		H_chedule[i].resize(time.size());
		for (unsigned j = 0; j < time.size(); j++)
		{
			H_chedule[i][j].value() = well_sche(i*time.size() + j);
			H_chedule[i][j].make_independent(i*time.size() + j);
		}
	}
	load_hschedule(H_control, ct);
	load_constrain();
}

void DiscreteProblem::load_constrain()
{
	std::ifstream input("input/wellconstrain.dat");
	H_constrain.resize(nct);
	for (unsigned w = 0; w < nw; w++)
	{
		for (unsigned i = 0; i < nct; i++)
		{
			H_constrain[i].resize(nw);
			input >> H_constrain[i][w].BHP_i;
			input >> H_constrain[i][w].qt_i;
			input >> H_constrain[i][w].qo_i;
			input >> H_constrain[i][w].qw_i;
			input >> H_constrain[i][w].wct_i;
		}
	}
	compute_Num_cs();
}

void DiscreteProblem::compute_Num_cs()
{
	Num_cs = 0; 
	for (unsigned w = 0; w < nw; w++)
	{
		for (unsigned i = 0; i < nct; i++)
		{
			(H_constrain[i][w].BHP_i) ? Num_cs++ : Num_cs;
			(H_constrain[i][w].qt_i) ? Num_cs++ : Num_cs;
			(H_constrain[i][w].qo_i) ? Num_cs++ : Num_cs;
			(H_constrain[i][w].qw_i) ? Num_cs++ : Num_cs;
			(H_constrain[i][w].wct_i) ? Num_cs++ : Num_cs;
		}
	}
}

void DiscreteProblem::loadwell_data(const std::vector<double> &time,const dlib::matrix<double> &well_sche)
{
	dlib::matrix<unsigned> ct, cs;  // input well constrain and well control
	double min_inj_rate(-0.001); // minimum injecting rate for inejector
	read_from_file_unsigned("input/wellcontrol.dat", ct);
	check_schedule_size(nw, time.size(), well_sche.size());
	schedule = time;
	nwc = ct.size();
	nct = time.size();
	load_hschedule(H_control, ct);
	load_constrain();
	for (unsigned i = 0; i < Np; i++)
	{
		H_chedule[i].resize(time.size());
		for (unsigned j = 0; j < time.size(); j++)
		{
			H_chedule[i][j] = well_sche(i*time.size() + j);
		}
	}
	for (unsigned i = Np; i < nw; i++)
	{
		H_chedule[i].resize(time.size());
		for (unsigned j = 0; j < time.size(); j++)
		{
			//set injecting rate for injector small than 0
			if (H_chedule[i][j]>0 && H_control[i][j]==3){
				H_chedule[i][j] = -min_inj_rate;
				dlib::matrix<double>& w_sche= const_cast<dlib::matrix<double> &> (well_sche);
				w_sche(i*time.size() + j) = -min_inj_rate;
			}
			else{
				H_chedule[i][j] = well_sche(i*time.size() + j);
			}
		}
	}
}

void DiscreteProblem::check_schedule_size(unsigned nw, unsigned nct, unsigned nwc)
{
	if (nw*nct != nwc)
	{
		std::cout << "The schedule's size is not match" << std::endl;
		system("pause");
	}

}



void
DiscreteProblem::setup_wells(std::size_t NX, std::size_t NY, std::size_t NZ, 
			     double LX, double LY, double LZ,
			     const std::vector<double> &vKX, 
			     const std::vector<double> &vKY, 
			     const std::vector<double> &vKZ  )
{
	const double DX = LX / static_cast<double>(NX);
	const double DY = LY / static_cast<double>(NY);
	const double DZ = LZ / static_cast<double>(NZ);
	N_hp=fracsys.np;
	N_hj = fracsys.nj;
	H_chedule.resize(fracsys.N_w);
	Hwell_q.resize(fracsys.N_w);
	for (std::size_t i = 0; i <fracsys.N_w;i++)
	{
		if (fracsys.h_well[i].isvertical)
		{
			for (std::size_t w = 0; w <fracsys.h_well[i].gridindex.size(); w++)
			{
				unsigned W_loc = fracsys.h_well[i].gridindex[w];
				const double Kh = DZ * sqrt(vKY[W_loc] * vKX[W_loc]);
				const double r1 = vKY[W_loc] / vKX[W_loc];
				const double r2 = vKX[W_loc] / vKY[W_loc];
				const double ro = 0.28 * std::sqrt(std::sqrt(r1)*DX*DX + std::sqrt(r2)*DY*DY) / (std::pow(r1, 0.25) + std::pow(r2, 0.25));
				fracsys.h_well[i].WI.push_back(0.00708*Kh / log(ro / fracsys.h_well[i].rw));
			}
		}
		else
		{
			for (std::size_t w = 0; w <fracsys.h_well[i].gridindex.size(); w++)
			{
				unsigned W_loc = fracsys.h_well[i].gridindex[w];
				const double Kh = DX * sqrt(vKY[W_loc] * vKZ[W_loc]);
				const double r1 = vKY[W_loc] / vKZ[W_loc];
				const double r2 = vKZ[W_loc] / vKY[W_loc];
				const double ro = 0.28 * std::sqrt(std::sqrt(r1)*DZ*DZ + std::sqrt(r2)*DY*DY) / (std::pow(r1, 0.25) + std::pow(r2, 0.25));
				fracsys.h_well[i].WI.push_back(0.00708*Kh / log(ro / fracsys.h_well[i].rw));
			}
		}
	}
}

void 
DiscreteProblem::initialize_state( DiscreteProblem::StateVector &state )
{
   state.resize( N_c );
   std::ifstream strm("./input/inital_state.dat");
   double Pi, Sw, Pif, Swf = 0;
   strm >> Pi >> Sw >> Pif >> Swf;
   for (std::size_t c = 0; c < mMesh.size_cells(); ++c)
   {
	   state[c].status = StatusID::OW;
	   state[c].Po = Pi;
	   state[c].Po.make_independent(eqnID(c, PhaseID::O));
	   state[c].Sw = Sw;
	   state[c].Sw.make_independent(eqnID(c, PhaseID::W));
   }

   for (std::size_t f = 0; f < N_f; ++f)
   {
	   std::size_t c = ftoc(f);
	   state[c].status = StatusID::OW;
	   state[c].Po = Pif;
	   state[c].Po.make_independent(eqnID(c, PhaseID::O));
	   state[c].Sw = Swf;
	   state[c].Sw.make_independent(eqnID(c, PhaseID::W));
   }
   for (std::size_t w = 0; w < nw;w++)
   {
	   std::size_t c = wtoc(w);
	   Hwell_q[w].P= Pi;
	   Hwell_q[w].P.make_independent(c);
   }
}

void
DiscreteProblem::initialize_state(DiscreteProblem::StateVector &uold, DiscreteProblem::StateVector &unew,
DiscreteProblem::StateVector &aold, DiscreteProblem::StateVector &anew)
{
	for (unsigned c = 0; c< N_c; c++)
	{
		aold[c].status = uold[c].status;
		anew[c].status = unew[c].status;
		aold[c].Po = uold[c].Po.value();
		anew[c].Po = unew[c].Po.value();
		aold[c].Sw = uold[c].Sw.value();
		anew[c].Sw = unew[c].Sw.value();
	}
}



void 
DiscreteProblem::bind_to_old_state( const DiscreteProblem::StateVector &old_state )
{
   compute_cell_properties( old_state );
   compute_accumulation( );
   for ( std::size_t c=0; c<N_c; ++c )
   {
      std::size_t eqn1 = eqnID(c,PhaseID::O);
      std::size_t eqn2 = eqnID(c,PhaseID::W);
      mAccum_old[ eqn1 ] = mAccum[ eqn1 ].value();
      mAccum_old[ eqn2 ] = mAccum[ eqn2 ].value();
   }
}

bool
DiscreteProblem::discretize( const DiscreteProblem::StateVector &state, double DT )
{
   bool is_badvalue = false;
   mDT = DT;
   compute_cell_properties( state );
   compute_accumulation( );
   //compute_BHP();
   compute_flow( );
   for ( std::size_t c=0; c<N_c; ++c )
   {
      std::size_t eqn1 = eqnID(c,PhaseID::O);
      std::size_t eqn2 = eqnID(c,PhaseID::W);
      mResidual[eqn1] = (mAccum[eqn1] - mAccum_old[eqn1]) + mDT * mFlow[eqn1];
      mResidual[eqn2] = (mAccum[eqn2] - mAccum_old[eqn2]) + mDT * mFlow[eqn2];
	  if (!std::isfinite(mResidual[eqn1].value()))
	  {
		  is_badvalue = true;
	  }
      if ( !std::isfinite( mResidual[eqn2].value() ) ) is_badvalue = true;
   }   


   for (std::size_t w = 0; w<nw; ++w)
   {
	   std::size_t eqn1 = wtoc(w);
	   if (H_control[w][nsche]==0)
	   {
		   mResidual[eqn1] = H_chedule[w][nsche] - Hwell_q[w].P;
	   }
	   else if (H_control[w][nsche] == 1)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qt);
	   }
	   else if (H_control[w][nsche] == 2)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qo);
	   }
	   else if (H_control[w][nsche] == 3)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qw);
	   }
	   else if (H_control[w][nsche] == 4)
	   {
		   mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].WCT);
	   }
   }

   return is_badvalue;
}
bool
DiscreteProblem::discretize_m(const DiscreteProblem::StateVector &state, double DT)
{
	bool is_badvalue = false;
	mDT = DT;
	compute_cell_properties(state);
	compute_accumulation();
	compute_BHP();
	compute_flow();
	for (std::size_t c = 0; c<N_c; ++c)
	{
		std::size_t eqn1 = eqnID(c, PhaseID::O);
		std::size_t eqn2 = eqnID(c, PhaseID::W);
		mResidual[eqn1] = (mAccum[eqn1] - mAccum_old[eqn1]) + mDT * mFlow[eqn1];
		mResidual[eqn2] = (mAccum[eqn2] - mAccum_old[eqn2]) + mDT * mFlow[eqn2];
		if (!std::isfinite(mResidual[eqn1].value())) is_badvalue = true;
		if (!std::isfinite(mResidual[eqn2].value())) is_badvalue = true;
	}
	for (std::size_t w = 0; w<nw; ++w)
	{
		std::size_t eqn1 = wtoc(w);
		if (H_control[w][nsche] == 0)
		{
			mResidual[eqn1] = H_chedule[w][nsche] - Hwell_q[w].P;
			//std::cout << mResidual[eqn1]<< std::endl;
		}
		else if (H_control[w][nsche] == 1)
		{
			mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qt);
		}
		else if (H_control[w][nsche] == 2)
		{
			mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qo);
		}
		else if (H_control[w][nsche] == 3)
		{
			mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].qw);
			//std::cout << mResidual[eqn1] << std::endl;
		}
		else if (H_control[w][nsche] == 4)
		{
			mResidual[eqn1] = mDT*(H_chedule[w][nsche] - Hwell_q[w].WCT);
		}
	}

	return is_badvalue;
}


bool
DiscreteProblem::is_converged ( ConvergenceInfo & nrm )
{
   bool is_MATBAL_converged = true;
   bool is_NRMSAT_converged = true;
   bool is_well_converged = true;
   double tot_PV = 0.0;
   for ( std::size_t c=0; c<mMesh.size_cells( ); ++c )
      tot_PV += mCells[c].phi.value() * mMesh.cell_measure( c )*0.1781076;
   for (std::size_t f= 0; f<N_f; ++f)
   {
	   std::size_t c = ftoc(f);
	   tot_PV += mCells[c].phi.value() * fCells[f].vol*0.1781076;
   }

   double sum_well = 0;
   for (std::size_t w = 0; w < nw; ++w)
   {
	   std::size_t c = wtoc(w);
	   sum_well += abs(mResidual[c].value());
   }
   if (abs(sum_well)> 1.0e-5) is_well_converged = false;
	 for ( std::size_t phs = 0; phs < 2; ++phs )
	 {
		 double sum_R = 0.0;
		 double avg_B = 0.0;
		 double max_R_PV = 0.0;

		 for ( std::size_t c=0; c<mMesh.size_cells( ); ++c )
		  {
			sum_R    += mResidual[ eqnID( c, phs ) ].value();
			avg_B    += 1.0 / mCells[c].b[ phs ].value();
			double R_PV = std::abs(mResidual[ eqnID( c, phs ) ].value() / (mCells[c].phi.value() * mMesh.cell_measure( c )*0.1781076 ));
			if (R_PV > max_R_PV) max_R_PV = R_PV;
			}
		 for (std::size_t f = 0; f<N_f; ++f)
		 {
			 std::size_t c = ftoc(f);
			 sum_R += mResidual[eqnID(c, phs)].value();
			 avg_B += 1.0 / mCells[c].b[phs].value();
			 double R_PV = std::abs(mResidual[eqnID(c, phs)].value() / (mCells[c].phi.value() * fCells[f].vol*0.1781076));
			 if (R_PV > max_R_PV) max_R_PV = R_PV;
		 }
		avg_B /= N_c;
		nrm.MatBal[ phs ]  = std::abs( avg_B * sum_R / tot_PV );
		nrm.NormSat[ phs ] = avg_B * max_R_PV;
		if ( nrm.MatBal[ phs ] > 1.0e-7 )  is_MATBAL_converged = false;
		if ( nrm.NormSat[ phs ] > 0.001 )  is_NRMSAT_converged = false;
	}
	 return (is_MATBAL_converged && is_NRMSAT_converged && is_well_converged);
}

void 
DiscreteProblem::initialize_transmissibility( const std::vector<double> & KX,
					     const std::vector<double> & KY,
					     const std::vector<double> & KZ )
{
   Point KL, KR;
   for ( std::size_t f=0; f < mMesh.size_faces( ) ; ++f )
   {
      const std::size_t c1 = mMesh.face_adj_cell_1( f );
      const std::size_t c2 = mMesh.face_adj_cell_2( f );
      KL.p[0] = KX[c1];
      KL.p[1] = KY[c1];
      KL.p[2] = KZ[c1];
      KR.p[0] = KX[c2];
      KR.p[1] = KY[c2];
      KR.p[2] = KZ[c2];
      const Point nrml = mMesh.unit_normal( f );
      const double kl = dot_product( nrml, KL );
      const double kr = dot_product( nrml, KR );
      const double DX = norm( mMesh.cell_coord(c1) - mMesh.cell_coord( c2 ) );
      const double A  = mMesh.face_measure( f );
      const double beta = DX*0.5/A*( 1.0/kl + 1.0/kr );
      mFaces[f].T = 0.00112712 / beta;
   }
};

void 
DiscreteProblem::compute_cell_properties( const DiscreteProblem::StateVector &state )
{
   for (std::size_t c = 0; c < mMesh.size_cells(); ++c ) 
   {
      mPropCalc.calculate( mPhi_ref[c], state[c], mCells[c] );
   }
   for (std::size_t f = 0; f < N_f; ++f)
   {
	   std::size_t c = ftoc(f);
	   mPropCalc.calculate_f(fCells[f].phi, state[c], mCells[c]);
   }
}

void 
DiscreteProblem::compute_accumulation( )
{
   for (std::size_t c = 0; c < mMesh.size_cells(); ++c ) 
   {
      const double VOLUME = mMesh.cell_measure( c )*0.1781076;
      mAccum[ eqnID(c,PhaseID::O) ] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::O] * mCells[c].b[PhaseID::O];
      mAccum[ eqnID(c,PhaseID::W) ] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::W] * mCells[c].b[PhaseID::W];
   }
   for (std::size_t f = 0; f < N_f; ++f)
   {
	   std::size_t c = ftoc(f);
	   const double VOLUME = fCells[f].vol*0.1781076;
	   mAccum[eqnID(c, PhaseID::O)] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::O] * mCells[c].b[PhaseID::O];
	   mAccum[eqnID(c, PhaseID::W)] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::W] * mCells[c].b[PhaseID::W];
   }
}

void 
DiscreteProblem::compute_flow( )
{
   compute_face_properties( );

   for (std::size_t c = 0; c < (2*N_c); ++c ) 
      mFlow[c] = 0.0;

   for (std::size_t f = 0; f < mMesh.size_faces(); ++f ) 
   {
      std::size_t c1 = mMesh.face_adj_cell_1( f );
      std::size_t c2 = mMesh.face_adj_cell_2( f );

      mTmpVars[PhaseID::W] = mFaces[f].T * mFaces[f].L[PhaseID::W] * mFaces[f].Pot[PhaseID::W];
      mTmpVars[PhaseID::O] = mFaces[f].T * mFaces[f].L[PhaseID::O] * mFaces[f].Pot[PhaseID::O];
      
      mFlow[ eqnID(c1,PhaseID::O) ] -= mTmpVars[ PhaseID::O ];
      mFlow[ eqnID(c1,PhaseID::W) ] -= mTmpVars[ PhaseID::W ];

      mFlow[ eqnID(c2,PhaseID::O) ] += mTmpVars[ PhaseID::O ];
      mFlow[ eqnID(c2,PhaseID::W) ] += mTmpVars[ PhaseID::W ];
   }
   for (std::size_t f = 0; f < N_fm; ++f)
   {
	   std::size_t c1 = mfc[f].ia;
	   std::size_t c2 = ftoc(mfc[f].ib);
	   mTmpVars[PhaseID::W] = mfc[f].trans * fFaces[f].L[PhaseID::W] * fFaces[f].Pot[PhaseID::W];
	   mTmpVars[PhaseID::O] = mfc[f].trans * fFaces[f].L[PhaseID::O] * fFaces[f].Pot[PhaseID::O];

	   mFlow[eqnID(c1, PhaseID::O)] -= mTmpVars[PhaseID::O];
	   mFlow[eqnID(c1, PhaseID::W)] -= mTmpVars[PhaseID::W];

	   mFlow[eqnID(c2, PhaseID::O)] += mTmpVars[PhaseID::O];
	   mFlow[eqnID(c2, PhaseID::W)] += mTmpVars[PhaseID::W];

   }
   for (std::size_t f = 0; f < N_ff; ++f)
   {
	   std::size_t c1 = ftoc(ffc[f].ia);
	   std::size_t c2 = ftoc(ffc[f].ib);
	   mTmpVars[PhaseID::W] = ffc[f].trans * fFaces[f + N_fm].L[PhaseID::W] * fFaces[f + N_fm].Pot[PhaseID::W];
	   mTmpVars[PhaseID::O] = ffc[f].trans * fFaces[f + N_fm].L[PhaseID::O] * fFaces[f + N_fm].Pot[PhaseID::O];

	   mFlow[eqnID(c1, PhaseID::O)] -= mTmpVars[PhaseID::O];
	   mFlow[eqnID(c1, PhaseID::W)] -= mTmpVars[PhaseID::W];

	   mFlow[eqnID(c2, PhaseID::O)] += mTmpVars[PhaseID::O];
	   mFlow[eqnID(c2, PhaseID::W)] += mTmpVars[PhaseID::W];
   }
   compute_hwells();
   //compute_wells( );
}

void 
DiscreteProblem::compute_face_properties( )
{
   for (std::size_t f = 0; f < mMesh.size_faces(); ++f ) 
   {
      std::size_t c1 = mMesh.face_adj_cell_1( f );
      std::size_t c2 = mMesh.face_adj_cell_2( f );
      double      dz = ( mMesh.cell_coord(c2) - mMesh.cell_coord(c1) ).p[2];

      for ( std::size_t ph=0; ph < 2; ++ph )
      {
	 mFaces[f].Pot[ph] = (mCells[c2].P[ph] - mCells[c1].P[ph]) + 
	    0.00694 * 0.5 * (mCells[c2].Rho[ph] + mCells[c1].Rho[ph]) * dz;
	    
	 mFaces[f].L[ph] = 0.5 * (mCells[c2].bmu[ph] + mCells[c1].bmu[ph]);
	    
	 if ( mFaces[f].Pot[ph].value() > 0.0 )
	    mFaces[f].L[ph] *= mCells[c2].Kr[ph];
	 else
	    mFaces[f].L[ph] *= mCells[c1].Kr[ph];
      }
   }

   for (std::size_t f = 0; f < N_fm; ++f)
   {
	   std::size_t c1 = mfc[f].ia;
	   std::size_t c2 = ftoc(mfc[f].ib);

	   for (std::size_t ph = 0; ph < 2; ++ph)
	   {
		   fFaces[f].Pot[ph] = (mCells[c2].P[ph] - mCells[c1].P[ph]);
		  // std::cout << fFaces[f].Pot[ph].value() << std::endl;

		   fFaces[f].L[ph] = 0.5 * (mCells[c2].bmu[ph] + mCells[c1].bmu[ph]);

		   if (fFaces[f].Pot[ph].value() > 0.0)
			   fFaces[f].L[ph] *= mCells[c2].Kr[ph];
		   else
			   fFaces[f].L[ph] *= mCells[c1].Kr[ph];
	   }
   }

   for (std::size_t f = 0; f < N_ff; ++f)
   {
	   std::size_t c1 = ftoc(ffc[f].ia);
	   std::size_t c2 = ftoc(ffc[f].ib);

	   for (std::size_t ph = 0; ph < 2; ++ph)
	   {
		   fFaces[f+N_fm].Pot[ph] = (mCells[c2].P[ph] - mCells[c1].P[ph]);

		   fFaces[f+N_fm].L[ph] = 0.5 * (mCells[c2].bmu[ph] + mCells[c1].bmu[ph]);

		   if (fFaces[f + N_fm].Pot[ph].value() > 0.0)
			   fFaces[f + N_fm].L[ph] *= mCells[c2].Kr[ph];
		   else
			   fFaces[f + N_fm].L[ph] *= mCells[c1].Kr[ph];
	   }
   }
}

void DiscreteProblem::compute_total_mobility(unsigned w, PhaseID I, adetl::ADscalar<> &T, adetl::ADscalar<> &Tp)
{
	T = 0;
	Tp = 0;
	for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
	{
		unsigned c = fracsys.h_well[w].gridindex[m];
		double WI = fracsys.h_well[w].WI[m];
		T += WI * mCells[c].Kr[I] * mCells[c].bmu[I];
		Tp += WI * mCells[c].Kr[I] * mCells[c].bmu[I] * mCells[c].P[0];
	}


}

void DiscreteProblem::compute_BHP()
{
	adetl::ADscalar<>  To,Top,Tw,Twp, T,Tp;
	for (std::size_t w = 0; w < fracsys.N_w; ++w)
	{
		if (fracsys.h_well[w].isproducer)
		{
			if (H_control[w][nsche] == 0)
			{
				Hwell_q[w].P = H_chedule[w][nsche];

			}
			else if (H_control[w][nsche] == 1)
			{
				compute_total_mobility(w, PhaseID::O, To, Top);
				compute_total_mobility(w, PhaseID::W, Tw, Twp);
				T = To + Tw;
				Tp = Top + Twp;
				Hwell_q[w].P = (Tp - H_chedule[w][nsche]) / T;
			}
			else if (H_control[w][nsche] == 2)
			{
				compute_total_mobility(w, PhaseID::O, To, Top);
				Hwell_q[w].P = (Top - H_chedule[w][nsche]) / To;

			}
			else if (H_control[w][nsche] == 3)
			{
				compute_total_mobility(w, PhaseID::W, Tw, Twp);
				Hwell_q[w].P = (Twp - H_chedule[w][nsche]) / Tw;

			}
			else
			{
				std::cout << "Invalid well control for producer. " << std::endl;
				system("pause");
			}
		}
		else
		{
			if (H_control[w][nsche] == 0)
			{
				Hwell_q[w].P = H_chedule[w][nsche];
			}
			else if (H_control[w][nsche] == 3)
			{
				compute_total_mobility(w, PhaseID::W, Tw, Twp);
				Hwell_q[w].P = (Twp - H_chedule[w][nsche]) / Tw;
			}
			else
			{
				std::cout << "Invalid well control for injector. " << std::endl;
				system("pause");
			}
		}
	}
}


void DiscreteProblem::compute_oil_water_ratio(std::vector<adetl::ADscalar<>> &BHP,unsigned w)
{

}


void
DiscreteProblem::compute_hwells()
{
	//////////////////////function use to compute both horizontal well and vertical well//////////////////////////////////
	if (schedule[nsche] == sum_t)
	{
		nsche++;
		std::cout <<"The "<<nsche<<"th schedule at "<<sum_t<<" day." << std::endl;
	}
	for (std::size_t w = 0; w < fracsys.N_w; ++w)
	{
		adetl::ADscalar<> qo = 0, qw = 0;
		if (fracsys.h_well[w].isproducer)
		{
			std::vector<adetl::ADscalar<>> BHP;
			//compute_oil_water_ratio(BHP,w);
			for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
			{
				unsigned c = fracsys.h_well[w].gridindex[m];
				if (mCells[c].P[PhaseID::O].value()>=Hwell_q[w].P.value())
				{
					mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
			
				}
				else
				{
					mTmpVars[0].value() = 0;
				}
				double WI = fracsys.h_well[w].WI[m];
				qo += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
				qw += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
				mFlow[eqnID(c, PhaseID::O)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
				mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
			}
			Hwell_q[w].qo = qo;
			Hwell_q[w].qw = qw;
			Hwell_q[w].qt = qo + qw;
			Hwell_q[w].WCT = qw / Hwell_q[w].qt;

		}
		else
		{
			if (H_control[w][nsche] == 0)
			{
				for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
				{
					unsigned c = fracsys.h_well[w].gridindex[m];
					if (mCells[c].P[PhaseID::O].value()<=Hwell_q[w].P.value())
					{
						mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
					}
					else
					{
						mTmpVars[0].value() = 0;
					}
					double WI = fracsys.h_well[w].WI[m];
					qw += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
					mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
				}
				Hwell_q[w].qw = qw;
			}
			else if (H_control[w][nsche] == 3)
			{
				for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
				{
					unsigned c = fracsys.h_well[w].gridindex[m];
					double WI = fracsys.h_well[w].WI[m];
					if (mCells[c].P[PhaseID::O].value() <= Hwell_q[w].P.value())
					{
						mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
					}
					else
					{
						mTmpVars[0].value() = 0;
					}
					mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
					qw += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
					mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
				}
				Hwell_q[w].qw = qw;
			}
			else
			{
				std::cout << "Invalid well control. " << std::endl;
			}
		}
	}

}

double
DiscreteProblem::compute_total_WI(Hor_well &h_well)
{
	double total=0;
	for (std::size_t m = 0; m < h_well.gridindex.size(); m++)
	{
		unsigned c = h_well.gridindex[m];
		total += h_well.WI[m]*mCells[c].Kr[PhaseID::W].value() * mCells[c].bmu[PhaseID::W].value();
	}

	return total;
}


void DiscreteProblem::get_inj_index(std::vector<unsigned> &iidx)
{
	for (std::size_t w = 0; w < fracsys.N_w; ++w)
	{
		if (!fracsys.h_well[w].isproducer)
		{
			iidx.push_back(fracsys.h_well[w].gridindex[0]);
		}
	}
}

double
DiscreteProblem::safeguard_MAC( double upd )
{
   if ( std::abs( upd ) > 0.2 )
      return upd/std::abs(upd) * 0.2;
   else
      return upd;

}




std::ostream & operator << ( std::ostream & ostr, 
			     const DiscreteProblem::ConvergenceInfo & _out )
{
   ostr << _out.NormSat[0]  << "\t" 
	<< _out.NormSat[1]  << "\t" 
	<< _out.MatBal[0]  << "\t" 
	<< _out.MatBal[1]  << std::endl; 

   return ostr;
}


/////////////////////////////////////////////////////Code for Adjoint//////////////////////////////////////////////////////////////////////////

void DiscreteProblem::assign_wellrate(std::vector<wellrate> &temp_v)
{
	temp_v.resize(nw);
	for (unsigned i = 0; i < nw;i++)
	{
		temp_v[i].qo = Hwell_q[i].qo;
		temp_v[i].qw = Hwell_q[i].qw;
		temp_v[i].P = Hwell_q[i].P;
		temp_v[i].qt = Hwell_q[i].qt;
		temp_v[i].WCT = Hwell_q[i].WCT;
	}
}




void DiscreteProblem::extract_obj_der(unsigned nsche)
{
	unsigned ny = max_num_eqns();
	std::vector<double> obj_y(ny);
	adetl::ADscalar<> x;
	//assign_wellrate(temp_y);
	for (unsigned i = 0; i < nw;i++)
	{
		if (H_constrain[nsche][i].BHP_i)
		{			
			H_constrain[nsche][i].BHP = -Hwell_q[i].P;
		}
		if (H_constrain[nsche][i].qt_i)
		{
			H_constrain[nsche][i].qt = -Hwell_q[i].qt;
		}
		if (H_constrain[nsche][i].qo_i)
		{
			H_constrain[nsche][i].qo = -Hwell_q[i].qo;
			well_rate[i].qo = Hwell_q[i].qo.value();
			//std::cout << Hwell_q[i].qo.value() << std::endl;
		}
		if (H_constrain[nsche][i].qw_i)
		{
			H_constrain[nsche][i].qw = -Hwell_q[i].qw;
			well_rate[i].qw = Hwell_q[i].qw.value();
		}
		if (H_constrain[nsche][i].wct_i)
		{
			H_constrain[nsche][i].wct = -Hwell_q[i].WCT;
		}
	}
}



void DiscreteProblem::extract_F_der(DiscreteProblem::StateVector &oldstate, DiscreteProblem::StateVector &state, double DT, DiscreteProblem::Sparse_matrix &M)
{
	std::vector<double> r(N_c);
	CSR f_y;
	bind_to_old_state(oldstate);
	discretize_m(state, DT);
	extract_R_Jm(r, f_y, 0);
	//std::cout << f_y.value().size() << std::endl;
	CSR_to_Sparse(f_y, M);

}

void
DiscreteProblem::extract_acc(DiscreteProblem::Sparse_matrix &M)
{
	std::size_t offset = 0;
	std::vector<double> r(N_tc,0);
	CSR m;
	mAccum.extract_CSR_non_sysmetric(r, m.rowptr(), m.colind(), m.value(),N_tc);
	for (std::size_t i = 0; i< m.rowptr().size(); ++i) m.rowptr()[i] += offset;
	for (std::size_t i = 0; i< m.colind().size(); ++i) m.colind()[i] += offset;
	m.check_size();
	CSR_to_Sparse(m, M);
}


void
DiscreteProblem::extract_R_Jm(std::vector<double> &r, CSR &m, std::size_t offset)
{
	mResidual.extract_CSR_non_sysmetric(r, m.rowptr(), m.colind(), m.value(), nwc);
	for (std::size_t i = 0; i< m.rowptr().size(); ++i) m.rowptr()[i] += offset;
	for (std::size_t i = 0; i< m.colind().size(); ++i) m.colind()[i] += offset;
	m.check_size();
	//if (r.size() != mResidual.size()) std::cout << "BUG IN JACOBIAN\t zero row found" << r.size() << "!=" << mResidual.size() << std::endl;
}

void DiscreteProblem::CSR_to_Sparse(CSR & x, DiscreteProblem::Sparse_matrix &M)
{
	unsigned n = x.value().size();
	unsigned  old = 0;
	M.row.clear();
	M.col = x.colind();
	M.val = x.value();
	for (unsigned i = 0; i < x.rowptr().size(); i++)
	{
		if (x.rowptr()[i] != 0)
		{
			for (unsigned j = 0; j < (x.rowptr()[i] - old); j++)
			{
				M.row.push_back(i);
			}
			old = x.rowptr()[i];
		}
	}
}


////////////////////////////////////////////////////Code for Adjoint//////////////////////////////////////////////////////////////////////////