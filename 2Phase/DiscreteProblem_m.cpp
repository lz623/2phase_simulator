#include "DiscreteProblem_m.hpp"
#include "MeshBuilder.hpp"
#include "Point.hpp"
#include <fstream>

DiscreteProblem_m::DiscreteProblem_m(std::size_t NX, std::size_t NY, std::size_t NZ,
	double LX, double LY, double LZ,
	const std::vector<double> &vKX, const std::vector<double> &vKY, const std::vector<double> &vKZ,
	const std::vector<double> &Phi_ref) :
	mMesh(MeshBuilder<UniformCartesian>(NX, NY, NZ, LX, LY, LZ)),
	mPropCalc(),
	fracsys(NX, NY, NZ, LX, LY, LZ, vKX, Phi_ref),
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
	mResidualm.resize(2 * N_c, 0);
	initialize_transmissibility(vKX, vKY, vKZ);
	setup_wells(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ);
};


void DiscreteProblem_m::transfer_data(EDFM& model)
{
	unsigned n = model.inter_nn.coor.size();
	connection tmp;
	N_fm = model.n_frac.size();
	N_f = N_fm;
	N_ff = model.inter_nn.i_index.size();
	N_m = mMesh.size_cells();
	N_c = N_fm + N_m;
	fCells.resize(N_fm);
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
	for (unsigned i = 0; i < fracsys.h_well.size(); i++)
	{
		fracsys.h_well[i].ifbreakthrough = false;
	}

}



void DiscreteProblem_m::extract_inform(const std::vector<std::vector<unsigned>> &x, std::vector<std::vector<unsigned>> &inf)
{
	unsigned key;
	inf.resize(nw);
	for (unsigned i = 0; i < nw; i++)
	{
		inf[i].resize(5);
		std::fill(inf[i].begin(), inf[i].end(), 0);
		for (unsigned j = 0; j <x[i].size(); j++)
		{
			key = x[i][j];
			inf[i][key] = j;
		}
	}
}

void
read_from_file_unsigned(std::string filename, dlib::matrix<unsigned> &v)
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

void DiscreteProblem_m::loadwell_data(const std::vector<double> &time, const dlib::matrix<double> &well_sche)
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
	load_hschedule(H_constrain,cs);
}

void DiscreteProblem_m::check_schedule_size(unsigned nw, unsigned nct, unsigned nwc)
{
	if (nw*nct != nwc)
	{
		std::cout << "The schedule's size is not match" << std::endl;
		system("pause");
	}

}

void DiscreteProblem_m::load_hschedule(std::vector<std::vector<unsigned>> &hwells,  dlib::matrix<unsigned> &s)
{
	hwells.resize(nw);
	for (unsigned i = 0; i < nw; i++)
	{
		hwells[i].resize(nct);
		for (unsigned j = 0; j <nct; j++)
		{
			hwells[i][j] = s(i*nct + j);
		}
	}
	//for (unsigned i = Np; i < nw; i++)
	//{
	//	hwells[i].resize(nct);
	//	for (unsigned j = 0; j < nct; j++)
	//	{
	//		hwells[i][j] = s(i*nct + j);
	//	}
	//}
}

void
DiscreteProblem_m::setup_wells(std::size_t NX, std::size_t NY, std::size_t NZ,
double LX, double LY, double LZ,
const std::vector<double> &vKX,
const std::vector<double> &vKY,
const std::vector<double> &vKZ)
{
	/*nw = read_well();*/
	const double DX = LX / static_cast<double>(NX);
	const double DY = LY / static_cast<double>(NY);
	const double DZ = LZ / static_cast<double>(NZ);
	N_hp = fracsys.np;
	N_hj = fracsys.nj;
	H_chedule.resize(fracsys.N_w);
	nw = fracsys.N_w;

	//for (std::size_t i = 0; i < nw; i++)
	//{
	//	const double Kh = DZ * sqrt(vKY[mWells[i].loc] * vKX[mWells[i].loc]);
	//	const double r1 = vKY[mWells[i].loc] / vKX[mWells[i].loc];
	//	const double r2 = vKX[mWells[i].loc] / vKY[mWells[i].loc];
	//	const double ro = 0.28 * std::sqrt(std::sqrt(r1)*DX*DX + std::sqrt(r2)*DY*DY) / (std::pow(r1, 0.25) + std::pow(r2, 0.25));
	//	mWells[i].WI = 0.00708*Kh / log(ro / mWells[i].rw);
	//	fracsys.find_v_intersect(mWells[i].loc, mWells[i].f_index, mWells[i].f_trans,8);
	//}

	for (std::size_t i = 0; i <fracsys.N_w; i++)
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

	btt.resize(N_hp);
	Hwell_q.resize(fracsys.N_w);
	//Well_Q.resize(nw);

}

void
DiscreteProblem_m::initialize_state(DiscreteProblem_m::StateVector &uold, DiscreteProblem_m::StateVector &unew,
DiscreteProblem_m::StateVector &aold, DiscreteProblem_m::StateVector &anew)
{
	for (unsigned c = 0; c< N_c;c++)
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
DiscreteProblem_m::bind_to_old_state(const DiscreteProblem_m::StateVector &old_state)
{
	compute_cell_properties(old_state);
	compute_accumulation();
	for (std::size_t c = 0; c<N_c; ++c)
	{
		std::size_t eqn1 = eqnID(c, PhaseID::O);
		std::size_t eqn2 = eqnID(c, PhaseID::W);
		mAccum_old[eqn1] = mAccum[eqn1].value();
		mAccum_old[eqn2] = mAccum[eqn2].value();
	}
}


bool
DiscreteProblem_m::discretize(const DiscreteProblem_m::StateVector &state, double DT)
{
	bool is_badvalue = false;
	mDT = DT;
	compute_cell_properties(state);
	compute_accumulation();
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

	return is_badvalue;
}




bool
DiscreteProblem_m::is_converged(ConvergenceInfo & nrm)
{
	bool is_MATBAL_converged = true;
	bool is_NRMSAT_converged = true;

	double tot_PV = 0.0;
	for (std::size_t c = 0; c<mMesh.size_cells(); ++c)
		tot_PV += mCells[c].phi.value() * mMesh.cell_measure(c)*0.1781076;
	for (std::size_t f = 0; f<N_f; ++f)
	{
		std::size_t c = ftoc(f);
		tot_PV += mCells[c].phi.value() * fCells[f].vol*0.1781076;
	}


	for (std::size_t phs = 0; phs < 2; ++phs)
	{
		double sum_R = 0.0;
		double avg_B = 0.0;
		double max_R_PV = 0.0;
		for (std::size_t c = 0; c<mMesh.size_cells(); ++c)
		{
			sum_R += mResidual[eqnID(c, phs)].value();
			avg_B += 1.0 / mCells[c].b[phs].value();
			double R_PV = std::abs(mResidual[eqnID(c, phs)].value() / (mCells[c].phi.value() * mMesh.cell_measure(c)*0.1781076));
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
		nrm.MatBal[phs] = std::abs(avg_B * sum_R / tot_PV);
		nrm.NormSat[phs] = avg_B * max_R_PV;
		if (nrm.MatBal[phs] > 1.0e-7)  is_MATBAL_converged = false;
		if (nrm.NormSat[phs] > 0.001)  is_NRMSAT_converged = false;
	}
	return (is_MATBAL_converged && is_NRMSAT_converged);
}

void
DiscreteProblem_m::initialize_transmissibility(const std::vector<double> & KX,
const std::vector<double> & KY,
const std::vector<double> & KZ)
{
	Point KL, KR;
	for (std::size_t f = 0; f < mMesh.size_faces(); ++f)
	{
		const std::size_t c1 = mMesh.face_adj_cell_1(f);
		const std::size_t c2 = mMesh.face_adj_cell_2(f);
		KL.p[0] = KX[c1];
		KL.p[1] = KY[c1];
		KL.p[2] = KZ[c1];
		KR.p[0] = KX[c2];
		KR.p[1] = KY[c2];
		KR.p[2] = KZ[c2];
		const Point nrml = mMesh.unit_normal(f);
		const double kl = dot_product(nrml, KL);
		const double kr = dot_product(nrml, KR);
		const double DX = norm(mMesh.cell_coord(c1) - mMesh.cell_coord(c2));
		const double A = mMesh.face_measure(f);
		const double beta = DX*0.5 / A*(1.0 / kl + 1.0 / kr);
		mFaces[f].T = 0.00112712 / beta;
	}
};

void
DiscreteProblem_m::compute_cell_properties(const DiscreteProblem_m::StateVector &state)
{
	for (std::size_t c = 0; c < mMesh.size_cells(); ++c)
	{
		mPropCalc.calculate(mPhi_ref[c], state[c], mCells[c]);
	}
	for (std::size_t f = 0; f < N_f; ++f)
	{
		std::size_t c = ftoc(f);
		mPropCalc.calculate_f(fCells[f].phi, state[c], mCells[c]);
	}
}

void
DiscreteProblem_m::compute_accumulation()
{
	for (std::size_t c = 0; c < mMesh.size_cells(); ++c)
	{
		const double VOLUME = mMesh.cell_measure(c)*0.1781076;
		mAccum[eqnID(c, PhaseID::O)] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::O] * mCells[c].b[PhaseID::O];
		mAccum[eqnID(c, PhaseID::W)] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::W] * mCells[c].b[PhaseID::W];
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
DiscreteProblem_m::compute_flow()
{
	compute_face_properties();

	for (std::size_t c = 0; c < (2 * N_c); ++c)
		mFlow[c] = 0.0;

	for (std::size_t f = 0; f < mMesh.size_faces(); ++f)
	{
		std::size_t c1 = mMesh.face_adj_cell_1(f);
		std::size_t c2 = mMesh.face_adj_cell_2(f);

		mTmpVars[PhaseID::W] = mFaces[f].T * mFaces[f].L[PhaseID::W] * mFaces[f].Pot[PhaseID::W];
		mTmpVars[PhaseID::O] = mFaces[f].T * mFaces[f].L[PhaseID::O] * mFaces[f].Pot[PhaseID::O];

		mFlow[eqnID(c1, PhaseID::O)] -= mTmpVars[PhaseID::O];
		mFlow[eqnID(c1, PhaseID::W)] -= mTmpVars[PhaseID::W];

		mFlow[eqnID(c2, PhaseID::O)] += mTmpVars[PhaseID::O];
		mFlow[eqnID(c2, PhaseID::W)] += mTmpVars[PhaseID::W];
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
DiscreteProblem_m::compute_face_properties()
{
	for (std::size_t f = 0; f < mMesh.size_faces(); ++f)
	{
		std::size_t c1 = mMesh.face_adj_cell_1(f);
		std::size_t c2 = mMesh.face_adj_cell_2(f);
		double      dz = (mMesh.cell_coord(c2) - mMesh.cell_coord(c1)).p[2];

		for (std::size_t ph = 0; ph < 2; ++ph)
		{
			mFaces[f].Pot[ph] = (mCells[c2].P[ph] - mCells[c1].P[ph]) +
				0.00694 * 0.5 * (mCells[c2].Rho[ph] + mCells[c1].Rho[ph]) * dz;

			mFaces[f].L[ph] = 0.5 * (mCells[c2].bmu[ph] + mCells[c1].bmu[ph]);

			if (mFaces[f].Pot[ph].value() > 0.0)
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
			fFaces[f + N_fm].Pot[ph] = (mCells[c2].P[ph] - mCells[c1].P[ph]);

			fFaces[f + N_fm].L[ph] = 0.5 * (mCells[c2].bmu[ph] + mCells[c1].bmu[ph]);

			if (fFaces[f + N_fm].Pot[ph].value() > 0.0)
				fFaces[f + N_fm].L[ph] *= mCells[c2].Kr[ph];
			else
				fFaces[f + N_fm].L[ph] *= mCells[c1].Kr[ph];
		}
	}
}

void
DiscreteProblem_m::compute_hwells()
{
	//////////////////////function use to compute both horizontal well and vertical well//////////////////////////////////
	if (schedule[nsche] == sum_t)
	{
		nsche++;
		std::cout << "The " << nsche << "th schedule at " << sum_t << " day." << std::endl;
	}

	for (std::size_t w = 0; w < fracsys.N_w; ++w)
	{
		adetl::ADscalar<> qo = 0, qw = 0, qinj = 0, To, T;
		mTmpVars[0] = 0;
		if (fracsys.h_well[w].isproducer)
		{
			bool icv_is_open(false);
			if (H_control[w][nsche] == 0)
			{
				for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
				{
					unsigned c = fracsys.h_well[w].gridindex[m];
					double WI = fracsys.h_well[w].WI[m];

					if (mCells[c].P[PhaseID::O].value()>H_chedule[w][nsche].value())
					{
						mTmpVars[0] = mCells[c].P[PhaseID::O] - H_chedule[w][nsche];
						icv_is_open = true;
					}
					else
					{
						mTmpVars[0].value() = 0;
					}

					qo += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
					qw += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];

					mFlow[eqnID(c, PhaseID::O)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
					mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
				}
				/*	std::cout << qo.value() << std::endl;*/
				Hwell_q[w].qo = qo;
				Hwell_q[w].qw = qw;
				Hwell_q[w].qt = qo + qw;
				Hwell_q[w].WCT = qw / Hwell_q[w].qt;
			}
			else if (H_control[w][nsche] == 1)
			{
				adetl::ADscalar<> BHP;
				To = 0, T = 0;
				for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
				{
					unsigned c = fracsys.h_well[w].gridindex[m];
					double WI = fracsys.h_well[w].WI[m];
					T += WI * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
					T += WI * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
					To += WI * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O] * mCells[c].P[PhaseID::O];
					To += WI * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W] * mCells[c].P[PhaseID::O];
				}
				Hwell_q[w].P = (To - H_chedule[w][nsche]) / T;
				for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
				{
					unsigned c = fracsys.h_well[w].gridindex[m];
					double WI = fracsys.h_well[w].WI[m];

					if (mCells[c].P[PhaseID::O].value()>Hwell_q[w].P.value())
					{
						mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
						icv_is_open = true;
					}
					else
					{
						mTmpVars[0].value() = 0;
					}

					qo += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
					qw += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];

					mFlow[eqnID(c, PhaseID::O)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
					mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];

				}
				Hwell_q[w].qo = qo;
				Hwell_q[w].qw = qw;
				Hwell_q[w].qt = qo + qw;
				Hwell_q[w].WCT = qw / Hwell_q[w].qt;
				//std::cout << "Don't know how to do it " << std::endl;
			}
			else if (H_control[w][nsche] == 2)
			{
				std::cout << "Invalid well control for specify qo " << std::endl;
			}
			else if (H_control[w][nsche] == 3)
			{
				std::cout << "Invalid well control for specify qw " << std::endl;
			}
			else
			{
				std::cout << "Invalid well control for specify wct " << std::endl;
			}
		}
		else
		{
			if (H_control[w][nsche] == 0)
			{
				for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
				{
					unsigned c = fracsys.h_well[w].gridindex[m];
					double WI = fracsys.h_well[w].WI[m];

					if (mCells[c].P[PhaseID::O].value()<H_chedule[w][nsche].value())
					{
						mTmpVars[0] = mCells[c].P[PhaseID::O] - H_chedule[w][nsche];
					}
					else
					{
						mTmpVars[0].value() = 0;
					}
					qw += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];

					mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
				}
				Hwell_q[w].qw = qw;
				Hwell_q[w].P = H_control[w][nsche];
			}
			else if (H_control[w][nsche] == 3)
			{
				To = 0, T = 0;
				/*	double WI_t = compute_total_WI(fracsys.h_well[w]);*/
				for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
				{
					unsigned c = fracsys.h_well[w].gridindex[m];
					To += fracsys.h_well[w].WI[m] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W] * mCells[c].P[PhaseID::O];
					T += fracsys.h_well[w].WI[m] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
					//qinj += WI / WI_t*H_chedule[w][nsche];
				}
				Hwell_q[w].P = (H_chedule[w][nsche] + To) / T;
				for (std::size_t m = 0; m < fracsys.h_well[w].gridindex.size(); m++)
				{
					unsigned c = fracsys.h_well[w].gridindex[m];
					double WI = fracsys.h_well[w].WI[m];
					mTmpVars[0] = mCells[c].P[PhaseID::O] - Hwell_q[w].P;
					qw += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
					mFlow[eqnID(c, PhaseID::W)] += WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
					//////////////////////////////////////////////////////////////////////////////////////////////////////
					//std::vector<double> x(72, 0);
					//double test;
					//H_chedule[i][j].extract(test, x.data(), 60);
					//std::cout << " start to output" << i*time.size() + j << std::endl;
					//for (unsigned ss = 0; ss < x.size(); ss++)
					//{

					//	std::cout << ss << "   " << x[ss] << std::endl;
					//}
				}
				Hwell_q[w].qw = H_control[w][nsche];

			}
			else
			{
				std::cout << "Invalid well control. " << std::endl;
			}
		}
	}
}


//void DiscreteProblem_m::get_inj_index(std::vector<unsigned> &iidx)
//{
//	for (std::size_t w = 0; w < fracsys.N_w; ++w)
//	{
//		if (!fracsys.h_well[w].isproducer)
//		{
//			iidx.push_back(fracsys.h_well[w].gridindex[0]);
//			//std::cout << fracsys.h_well[w].gridindex[0] << std::endl;
//		}
//	}
//}

double
DiscreteProblem_m::safeguard_MAC(double upd)
{
	if (std::abs(upd) > 0.2)
		return upd / std::abs(upd) * 0.2;
	else
		return upd;

}

void DiscreteProblem_m::CSR_to_Sparse(CSR & x, DiscreteProblem_m::Sparse_matrix &M)
{
	unsigned n=x.value().size();
	unsigned  old=0;
	M.row.clear();
	M.col=x.colind();
	M.val=x.value();
	for (unsigned i = 0; i < x.rowptr().size();i++)
	{
		if (x.rowptr()[i]!=0)
		{
			for (unsigned j = 0; j < (x.rowptr()[i] - old);j++)
			{
				M.row.push_back(i);
			}
			old = x.rowptr()[i];
		}
	}
}


void DiscreteProblem_m::extract_F_der(DiscreteProblem_m::StateVector &oldstate, DiscreteProblem_m::StateVector &state, double DT, DiscreteProblem_m::Sparse_matrix &M)
{
	std::vector<double> r(N_c);
	CSR f_y;
	bind_to_old_state(oldstate);
	discretize(state, DT);
	extract_R_J(r, f_y, 0);
	//std::cout << f_y.value().size() << std::endl;
	CSR_to_Sparse(f_y, M);

}

void
DiscreteProblem_m::extract_R_J(std::vector<double> &r, CSR &m, std::size_t offset)
{
	mResidual.extract_CSR_non_sysmetric(r, m.rowptr(), m.colind(), m.value(),nwc);
	for (std::size_t i = 0; i< m.rowptr().size(); ++i) m.rowptr()[i] += offset;
	for (std::size_t i = 0; i< m.colind().size(); ++i) m.colind()[i] += offset;
	m.check_size();
	//if (r.size() != mResidual.size()) std::cout << "BUG IN JACOBIAN\t zero row found" << r.size() << "!=" << mResidual.size() << std::endl;
}

void DiscreteProblem_m::extract_obj_der(unsigned nsche,
	std::vector<adetl::ADscalar<>> &objy)
{
	unsigned ny = max_num_eqns();
	std::vector<double> obj_y(ny);
	double tmp;
	adetl::ADscalar<> x;
	//assign_wellrate(temp_y);
	for (unsigned i = 0; i < nw; i++)
	{
		if (H_constrain[i][nsche] == 0)
		{
			objy[i*nct + nsche] = Hwell_q[i].P;
			//objv[i*nct + nsche] = temp_v[i].P;
			tmp = Hwell_q[i].P.value();
			//x = -temp_y[i].P;
			//std::cout << "extract at " << nsche << std::endl;
			//x.extract(tmp, obj_y.data(), ny);
			//std::cout << "finish extract " << nsche << std::endl;
			//temp_v[i].P.extract(tmp, obj_v.data(), ny);
		}
		else if (H_constrain[i][nsche] == 3)
		{
			objy[i*nct + nsche] = Hwell_q[i].qw;
			//objv[i*nct + nsche] = temp_v[i].qw;
			tmp = Hwell_q[i].qw.value();
			//x = -temp_y[i].qw;
			//x.extract(tmp, obj_y.data(), ny);
			//temp_v[i].qw.extract(tmp, obj_v.data(), ny);
		}
		else if (H_constrain[i][nsche] == 2)
		{
			objy[i*nct + nsche] = Hwell_q[i].qo;
			//objv[i*nct + nsche] = temp_v[i].qo;
			tmp = Hwell_q[i].qo.value();
			//x = -temp_y[i].qo;
			//x.extract(tmp, obj_y.data(), ny);
			//temp_v[i].qo.extract(tmp, obj_v.data(), ny);
		}
		else if (H_constrain[i][nsche] == 1)
		{
			objy[i*nct + nsche] = Hwell_q[i].qt;
			//objv[i*nct + nsche] = temp_v[i].qt;
			tmp = Hwell_q[i].qt.value();
			//x = -temp_y[i].qt;
			//x.extract(tmp1, obj_y.data(), ny);
			//temp_v[i].qt.extract(tmp, obj_v.data(), ny);
			//std::cout << tmp <<"   "<<std::endl;
		}
		else if (H_constrain[i][nsche] == 4)
		{
			objy[i*nct + nsche] = Hwell_q[i].WCT;
			//objv[i*nct + nsche] = temp_v[i].WCT;
			tmp = Hwell_q[i].WCT.value();
			//x = -temp_y[i].WCT;
			//x.extract(tmp, obj_y.data(), ny);
			//temp_v[i].WCT.extract(tmp, obj_v.data(), ny);
		}
		else
		{
			std::cout << "Invalid input well constrain" << std::endl;
			system("pause");
		}

	}
}
