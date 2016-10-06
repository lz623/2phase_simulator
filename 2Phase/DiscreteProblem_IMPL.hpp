#pragma once

template< typename V, typename M >
void 
DiscreteProblem::extract_R_J( V &r, M &m, std::size_t offset )
{
   mResidual.extract_CSR( r, m.rowptr(), m.colind(), m.value() );
   for ( std::size_t i=0 ; i< m.rowptr().size() ; ++ i ) m.rowptr()[i] += offset;
   for ( std::size_t i=0 ; i< m.colind().size() ; ++ i ) m.colind()[i] += offset;
   m.check_size();
   if ( r.size() != mResidual.size() ) std::cout << "BUG IN JACOBIAN\t zero row found" << r.size() << "!=" << mResidual.size() << std::endl;
}

template< typename V, typename M >
void
DiscreteProblem::extract_R_J_nc(V &r, M &m, std::size_t offset)
{
	mResidual.extract_CSR_non_sysmetric(r, m.rowptr(), m.colind(), m.value(),nw*nct);
	for (std::size_t i = 0; i< m.rowptr().size(); ++i) m.rowptr()[i] += offset;
	for (std::size_t i = 0; i< m.colind().size(); ++i) m.colind()[i] += offset;
	m.check_size();
	if (r.size() != mResidual.size()) std::cout << "BUG IN JACOBIAN\t zero row found" << r.size() << "!=" << mResidual.size() << std::endl;
}



template< typename V, typename M >
void
DiscreteProblem::extract_acc(V &r, M &m, std::size_t offset)
{
	mAccum.extract_CSR_non_sysmetric(r, m.rowptr(), m.colind(), m.value());
	for (std::size_t i = 0; i< m.rowptr().size(); ++i) m.rowptr()[i] += offset;
	for (std::size_t i = 0; i< m.colind().size(); ++i) m.colind()[i] += offset;
	m.check_size();
	if (r.size() != mResidual.size()) std::cout << "BUG IN OLD JACOBIAN\t zero row found" << r.size() << "!=" << mAccum.size() << std::endl;

}


template< typename V >
void 
DiscreteProblem::extract_dR_dDT( V &r )
{
   r.resize( mFlow.size() );
   for ( std::size_t i=0 ; i< mFlow.size() ; ++ i ) r[i] = mFlow[i].value();
}

template< typename R >
void 
DiscreteProblem::update_state( StateVector &state, const R & update, bool do_safeguard )
{
   for ( std::size_t c = 0; c < N_c; ++c )
   {
      state[c].Po  += update[ eqnID(c,PhaseID::O) ];
      if (do_safeguard)
         state[c].Sw  += safeguard_MAC( update[ eqnID(c,PhaseID::W) ] );
      else
         state[c].Sw  += update[ eqnID(c,PhaseID::W) ];
   }
   for (std::size_t w = 0; w < nw;++w)
   {
	   std::size_t c = wtoc(w);
	   Hwell_q[w].P += update[c];
   }
}
   
template<typename x>
void DiscreteProblem::load_hschedule(std::vector<std::vector<x>> &hwells,  dlib::matrix<x> &s)
{
	hwells.resize(nw);
	for (unsigned i = 0; i < Np; i++)
	{
		hwells[i].resize(nct);
		for (unsigned j = 0; j <nct; j++)
		{
			hwells[i][j] = s(i*nct + j);
		}
	}
	for (unsigned i = Np; i < nw; i++)
	{
		hwells[i].resize(nct);
		for (unsigned j = 0; j < nct; j++)
		{
			hwells[i][j] = s(i*nct + j);
		}
	}
}