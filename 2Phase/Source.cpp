#include <dlib/matrix.h>
#include <fstream>
#include "simulator.cpp"


typedef struct
{
	bool        is_producer;
	std::size_t loc;
	bool control_type;
	double      WI;
	double      Pbh;
	double      rw;
	double      QINJ[3];
	double      Rso;
}                                            Well;
std::vector<Well>   mWells;
std::size_t Np, Nj, nw;

void
read_from_file(const char * filename, dlib::matrix<double> &v)
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


std::size_t read_well()
{
	std::ifstream strm("input/well_inform.dat");

	Np = 0;
	Nj = 0;
	strm >> nw;
	mWells.resize(nw);
	for (unsigned i = 0; i < nw; i++)
	{
		strm >> mWells[i].is_producer;
		(mWells[i].is_producer) ? Np++ : Nj++;
		strm >> mWells[i].loc;
		strm >> mWells[i].rw;
		strm >> mWells[i].control_type;
		if (mWells[i].control_type)
		{
			strm >> mWells[i].Pbh;
		}
		else
		{
			strm >> mWells[i].QINJ[0];
			strm >> mWells[i].QINJ[1];
			strm >> mWells[i].QINJ[2];
		}
		mWells[i].Rso = 0;
	}
	return nw;
}



int main()
{
	dlib::matrix<double> m, d;
	read_from_file("input/m_pri.dat", m);
	read_well();
	simulator orig(m, d);

	return -1;
}