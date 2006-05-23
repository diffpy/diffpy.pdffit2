// Up to date with 1.3.10 Fortran version

// in Fortran this was fourier.f

#include "pdffit.h"
#include <iomanip>
#include <sstream>

extern "C"
{
    void anoma_(char*,char*,int*);
    void fsc_(char*,int*);

    //common /anomal/ delfr,delfi
    extern struct {
            float delfr, delfil;
    } anomal_;

    //common /neu/    fneu
    extern struct {
        float fneu;    // neutron scattering length
    } neu_;

    //common /scat/   fa,fb,fc
    extern struct {
        // float fa(4,8),fb(4,8),fc(8);  FORTRAN
        float fa[8][4], fb[8][4], fc[8];
    } scat_;
}

/**************************************************************************
        This routine dlinks between the PDFFIT package and the routines
        that generate the atomic form factors. These routines have been
        adapted from the program LAZY.
**************************************************************************/
void Phase::dlink(matrix<double> &scat, bool lxray)
{
    for (int i=0; i<nscat; i++)
        get_scat(i, scat[i], lxray);
}

// get the scattering lengths for element <i>
void Phase::get_scat(int i, vector<double> &scat, bool lxray)
{
    char symwl[4];
    char element[4];
    fill(symwl, symwl+4, ' ');
    fill(element, element+4, ' ');
    const char* p = at_lis[i].c_str();
    size_t element_length = min(size_t(4), at_lis[i].length());
    copy(p, p+element_length, element);

    if(!lxray)
    {
        // check if scattering length has been set by the user
        if (Nscatlen[i].empty())
        {
            // neutron scattering
	    int kodlp = 2;
            element[2] = element[3] = ' ';
            anoma_(symwl, element, &kodlp);  // FORTRAN call
            if (element[0])
                scat[0] = neu_.fneu;
        }
        else
            scat[0] = Nscatlen[i][0];

        for(int j=1; j<9; j++)
            scat[j] = 0.0;
    }
    else
    {
        // check if scattering length has been set by the user
        if (Xscatfac[i].empty())
        {
            // x-ray scattering
            int one = 1;
            fsc_(element, &one);   // FORTRAN call
            if (element[0])
            {
                scat[0] = scat_.fc[0];
                for (int j=0; j<4; j++)
                {
                    scat[2*j+1] = scat_.fa[0][j];
                    scat[2*j+2] = scat_.fb[0][j];
                }
            }
        }
        else
            for (int j=0; j<9; j++)
                scat[j] = Xscatfac[i][j];
    }
}


/*****************************************************************
* calculates scattering factor for x-ray or neutron radiation
******************************************************************/
double Phase::scatteringFactor(vector<double>& scat, bool lxray)
{
    double f;
    f = scat[0];
    if (lxray)
    {
	// for x-rays [Acta Cryst. (1968) A24, p321
	//     f = sum(a_i * exp( - b_i * (sin(th)/lambda)^2 ) ) + c
	// However we assume Q = 0, thus f sums to just about atomic number
	for (size_t aidx = 1; aidx < 9; aidx += 2)
	{
	    f += scat[aidx];
	}
    }
    return f;
}

string Phase::show_scat(Sctp type)
{
    stringstream sout;
    for (int i=0; i<nscat; i++)
        sout << show_scat(type, i);

    return sout.str();
}

// local funtion (takes itype between 0 and nscat-1
string Phase::show_scat(Sctp type, int i)
{
    vector<double> scat(9);

    stringstream sout;
    get_scat(i, scat, (type==N) ? false : true);

    if(type == N)
    {
        sout << "Neutron scattering length for " << at_lis[i] << " :  ";
        sout << scat[0] << endl;
    }
    else if(type == X)
    {
        sout << "X-ray scattering factors for " << at_lis[i] << " :  \n";
        sout << "   a : ";
        for (int j=1; j<9; j+=2)
            sout << scat[j] << " ";
        sout << endl;

        sout << "   b : ";
        for (int j=2; j<9; j+=2)
            sout << scat[j] << " ";
        sout << endl;

        sout << "   c : " << scat[0] << endl << endl;
    }
    cout << sout.str();
    return sout.str();

}

void Phase::reset_scat(Sctp type, int itype)
{
    if ((itype<1)||(itype>nscat))
    {
        //warning("Incorrect atom type number");
        throw ValueError("Incorrect atom type number");
        return;
    }
    itype--;
    if (type == N)
        Nscatlen[itype].clear();
    else if (type == X)
        Xscatfac[itype].clear();

    show_scat(type, itype);
}

void Phase::set_scat(Sctp type, int itype, double len)
{
    if (type != N)
    {
        //warning("Incorrect arguments");
        throw ValueError("Incorrect arguments");
        return;
    }
    if ((itype<1)||(itype>nscat))
    {
        //warning("Incorrect atom type number");
        throw ValueError("Incorrect atom type number");
        return;
    }
    itype--;
    Nscatlen[itype].push_back(len);
    show_scat(type, itype);
}

void Phase::set_scat(Sctp type, int itype, double a1, double b1, double a2, double b2,
    double a3, double b3, double a4, double b4, double c)
{
    if (type != X)
    {
        //warning("Incorrect arguments");
        throw ValueError("Incorrect arguments");
        return;
    }
    if ((itype<1)||(itype>nscat))
    {
        //warning("Incorrect atom type number");
        throw ValueError("Incorrect atom type number");
        return;
    }
    itype --;
    Xscatfac[itype].resize(9);
    Xscatfac[itype][0] = c;
    Xscatfac[itype][1] = a1;
    Xscatfac[itype][2] = b1;
    Xscatfac[itype][3] = a2;
    Xscatfac[itype][4] = b2;
    Xscatfac[itype][5] = a3;
    Xscatfac[itype][6] = b3;
    Xscatfac[itype][7] = a4;
    Xscatfac[itype][8] = b4;
    show_scat(type, itype);
}
