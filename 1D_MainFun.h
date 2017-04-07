# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <algorithm>

# include "Initials.h"
//# include "1Dmesh.h"
//# include "1DPoisson.h"
//# include "1DTransport.h"
# include "Chemistry_calc.h"
# include "EEDF_calc.h"
# include "Gas_calc.h"


#ifndef MAINFUN_H
#define MAINFUN_H

# define LEN 20  //dots per axis
# define NEmax 1000
# define Nmax 100
# define CSmax 100
# define NRmax 1000

const double
	pi = 3.141592653589,
	c = 2.997924562e+10,
	ma = 1.67e-24,//[�]
	me = 9.1e-28,//[�]
	e = 4.8e-10,//[���]
	E0 = 300,//E[�/��]=300E[���.���]
	Na = 6.022e+23,
	kb = 1.38e-16,
	eV_K = 11605,//1 �� � ��������� �
	p0 = 1333.22, // ����-� �������� ���� --> ��� [���/��^3]
	exact = 1.0e-5,//�������� �����
    cm_eV = 8065.5447;//����-� �������� cm-1 --> ��


extern int NR,Ndots;//N,Nt,Nte,Nchem,Nedf,Ndots,
extern char Spec[Nmax][10],Spec_R[Nmax][10],Geom[10];
extern double Mi[Nmax],HCpSi[3][Nmax],CXi[Nmax][2][8];
extern double Emax,dE,dEev;


//extern double tau,dt;//dte;

/*
old declaration:

extern int N,Nt,Nte,Nchem,Nedf,Ndots,NR;

extern double E,E_N,Nel,Ee,Te,Tv,Muel,Jel,Qel,QE;
extern double dTgas,dTe,dNel;
extern double Len,Tw,Lam;
extern double tau,dt,dte;

extern double Kel[CSmax];
extern double HCpSi[3][Nmax];

extern char Rtype[NRmax][10];

extern int v0,vlen;
extern double V0,Vlen;

*/


#endif
