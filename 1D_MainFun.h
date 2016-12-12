# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>

#ifndef 1D_MAINFUN_H
#define 1D_MAINFUN_H

# define I 100  //dots per axis

const double
	pi = 3.141592653589,
	c = 2.997924562e+10,
	ma = 1.67e-24,//[г]
	me = 9.1e-28,//[г]
	e = 4.8e-10,//[СГС]
	E0 = 300,//E[В/см]=300E[абс.СГС]
	Na = 6.022e+23,
	kb = 1.38e-16,
	eV = 11605,//1 эВ в кельвинах К
	p0 = 1333.22; // коэф-т перевода Торр --> СГС [эрг/см^3]

extern int N,Nt,Nte,Nchem,Nedf,Ndots,NR;

extern double Ne[NEmax],Ni[Nmax],Mi[Nmax],LJi[Nmax][2],Roi[Nmax],Pgas,Tgas,Ngas,Hgas,Rogas;
extern double E,E_N,Nel,Ee,Te,Tv,Vdr,Muel,Jel,Qel,QE;
extern double dTgas,dTe,dNel;
extern double Len,Tw,Lam;
extern double tau,dt,dte;
extern double Emax,dE,dEev;
extern double CXi[Nmax][2][8];
extern double Kel[CSmax];
extern double HCpSi[3][Nmax];

extern char Rtype[NRmax][10];
extern char Spec[Nmax][10],Spec_R[Nmax][10],Geom[10];

extern int v0,vlen;
extern double V0,Vlen;

#endif
