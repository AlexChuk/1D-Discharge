# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <algorithm>

#ifndef MAINFUN_H
#define MAINFUN_H

# define I 100  //dots per axis
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
	eV = 11605,//1 �� � ��������� �
	p0 = 1333.22; // ����-� �������� ���� --> ��� [���/��^3]

extern int N,Nt,Nte,Nchem,Nedf,Ndots,NR;

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

extern double Fi[I+2];
extern double Di[I+2],Mui[I+2];
extern double Ne[NEmax],Ni[Nmax][I+2],Mi[Nmax],LJi[Nmax][2],Roi[Nmax][I+2],
              Pgas[I+2],Tgas[I+2],Ngas[I+2],Hgas[I+2],Rogas[I+2];


#endif
