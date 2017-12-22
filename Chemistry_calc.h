//Защита от повторного включения заголовка
#ifndef CHEMISTRY_CALC_H
#define CHEMISTRY_CALC_H

int chem_make_react(int);
void chem_read_react(int,int);
void chem_const(double *,double *,int,int,int,double,double,double);
void chem_runge_kutta4(double *,int,double*,int,double*,double,double,int);
void chem_spec_contrib(int,int,double);
int chem_VV_VT_make(int,char *,char *,char *,double *,int);
int chem_VV_VT_const(double *,int,int,char *,double,double);

#endif
