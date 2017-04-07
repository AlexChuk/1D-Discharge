//Защита от повторного включения заголовка
#ifndef CHEMISTRY_CALC_H
#define CHEMISTRY_CALC_H

int chem_make_react(int);
void chem_read_react(int);
void chem_const(int,double *,double *,double,double,double);
void chem_runge_kutta4(double *,int,double,double,int);
void chem_spec_contrib(int, double);
int chem_VV_VT_make(int,char *,char *,char *,double *,int);
int chem_VV_VT_const(int,int,char *,double,double);

#endif
