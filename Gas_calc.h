//Защита от повторного включения заголовка
#ifndef GAS_CALC_H
#define GAS_CALC_H
void gas_HCpSi_calc(double,int);
void gas_TP_calc(double *,double *,int,double *,double *,double *,double *,double *,double *,double *,double *,double *,double,double *);
void gas_LenPrint(double *,int,double *,double *,double *,double *,double *,double *,double *,double *,double,double*,double,double*,double);
void gas_TimePrint(double *,int,double,double,double,double,double,double,double,double,double,double);
void gas_SavePrint(double *,int,double *,double *,double *,double *,double,double);
#endif
