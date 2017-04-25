//Защита от повторного включения заголовка
#ifndef EEDF_CALC_H
#define EEDF_CALC_H
int EEDF_read_CS(int);
void EEDF_calc(double *,double *,int,double *,double *,double,double,double,double,double,int);//решение уравнения Больцмана
void EEDF_const_calc(double *,int,double *,int,double,double);
void EEDF_print(double *,double,double,double,double);
#endif
