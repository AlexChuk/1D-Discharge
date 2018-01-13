//Защита от повторного включения заголовка
#ifndef POISSON_CALC_H
#define POISSON_CALC_H

void Poisson_SORsolve(double *,double *,double *,int,int);
void Poisson_boundary(double *,int,double,int,double);

#endif
