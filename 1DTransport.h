#ifndef TRANSPORT_H
#define TRANSPORT_H

void Transport_GFcalc(char *);
void Trasport_coefs_calc(double *,double *,double *);
void Transport_SWEEPsolve(double *,int,double *,double *,double *,double *,double *,double,double);
void Transport_boundary(double *,int ,double *,double *,double *,double *);
void HeatTransport_SWEEPsolve(double *,int);
void HeatTransport_boundary();


#endif // TRANSPORT_H
