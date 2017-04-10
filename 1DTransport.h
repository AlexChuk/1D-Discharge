#ifndef TRANSPORT_H
#define TRANSPORT_H

void Trasport_coefs_calc();
int Trasport_GFcalc(char);
void Transport_SWEEPsolve(int,double *,int);
void Transport_boundary(int);
void HeatTransport_SWEEPsolve(double *,int);
void HeatTransport_boundary();


#endif // TRANSPORT_H
