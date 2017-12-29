#ifndef TRANSPORT_H
#define TRANSPORT_H

void Transport_GFcalc(char *);
void Trasport_coefs_calc(int,double *,double *,double *,double *,double *,double *,double *,double *,double *);
double* Transport_SWEEPsolve(double *,double *,double *,double *,double *,double *,double);
double* Transport_SWEEPsolve_mod(double *,double *,double *,double *,double *,double *,double *,double);
void SpecTransport(int,double *,double *,double *,double *,double *,double *,double *,double,double);
void HeatTransport(double *,double *,double *,double *,double *,double *,double *,int,double,double *,double *,double);
void TeTransport(double *,double *,double *,double *,int,double);
void TransportBoundary(int,double *,double *,double *,double *,double *,double *,double);
void TransportBoundary_mod(int,double *,double *,double *,double *,double *,double *,double *,double);

#endif // TRANSPORT_H
