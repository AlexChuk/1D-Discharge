#ifndef TRANSPORT_H
#define TRANSPORT_H

void 1DTrasport_coefs_calc();
int 1DTrasport_GFcalc(char);
void 1DTransport_SWEEPsolve(int,double *,int);
void 1DTransport_boundary(int);
void 1DHeatTransport_SWEEPsolve(double *,int);
void 1DHeatTransport_boundary();


#endif // TRANSPORT_H
