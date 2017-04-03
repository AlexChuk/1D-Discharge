#ifndef TRANSPORT_H
#define TRANSPORT_H

void 1DTrasport_coefs_calc();
int 1DTrasport_GFcalc(char);
void 1DTransport_SWEEPsolve(double *,int);
void 1DTransport_boundary();
void 1DHeatTransport_SWEEPsolve(double *,int);
void 1DHeatTransport_boundary();


#endif // TRANSPORT_H
