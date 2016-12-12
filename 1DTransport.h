#ifndef 1DTRANSPORT_H_INCLUDED
#define 1DTRANSPORT_H_INCLUDED

void 1DTrasport_coefs_calc();
int 1DTrasport_GFcalc(char);
void 1DTransport_SWEEPsolve(double *,double);
void 1DTransport_boundary();

#endif // 1DTRANSPORT_H_INCLUDED
