# include "1D_MainFun.h"
# include "1DTransport.h"

double GF_C[I+2],GF_L[I+2],GF_R[I+2];
double Di[I+2],Mui[I+2];
int Gf = 0;

void Trasport_coefs_calc()
{
    int i;
    for(i=0;i<I+2;i++)
    {
        Di[i] = 0.0;
        Mui[i] = 0.0;//e/me/Vm_av;
        //from EEDF_calc;
    }
}
int Trasport_GFcalc(char Geom)
{
    //—етка по длине:
	/*            left wall                                                               right wall
                  |                                                                       |
    Ni[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
                  |                                                                       |
                  |                                                                       |
    */

    double lC,lR,lL; //координаты центров €чеек
    double dlC,dlR,dlL; //разности граней [i]-й €чейки, между центрами €чеек слева ([i-1] и [i]) и справа ([i] и [i+1])

    int i;

    if(!strcmp(Geom,"radial"))
    {
        for(i=1;i<I+1;i++)
        {
            lL = 0.5*(l[i-1]+l[i]);
            lC = 0.5*(l[i]+l[i+1]);
            lR = 0.5*(l[i+1]+l[i+2]);

            dlL = lC-lL;
            dlC = l[i+1]-l[i];
            dlR = lR-lC;

            GF_L[i] = l[i]/(lC*dlC*dlL);
            GF_R[i] = l[i+1]/(lC*dlC*dlR);
            GF_C[i] = -(GF_R[i]+GF_L[i]);
        }

    }
    else if(!strcmp(Geom,"axis"))
    {
        for(i=1;i<I+1;i++)
        {
            lL = 0.5*(l[i-1]+l[i]);
            lC = 0.5*(l[i]+l[i+1]);
            lR = 0.5*(l[i+1]+l[i+2]);

            dlL = lC-lL;
            dlC = l[i+1]-l[i];
            dlR = lR-lC;

            GF_L[i] = 1.0/dlC*dlL);
            GF_R[i] = 1.0/dlC*dlR);
            GF_C[i] = -(GF_R[i]+GF_L[i]);
        }
    }
}
void Transport_SWEEPsolve(double Ni, int Gf)
{
	/*
	**************************************************************

	Transport equation (in Drift-Diffusion approximation ) solution with
	implicit SWEEP/Shuttle Method

	Drift-Diffusion Term (DDT) is implemented using exact solution for steady 1D DDT problem:

	d(RoV*F)/dx = dd(D*F)/dx*dx         (*),

	(F(x)-F(0))/(F(L)-F(0))=(exp(Pe*x/L)-1)/(exp(Pe)-1), Pe=RoV*L/D - Peclet Number

    In our Case:

    J[i+1/2] = Vdr[i+0.5]*(N[i]-(N[i+1]-N[i])/(exp(Pe[i+0.5])-1));

    "Power Law Scheme" (Pantakar_1980) used for exp() approximation

    Linearization of Equation (*):

       [i-1]  [i]  [i+1]
    --|--x--|--x--|--x--|--

    A[i]*N[i] = A[i+1]*N[i+1] + A[i-1]*N[i-1]

    Coefficients:
    A[i+1] = D[i+0.5]*F(|Pe[i+0.5]|)+Max(-Vdr[i+0.5];0)
    A[i-1] = D[i-0.5]*F(|Pe[i-0.5]|)+Max(Vdr[i-0.5];0)
    A[i] = A[i+1]+A[i-1]+(Vdr[i+0.5]-Vdr[i-0.5])

    where:
    F(|Pe|) = Max(0;(1-0.1*|Pe|)^5)

    In "RADIAL" case:

    AA[i+1] = A[i+1]*r[i+0.5]/(r[i]*(r[i+0.5]-r[i-0.5]))
    AA[i-1] = A[i-1]*r[i-0.5]/(r[i]*(r[i+0.5]-r[i-0.5]))
    AA[i] = AA[i+1]+AA[i-1]+(r[i+0.5]*Vdr[i+0.5]-r[i-0.5]*Vdr[i-0.5])/(r[i]*(r[i+0.5]-r[i-0.5]))

    Chemistry part considered as explicit (source term)

	**************************************************************
	*/

	double


	if(Gf==0)
        Gf = Transport_GFcalc(char Geom);


    //SWEEP-SHUTTLE_CICLE*****************************************

    //Defining_sweep_coefficients*********************************
    for(n=0;n<N;n++)
    {
        for(i=0;i<=I+1;i++)
        {
            //Diffusion
            if(i==0)
                D_L = D[i];
            else
                D_L = 0.5*(D[i]+D[i-1]);

            if(i==I+1)
                D_R = D[i];
            else
                D_R = 0.5*(D[i+1]+D[i]);

            //Drift-part
            Mui[i]

            lL = 0.5*(l[i-1]+l[i]);
            lC = 0.5*(l[i]+l[i+1]);
            lR = 0.5*(l[i+1]+l[i+2]);

            //DDT_Coefficients:
            Pe_R = fabs(Vdr[i]*(lR-lC)/D_R);
            fPe = pow((1.0-0.1*Pe_R),5.0);
            A_R = D_R*max(0.0,fPe)+max(-Vdr[i],0.0);

            Pe_R = fabs(Vdr[i-1]*(lC-lL)/D_L);
            fPe = pow((1.0-0.1*Pe_L),5.0);
            A_L = D_L*max(0.0,fPe)+max(Vdr[i-1],0.0);

            //Accounting for Geometry Factors:
            AA_R = A_R*GF_R[i];
            AA_L = A_R*GF_L[i];
            AA_C = AA_R+AA_L+GF_R[i]*Vdr_R-GF_L[i]*Vdr_L;

            //SWEEP Coefficients:
            A = -AA_L;//i-1(left)
            B = -AA_R;//i+1(right)
            C = AA_C+1.0/dt;//i(center)

            //RHS-part
            F = Ni[n][i]/dt + Rchem[n][i];

            if(k==0)
            {
                al[i+1] = -B/C;
                bet[i+1] = F/C;
            }
            else if(i==I+1)
            {}
            else
            {
                den = 1.0/(A*al[i]+C);
                al[i+1] = -B*den;
                bet[i+1] = (F-A*bet[i])*den;
            }
        }

        //Boundary_conditions************************************************
        // Ni[n][I+1] = 0.0;
        Ni[n][I+1] = (F-A*bet[I+1])/(A*al[I+1]+C);//уточнить!!!!

        //Reverse_sweep_cycle************************************************
        for(i=I;i>=0;i--)
        {
            Ni[n][i] = al[i+1]*Ni[n][i+1]+bet[i+1];
            if(Ni[n][i]<1.e-30)
                Ni[n][i] = 0.0;
        }

    }

    return 0;

}
void Transport_boundary();
{

}
