# include "1D_MainFun.h"
# include "1DTransport.h"

double GF_L[I+2],GF_R[I+2];
double Di[I+2],Mui[I+2];
double al_bound[N][2];
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
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    double lC,dlC; //координаты центров €чеек,разности граней [i]-й €чейки
    int i;

    if(!strcmp(Geom,"radial"))
    {
        for(i=1;i<I+1;i++)
        {
            lC = 0.5*(l[i]+l[i+1]);

            dlC = l[i+1]-l[i];

            GF_L[i] = l[i]/(lC*dlC);
            GF_R[i] = l[i+1]/(lC*dlC);
        }

    }
    else if(!strcmp(Geom,"axis"))
    {
        for(i=1;i<I+1;i++)
        {
            dlC = l[i+1]-l[i];

            GF_L[i] = 1.0/dlC;
            GF_R[i] = 1.0/dlC;
        }
    }
}
void Transport_SWEEPsolve(int n,double Ni,int Gf)
{
	/*
    //-----------------------------------------------------------------------------------------------------------------------------

	Transport equation (in Drift-Diffusion approximation ) solution with implicit SWEEP/Shuttle Method:

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

    Coefficients (linear geometry):
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

	//-----------------------------------------------------------------------------------------------------------------------------
	*/

	double  D_L,D_R,Vdr_L,Vdr_R,dL,
            Pe,fPe,
            A_L,A_R,A_C;

    double  A,B,C,F,den,
            al[I+1],bet[I+1];

    int i,ii;

	if(Gf==0)
        Gf = Transport_GFcalc(char Geom);

    //Boundary_conditions*****************************************
    Transport_boundary(n);

    //SWEEP-SHUTTLE_CICLE*****************************************

    //Defining_sweep_coefficients
    for(i=1;i<=I;i++)
    {
        //Diffusion
        D_L = D_R;
        if(i==1)
            D_L = 0.5*(D[n][i]+D[n][i-1]);
        D_R = 0.5*(D[n][i+1]+D[n][i]);

        //Drift
        Vdr_L = Vdr_R;
        if(i==1)
            Vdr_L = 0.5*(Mui[n][i]+Mui[n][i-1])*E[i-1];
        Vdr_R = 0.5*(Mui[n][i+1]+Mui[n][i])*E[i];

        //DDT_Coefficients:
        //right-edge
        dL = 0.5*(l[i+2]-l[i]);
        Pe = fabs(Vdr_R*dL/D_R);
        fPe = pow((1.0-0.1*Pe),5.0);
        A_R = D_R*max(0.0,fPe)+max(-Vdr_R,0.0);

        //left-edge
        dL = 0.5*(l[i+1]-l[i-1]);
        Pe = fabs(Vdr_L*dL/D_L);
        fPe = pow((1.0-0.1*Pe),5.0);
        A_L = D_L*max(0.0,fPe)+max(Vdr_L,0.0);

        //Accounting for Geometry Factors:
        A_R = A_R*GF_R[i];
        A_L = A_R*GF_L[i];
        A_C = A_R+A_L+GF_R[i]*Vdr_R-GF_L[i]*Vdr_L;

        //SWEEP Coefficients:
        A = -A_L;//[i-1](left_cell)
        B = -A_R;//[i+1](right_cell)
        C = A_C+1.0/dt;//[i](center_cell)
        F = Ni[n][i]*1.0/dt + Rchem[n][i];//RHS-part

        if(i==1)//see Transport_boundary();
        {
            al[i-1] = al_bound[n][0];
            bet[i-1] = 0.0;
        }

        den = 1.0/(A*al[i-1]+C);
        al[i] = -B*den;
        bet[i] = (F-A*bet[i-1])*den;

        if(i==I)//see Transport_boundary();]
        {
            al[i] = al_bound[n][1];
            bet[i] = 0.0;
        }

    }

    //Reverse_sweep_cycle************************************************
    for(i=I;i>=0;i--)
    {
        Ni[n][i] = al[i]*Ni[n][i+1]+bet[i];
        if(Ni[n][i]<1.e-30)
            Ni[n][i] = 0.0;
    }

    return 0;
}
void Transport_boundary(int n)
{
    double Pe,D,Vd,Vt;

    //left_boundary*************************************************
    if(Gamma[n][0] == 0.0)//Ni[n][0] = Ni[n][1];
        al_bound[n][0] = 1.0;
    else//equation for DDT with kinetic wall flux
    {
        D = 0.5*(D[n][1]+D[n][0]);
        Vd = 0.5*(Mui[n][1]+Mui[n][0])*E[0];

        Pe = Vd*0.5*(l[2]-l[0])/D;
        Vt = sqrt(8*kb*Tgas[0]/(pi*Mi[n]));

        al_bound[n][0] = (1.0+0.25*Gamma[n][1]*Pe*Vt/Vd)/(1.0+Pe);
    }

    //right_boundary************************************************
    if(Gamma[n][1] == 0.0)//Ni[n][I+1] = Ni[n][I];
        al_bound[n][1] = 1.0;
    else//equation for DDT with kinetic wall flux
    {
        D = 0.5*(D[n][I+1]+D[n][I]);
        Vd = 0.5*(Mui[n][I]+Mui[n][I+1])*E[I];

        Pe = Vd*0.5*(l[I+2]-l[I])/D;
        Vt = sqrt(8*kb*Temp[I+1]/(pi*Mi[n]));

        al_bound[n][1] = (1.0+0.25*Gamma[n][1]*Pe*Vt/Vd)/(1.0+Pe);
    }

}
void HeatTransport_SWEEPsolve(double Ni, int Gf)
{
    //—етка по длине:
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    //”чет нагрева электронами в упругих соударени€х + нагрев газа полем + теплопроводность на стенку трубки
    Lam = 1.0e+7*(2.714+5.897e-3*Tin)*1.0e-4; //[Ёрг/(с*см^2*K)]Ref. - Nuttall_JRNBS_1957 (N2)
    Hin += (Qel)*dt - Lam*(Tin-Tw)/pow(Rad/2.4,2.0)*dt;//w/o QE


}
void HeatTransport_boundary()
{

}
