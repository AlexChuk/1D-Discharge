# include "1D_MainFun.h"
# include "1DTransport.h"

double GF_C[LEN+2],GF_L[LEN+2],GF_R[LEN+2];
double Di[LEN+2],Mui[LEN+2];
double al_bound[2];
int Gf = 0;

void Transport_GFcalc(char *Geom)
{
    //—етка по длине:
	/*           left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    int i;
    double lC,lR,lL; //координаты центров €чеек
    double dlC,dlR,dlL; //разности граней [i]-й €чейки, между центрами €чеек слева ([i-1] и [i]) и справа ([i] и [i+1])

    if(!strcmp(Geom,"radial"))
    {
        for(i=1;i<LEN+1;i++)
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
        for(i=1;i<LEN+1;i++)
        {
            lL = 0.5*(l[i-1]+l[i]);
            lC = 0.5*(l[i]+l[i+1]);
            lR = 0.5*(l[i+1]+l[i+2]);

            dlL = lC-lL;
            dlC = l[i+1]-l[i];
            dlR = lR-lC;

            GF_L[i] = 1.0/(dlC*dlL);
            GF_R[i] = 1.0/(dlC*dlR);
            GF_C[i] = -(GF_R[i]+GF_L[i]);
        }
    }

}
void Trasport_coefs_calc(double *Ni,double *Di,double *Mui)
{
    for(int i=0;i<LEN+2;i++)
    {
        Di[i] = 100.0;
        Mui[i] = 1.0;//e/me/Vm_av;
        //from EEDF_calc;
    }
}
void Transport_SWEEPsolve(double *Ni,int n,double *Di,double *Mui,double *Ngas,double *Tgas,double *E,double dt,double tic)
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
            al[LEN+1],bet[LEN+1];

    int i,ii;

	//if(Gf==0)
        //Gf = Transport_GFcalc(char Geom);

    //Boundary_conditions*****************************************
    Transport_boundary(Ni,n,Di,Mui,Tgas,E);

    //SWEEP-SHUTTLE_CICLE*****************************************

    //Defining_sweep_coefficients
    for(i=1;i<=LEN;i++)
    {
        //Diffusion
        D_L = D_R;
        if(i==1)
            D_L = 0.5*(Di[i]+Di[i-1]);
        D_R = 0.5*(Di[i+1]+Di[i]);

        //Drift
        Vdr_L = Vdr_R;
        if(i==1)
            Vdr_L = 0.5*(Mui[i]+Mui[i-1])*E[i-1];
        Vdr_R = 0.5*(Mui[i+1]+Mui[i])*E[i];

        //DDT_Coefficients:
        //right-edge
        dL = 0.5*(l[i+2]-l[i]);
        Pe = fabs(Vdr_R*dL/D_R);
        fPe = pow((1.0-0.1*Pe),5.0);
        A_R = D_R*std::max(0.0,fPe)+std::max(-Vdr_R,0.0);

        //left-edge
        dL = 0.5*(l[i+1]-l[i-1]);
        Pe = fabs(Vdr_L*dL/D_L);
        fPe = pow((1.0-0.1*Pe),5.0);
        A_L = D_L*std::max(0.0,fPe)+std::max(Vdr_L,0.0);

        //Accounting for Geometry Factors:
        A_R = A_R*GF_R[i];
        A_L = A_R*GF_L[i];
        A_C = A_R+A_L+GF_R[i]*Vdr_R-GF_L[i]*Vdr_L;

        //SWEEP Coefficients:
        A = -A_L;//[i-1](left_cell)
        B = -A_R;//[i+1](right_cell)
        C = A_C+1.0/dt;//[i](center_cell)
        F = Ni[i]*1.0/dt;//+ Rchem[i];//RHS-part

        if(i==1)//see Transport_boundary();
        {
            al[i-1] = al_bound[0];
            bet[i-1] = 0.0;
        }

        den = 1.0/(A*al[i-1]+C);
        al[i] = -B*den;
        bet[i] = (F-A*bet[i-1])*den;

        if(i==LEN)//see Transport_boundary();]
        {
            al[i] = al_bound[1];
            bet[i] = 0.0;
        }

    }

    //Reverse_sweep_cycle************************************************
    for(i=LEN;i>=0;i--)
    {
        Ni[i] = al[i]*Ni[i+1]+bet[i];
        if(Ni[i]<1.e-30)
            Ni[i] = 0.0;
    }
}
void Transport_boundary(double *Ni,int n,double *Di,double *Mui,double *Tgas,double *E)
{
    double Pe,D,Vd,Vt;

    //left_boundary*************************************************
    if(Gamma[n][0] == 0.0)//Ni[n][0] = Ni[n][1];
    {//Ni[n][0] = Ni[n][1];
        al_bound[0] = 1.0;

        Ni[0] = Ni[1];
    }
    else//equation for DDT with kinetic wall flux
    {
        D = 0.5*(Di[1]+Di[0]);
        Vd = 0.5*(Mui[1]+Mui[0])*E[0];

        Pe = Vd*0.5*(l[2]-l[0])/D;
        Vt = sqrt(8*kb*Tgas[0]/(pi*Mi[n]));

        al_bound[0] = (1.0+0.25*Gamma[n][0]*Pe*Vt/Vd)/(1.0+Pe);

        Ni[0] = Ni[1]/al_bound[0];
    }

    //right_boundary************************************************
    if(Gamma[n][1] == 0.0)//Ni[n][I+1] = Ni[n][I];
    {
        al_bound[1] = 1.0;

        Ni[LEN+1] = Ni[LEN];
    }
    else//equation for DDT with kinetic wall flux
    {
        D = 0.5*(Di[LEN+1]+Di[LEN]);
        Vd = 0.5*(Mui[LEN]+Mui[LEN+1])*E[LEN];

        Pe = Vd*0.5*(l[LEN+2]-l[LEN])/D;
        Vt = sqrt(8*kb*Tgas[LEN+1]/(pi*Mi[n]));

        al_bound[1] = (1.0+0.25*Gamma[n][1]*Pe*Vt/Vd)/(1.0+Pe);//al_bound=Ni[LEN]/Ni[LEN+1]

        Ni[LEN+1] = Ni[LEN]/al_bound[1];
    }

}
void HeatTransport_SWEEPsolve(double Ni, int Gf)
{
    double Lam,Hin;

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
    //Lam = 1.0e+7*(2.714+5.897e-3*Tin)*1.0e-4; //[Ёрг/(с*см^2*K)]Ref. - Nuttall_JRNBS_1957 (N2)
    //Hin += (Qel)*dt - Lam*(Tin-Tw)/pow(Rad/2.4,2.0)*dt;//w/o QE


}
void HeatTransport_boundary()
{

}
