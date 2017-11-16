# include "1D_MainFun.h"
# include "1DTransport.h"
# include "Gas_calc.h"

double GF_C[LEN+2],GF_L[LEN+2],GF_R[LEN+2];
double al_bound[2];
int Gf = 0;

void Transport_GFcalc(char *Geom)
{
    //Сетка по длине:
	/*           left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    int i;
    double lC,lR,lL; //координаты центров ячеек
    double dlC,dlR,dlL; //разности граней [i]-й ячейки, между центрами ячеек слева ([i-1] и [i]) и справа ([i] и [i+1])

    if(!strcmp(Geom,"axial"))//"radial"
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
    else if(!strcmp(Geom,"cartesian"))
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
void Trasport_coefs_calc(int N,double *Ni,double *Del,double *Muel,double *Di,double *Mui,double *Lam,double *Pgas,double *Tgas,double *Te)
{
    double Damb[LEN+2],PTi,PTn;
    PTi = 760*p0/(300*300);
    PTn = 760*p0/pow(300,1.5);

    for(int n=0;n<N;n++)
    {
        for(int i=0;i<LEN+2;i++)
        {
            /*if(n==0)
            {
                Di[n*(LEN+2)+i] = 0.001*Del[i];
                Mui[n*(LEN+2)+i] = 0.0;//1e-5*Muel[i];//1.0;//e/me/Vm_av;

            }*/
            //Ambipolar_Case:**********************************************
            if(n==0)
                Damb[i] = 0.0;
            else if(n<=3)
            {
                if(n==1 || n==3)
                    Di[n*(LEN+2)+i] = 0.07;///[cm2/s]--из справочника "Физ.Величины"_стр.433
                if(n==2)
                    Di[n*(LEN+2)+i] = 0.058;///[cm2/s]--из справочника "Физ.Величины"_стр.433
                //if()
                    //Di[n*(LEN+2)+i] *= Ni[n*(LEN+2)+i]/Ni[i];
                Damb[i] += Di[n*(LEN+2)+i];

                if(n==3)
                {
                    Damb[i] *= PTi*Tgas[i]*Tgas[i]/Pgas[i]*(1.0+Te[i]*eV_K/Tgas[i]);
                    for(int k=0;k<=n;k++)
                        Di[k*(LEN+2)+i] = Damb[i]/n;
                }

                Mui[n*(LEN+2)+i] = 0.0;//1e-5*Muel[i];//1.0;//e/me/Vm_av;
            }
            else
            {
                if(n>=8 && n<=10)
                    Di[n*(LEN+2)+i] = 0.06;///[N,N(2D),N(2P)]_MankModel-data
                else
                Di[n*(LEN+2)+i] = 0.028;///[N2...]_[cm2/s]_MankModel-data

                Di[n*(LEN+2)+i] *= PTn*pow(Tgas[i],1.5)/Pgas[i];

                Mui[n*(LEN+2)+i] = 0.0;//1.0;//e/me/Vm_av;
            }

            ///Lam = 1.0e+7*(2.714+5.897e-3*Tin)*1.0e-4; //[Эрг/(с*см^2*K)]Ref. - Nuttall_JRNBS_1957 (N2)
            if(n==N-1)
                Lam[i] = 1.0e+7*(2.714+5.897e-3*Tgas[i])*1.0e-4;
        }
    }
}
double* Transport_SWEEPsolve(double *NNi,double *Di,double *Vi,double *Si,double dt)
{
    /*
    //-----------------------------------------------------------------------------------------------------------------------------

	Transport equation (in Drift-Diffusion approximation ) solution with implicit SWEEP/Shuttle Method.

	dNi/dt = d(D*(dNi/dx))/dx - d(Vi*Ni)dx + S , S - source term.

	Drift-Diffusion Term (DDT) is implemented using exact solution for steady 1D DDT problem:

	d(RoV*F)/dx = dd(D*F)/dx*dx         (*),

	(F(x)-F(0))/(F(L)-F(0))=(exp(Pe*x/L)-1)/(exp(Pe)-1), Pe=RoV*L/D - Peclet Number

    In our Case:

    1)Species:
    J[i+1/2] = Vdr[i+0.5]*(N[i]-(N[i+1]-N[i])/(exp(Pe[i+0.5])-1));

    "Power Law Scheme" (Pantakar_1980) used for exp() approximation

    Linearization of Equation (*):

    J[i+1/2] = Vdr[i+0.5]*(N[i]-(1-0.1*Pe)^5(N[i+1]-N[i])/Pe), for 0<Pe<=5;

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

	double  D_L,D_R,V_L,V_R,dL,
            Pe,fPe,
            A_L,A_R,A_C;

    double  A,B,C,F,den,
            al[LEN+1],bet[LEN+1];

    //SWEEP-SHUTTLE_CICLE*****************************************

    //Defining_sweep_coefficients
    int i;
    for(i=1;i<=LEN;i++)
    {
        //Defining_sweep_coefficients
        //Diffusion
        D_L = D_R;
        if(i==1)
            D_L = 0.5*(Di[i]+Di[i-1]);
        D_R = 0.5*(Di[i+1]+Di[i]);

        //Drift
        V_L = V_R;
        if(i==1)
            V_L = Vi[i-1];//0.5*(Vi[i]+Vi[i-1]);
        V_R = Vi[i];//0.5*(Vi[i+1]+Vi[i]);

        //DDT_Coefficients:
        //right-edge
        dL = 0.5*(l[i+2]-l[i]);
        Pe = fabs(V_R*dL/D_R);
        fPe = pow((1.0-0.1*Pe),5.0);
        //A_R = D_R*std::max(0.0,fPe)+std::max(-V_R,0.0);
        A_R = D_R*fmax(0.0,fPe)+fmax(-V_R,0.0);

        //left-edge
        dL = 0.5*(l[i+1]-l[i-1]);
        Pe = fabs(V_L*dL/D_L);
        fPe = pow((1.0-0.1*Pe),5.0);
        //A_L = D_L*std::max(0.0,fPe)+std::max(V_L,0.0);
        A_L = D_L*fmax(0.0,fPe)+fmax(V_L,0.0);

        //Accounting for Geometry Factors:
        A_R = A_R*GF_R[i];
        A_L = A_L*GF_L[i];
        A_C = A_R+A_L+GF_R[i]*V_R-GF_L[i]*V_L;

        //SWEEP Coefficients:
        A = -A_L;//[i-1](left_cell)
        B = -A_R;//[i+1](right_cell)
        C = A_C+1.0/dt;//[i](center_cell)
        F = NNi[i]*1.0/dt;//+ Rchem[i];//RHS-part

        if(i==1)//see SpecTransport();
        {
            al[i-1] = al_bound[0];
            bet[i-1] = 0.0;
        }

        den = 1.0/(A*al[i-1]+C);
        al[i] = -B*den;
        bet[i] = (F-A*bet[i-1])*den;

        if(i==LEN)//see SpecTransport();
        {
            al[i] = al_bound[1];
            bet[i] = 0.0;
        }

    }

    //Reverse_sweep_cycle************************************************
    for(i=LEN;i>=0;i--)
    {
        NNi[i] = al[i]*NNi[i+1]+bet[i];
        if(NNi[i]<1.e-30)
            NNi[i] = 0.0;
    }

    return &NNi[0];
}
void SpecTransport(int n,double *Ni,double *Di,double *Mui,double *Ngas,double *Tgas,double *Te,double E,double dt)
{
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    double Vdr[LEN+1],Si[LEN+2];
    double Pe,D,Vd,Vt;

    //Boundary_conditions*****************************************

    //left_boundary***********************************************
    if(Gamma[n][0] == 0.0)//Ni[n][0] = Ni[n][1];
    {//Ni[n][0] = Ni[n][1];
        al_bound[0] = 1.0;

        Ni[0] = Ni[1];
    }
    else//equation for DDT with kinetic wall flux
    {
        D = 0.5*(Di[1]+Di[0]);

        if(n==0)
            Vt = sqrt(8*kb*Te[0]*eV_K/(pi*Mi[n]));
        else
            Vt = sqrt(8*kb*Tgas[0]/(pi*Mi[n]));

        ///if(Vd!=0)
        ///Vd = 0.5*(Mui[1]+Mui[0])*E;
        ///Pe = Vd*0.5*(l[2]-l[0])/D;
        ///al_bound[0] = (1.0+0.25*Gamma[n][0]*Pe*Vt/Vd)/(1.0+Pe);

        Pe = 0.5*(l[2]-l[0])/D;
        al_bound[0] = (1.0+0.25*Gamma[n][0]*Vt*Pe);

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

        if(n==0)
            Vt = sqrt(8*kb*Te[LEN+1]*eV_K/(pi*Mi[n]));
        else
            Vt = sqrt(8*kb*Tgas[LEN+1]/(pi*Mi[n]));

        ///if(Vd!=0)
        ///Vd = 0.5*(Mui[LEN]+Mui[LEN+1])*E;
        ///Pe = Vd*0.5*(l[LEN+2]-l[LEN])/D;
        ///al_bound[1] = (1.0+0.25*Gamma[n][1]*Pe*Vt/Vd)/(1.0+Pe);//al_bound=Ni[LEN]/Ni[LEN+1]

        Pe = 0.5*(l[LEN+2]-l[LEN])/D;
        al_bound[1] = (1.0+0.25*Gamma[n][1]*Vt*Pe);//al_bound=Ni[LEN]/Ni[LEN+1]

        Ni[LEN+1] = Ni[LEN]/al_bound[1];
    }

    //Defining_Drift-coefficient and Source_term******************
    int i;
    for(i=1;i<=LEN;i++)
    {
        if(i==1)
        {
            Vdr[i-1] = 0.5*(Mui[i]+Mui[i-1])*0.0;//E[i-1];
            Si[i-1] = 0.0;
            Si[LEN+1] = 0.0;
        }

        Vdr[i] = 0.5*(Mui[i+1]+Mui[i])*0.0;//*E[i];
        Si[i] = 0.0;
    }

    //Transport_equation_solve************************************
    double *res;
    res = Transport_SWEEPsolve(Ni,Di,Vdr,Si,dt);

    for(i=1;i<=LEN;i++)
        Ni[i] = *(res+i);
}
void HeatTransport(double *Hgas,double *Tgas,double *Lam,double *Ni,double *Di,int N,double E,double *J,double Twall,double dt)
{
    //Сетка по длине:
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    //Учет нагрева электронами в упругих соударениях + нагрев газа полем + теплопроводность на стенку трубки
    //Lam = 1.0e+7*(2.714+5.897e-3*Tin)*1.0e-4; //[Эрг/(с*см^2*K)]Ref. - Nuttall_JRNBS_1957 (N2)
    //Hin += (Qel)*dt - Lam*(Tin-Tw)/pow(Rad/2.4,2.0)*dt;//w/o QE

    double  Lm,QL,QR,JHL,JHR,Wrad,T0,Dn;
    int i,n;

    //Boundary_conditions*****************************************

    //left_boundary***********************************************
    if(!strcmp(Geom,"axial"))
        Tgas[0] = Tgas[1];
    else if(!strcmp(Geom,"cartesian"))
        Tgas[0] = Twall;

    //right_boundary***********************************************
    Tgas[LEN+1] = Twall;

    //Calculation_of_enthalpy*************************************
    for(i=1;i<LEN+1;i++)
    {
        //Учет теплопроводности******************

        QL = QR;
        if(i==1)
        {
            Lm = 0.5*(Lam[i]+Lam[i-1]);
            QL = - Lm*(Tgas[i]-Tgas[i-1]);
        }
        Lm = 0.5*(Lam[i+1]+Lam[i]);
        QR = - Lm*(Tgas[i+1]-Tgas[i]);

        //Учет диффузионного переноса************

        JHL = JHR;
        if(i==1)
        {
            T0 = (Tgas[i]+Tgas[i-1])*0.5;
            gas_HCpSi_calc(T0,N);

            JHL = 0.0;
            for(n=1;n<N;n++)
            {
                Dn = 0.5*(Di[n*(LEN+2)+i]+Di[n*(LEN+2)+i-1]);
                JHL += - HCpSi[0][n]*Dn*(Ni[n*(LEN+2)+i]-Ni[n*(LEN+2)+i-1]);
            }
        }
        JHR = 0.0;
        T0 = (Tgas[i+1]+Tgas[i])*0.5;
        gas_HCpSi_calc(T0,N);
        for(n=1;n<N;n++)
        {
            Dn = 0.5*(Di[n*(LEN+2)+i+1]+Di[n*(LEN+2)+i]);
            JHR += - HCpSi[0][n]*Mi[n]*Dn*(Ni[n*(LEN+2)+i+1]-Ni[n*(LEN+2)+i]);
        }

        //Учет Радиационных потерь тепла*********
        Wrad = 0.0;//

        Hgas[i] += -dt*(GF_R[i]*(QR+JHR)-GF_L[i]*(QL+JHL))+dt*(J[i]*E-Wrad);///{J*E}=[erg/cm3/c]
    }
}
void TeTransport(double *Ne,double *Te,double *Lel,double *Vel,double dt)
{
    double NeTe[LEN+2],De[LEN+2],Ve[LEN+1],Se[LEN+2];

    //Сетка по длине:
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    //Boundary_conditions*****************************************
    //left_boundary**********************
    Te[0] = Te[1];
    NeTe[0] = 1.5*Ne[0]*Te[0];
    al_bound[0] = 1.0;

    //right_boundary**********************
    Te[LEN+1] = Te[LEN];
    NeTe[LEN+1] = 1.5*Ne[LEN+1]*Te[LEN+1];
    al_bound[1] = 1.0;

    //Defining_Drift-coefficient and Source_term******************
    int i;
    double dL;
    for(i=1;i<=LEN;i++)
    {
        NeTe[i] = 1.5*Ne[i]*Te[i];

        if(i==1)
        {
            De[0] = 2.0*Lel[0]/Ne[0]/3.0;
            De[LEN+1] = 2.0*Lel[LEN+1]/Ne[LEN+1]/3.0;
        }
        De[i] = 2.0*Lel[i]/Ne[i]/3.0;//0.666667

        dL = 0.5*(l[i+1]-l[i-1]);

        Ve[i-1] = 5.0*Vel[i-1]/3.0+0.5*(De[i]+De[i-1])*(Ne[i]-Ne[i-1])/dL;   //1.666667
        if(i==LEN)
        {
            dL = 0.5*(l[i+2]-l[i]);
            Ve[i] = 5.0*Vel[i]/3.0+0.5*(De[i]+De[i+1])*(Ne[i+1]-Ne[i])/dL;
        }

        Se[i] = 0.0;//(Jel[i]*E[i])-Ne*Sum(ki*Ni*dEei);
    }

    //проверить_геометрические факторы в поправке для сноса!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //Transport_equation_solve************************************
    double *res;

    res = Transport_SWEEPsolve(NeTe,De,Ve,Se,dt);
    for(i=1;i<=LEN;i++)
        Te[i] = *(res+i)/(1.5*Ne[i]);
}
