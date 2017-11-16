# include "1D_MainFun.h"

double Field_correction(double *,double,double *,double,double,double,char *,double *);

int main(void)
{
    extern double   Ni[Nmax][LEN+2],Pgas[LEN+2],Tgas[LEN+2],Ngas[LEN+2],Rogas[LEN+2],Hgas[LEN+2],
                    Ne[LEN+2][NEmax],Nel[LEN+2],Te[LEN+2],Tv[LEN+2],
                    Ez,Er[LEN+1],Iexp;

    extern double tau,dt;//dte;
    extern int N;

    int Nedf,Nchem;
    int i,nt,dot,dot1,Nt,n;
    double tic=0.0;

    init_read();
	init_gasDBparser();
	init_data();

	mesh_calc(Len);
	Transport_GFcalc(Geom);

    Nedf = EEDF_read_CS(N);
	Nchem = chem_make_react(Nedf);
	chem_read_react(Nchem,N);

    double Kel[LEN+2][Nedf],Kch[LEN+2][Nchem];
    double dTe[LEN+2],dTgas[LEN+2],dNel[LEN+2],dEzN[LEN+2],Icalc=0.0;
    double Di[N][LEN+2],Mui[N][LEN+2],Lam[LEN+2];
    double Del[LEN+2],Muel[LEN+2],Jel[LEN+2];

    gas_LenPrint(&Ni[0][0],N,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,Ez,Jel,Icalc,&Ne[0][0],0.0);
    gas_TimePrint(&Ni[0][1],N,Pgas[1],Tgas[1],Ngas[1],Nel[1],Te[1],Tv[1],Ez,0.0,Icalc,0.0);

	Nt = int(tau/dt);

	dot = 0, dot1 = 0;
	int Ndot1 = int(Ndots/5);
	//nt = 0;
	do//for(nt=0;nt<Nt;nt++)//Nt
	{
		dot += 1;
		dot1 += 1;

		//1DPoisson_SORsolve(Fi,&Ni[0][0]);

        for(i=0;i<=LEN+1;i++)
        {
            if((tic==0.0) || (dot1==Ndot1) || (fabs(dEzN[i])>0.5))
            {
                EEDF_calc(&Ne[i][0],&Ni[0][i],N,&Te[i],&dTe[i],Ez,Tgas[i],Nel[i],&Del[i],&Muel[i],&Jel[i],1.e-12,tic,dot);
                EEDF_const_calc(&Ne[i][0],N,&Kch[i][0],Nedf,Nel[i],tic);//&Kel[i][0]
            }

            if((tic==0.0) || (dTgas[i]>10.0) || (dTe[i]>0.05) || (dot==Ndots))//Ndots|| (dTe[i]>0.01)
                chem_const(&Kch[i][0],&Kel[i][0],Nedf,Nchem,N,Te[i],Tgas[i],tic);
            chem_runge_kutta4(&Ni[0][i],N,&Kch[i][0],Nchem,dt,tic,dot1);
        }

        if((tic==0.0) || (dot==Ndots) || (dTgas[i]>10.0) || (dTe[i]>0.05))
            Trasport_coefs_calc(N,&Ni[0][0],Del,Muel,&Di[0][0],&Mui[0][0],Lam,Pgas,Tgas,Te);

        for(n=0;n<N;n++)
            SpecTransport(n,&Ni[n][0],&Di[n][0],&Mui[n][0],Ngas,Tgas,Te,Ez,dt);

        HeatTransport(Hgas,Tgas,Lam,&Ni[0][0],&Di[0][0],N,Ez,Jel,Tw,dt);

        //TeTransport(&Ni[0][0],Te,Lel,Vel,dt);

        for(i=1;i<=LEN;i++)
            gas_TP_calc(&Ni[0][i],N,&Pgas[i],&Tgas[i],&dTgas[i],&Ngas[i],&Rogas[i],&Hgas[i],&Nel[i],&dNel[i],&Tv[i],Ez,&dEzN[i]);

        tic += dt;

        //Writing_data***********************************************
        if(dot==Ndots)
        {
            gas_LenPrint(&Ni[0][0],N,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,Ez,Jel,Icalc,&Ne[0][0],tic);
            gas_TimePrint(&Ni[0][1],N,Pgas[1],Tgas[1],Ngas[1],Nel[1],Te[1],Tv[1],Ez,Jel[1],Icalc,tic);
            dot = 0;
        }

        //Control_field_correction***********************************
        if((dot==Ndots-1) || (dot1==Ndot1) || (dEzN[1]>0.5 && dot1==int(Ndot1/10)))//(fabs(1.0-(Icalc/Iexp))>0.05
            Ez = Field_correction(&Icalc,Iexp,Jel,Ez,Ngas[1],dEzN[1],Geom,&dt);

        if(dot1==Ndot1)
            dot1 = 0;

	}while(tic<tau);

	return 0;
}
double Field_correction(double *Icalc,double Iexp,double *Jel,double Eold,double Ngas,double dEzN,char *Geom,double *dt)
{
    //Ñåòêà ïî äëèíå:
	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    //Regime_Iexp=const***********************************

    double Enew,Jsum = 0.0;
    if(!strcmp(Geom,"axial"))
    {
        for(int i=1;i<=LEN;i++)
            Jsum += Jel[i]*0.5*(l[i+1]+l[i])*(l[i+1]-l[i]);
        Jsum *= 2*pi;//[ÑÃÑ/ñ]
    }
    else if(!strcmp(Geom,"cartesian"))
    {
        for(int i=1;i<=LEN;i++)
            Jsum += Jel[i]*(l[i+1]-l[i]);
        Jsum *= Hght;//[ÑÃÑ/ñ]
    }
    Jsum *= eKl/e; //[Êë/ñ=A]

    double th;
    double ENmin = 10.0;///[Td]

    Enew = Eold;
    if(fabs(1.0-(*Icalc/Jsum))<=0.005)///Stable_case
    {
        th = 1.005;
        if(Iexp/Jsum>th)
        {
            Enew = Eold*th;
            printf("Stable_Ez-correction\n");
        }
        if(Jsum/Iexp>th) //&& (fabs(1.0-(*Icalc/Jsum))<0.001))
        {
            Enew = Eold/th;
            printf("Stable_Ez-correction\n");
        }
    }
    if((fabs(1.0-(*Icalc/Jsum))>0.005) && Eold*Eabs*1e17/Ngas>ENmin)///Unstable_case_1(fabs(1.0-(*Icalc/Jsum))>0.025)
    {
        th = 1.005;
        if(Jsum/Iexp>1.5)
        {
            //Enew = Eold/th;
            Enew = Eold*Iexp/Jsum;

            //*dt *= 1/th;
            printf("Unstable_Ez-correction_1\n");
        }
        else if(Jsum/Iexp>th)
        {
            Enew = Eold/th;
            //Enew = Eold*Iexp/Jsum;

            //*dt *= 1/th;
            printf("Unstable_Ez-correction_1\n");
        }
        else if(Iexp/Jsum>th)
        {
            Enew = Eold*th;

            printf("Unstable_Ez-correction_1\n");
        }
    }
    if(dEzN>0.5 && Eold*Eabs*1e17/Ngas>ENmin)///Unstable_case_2
    {
        th = 1.1;
        if(Jsum/Iexp>th)
        {
            //Enew = Eold/th;

            Enew = Eold*Iexp/Jsum;

            //*dt *= 1/th;
            printf("Unstable_Ez-correction_2\n");
        }
    }

    //Enew = Eold*Iexp/Jsum;

    *Icalc = Jsum;///[A]

    if(Enew!=Eold)
        printf("Ez/N = %.1lf[Td]\tEz = %.1f[B/cm]\n",Enew*Eabs*1e17/Ngas,Enew*Eabs);

    return Enew;
}
