# include "1D_MainFun.h"

double Field_correction(double *,double,double *,double,double,double,char *,double *,int);

int main(void)
{
    extern double   Ni[Nmax][LEN+2],Xi[Nmax][LEN+2],Pgas[LEN+2],Tgas[LEN+2],Ngas[LEN+2],Rogas[LEN+2],Hgas[LEN+2],
                    Ne[LEN+2][NEmax],Nel[LEN+2],Te[LEN+2],Tv[LEN+2],
                    Ez,Er[LEN+1],Fir[LEN+2],Iexp;

    extern double tau,dt;//dte;
    extern int N,Nneg,Npos;
    extern bool EQP;

    int Nedf,Nchem;
    int i,n,nt,dot,dot1,Ndot1;
    double tic;

    init_read();
	init_gasDBparser();
	tic = init_data();

	Mesh_calc(Len);
	Mesh_GFcalc(Geom);

    Nedf = EEDF_read_CS(N);
	Nchem = chem_make_react(Nedf);
	chem_read_react(Nchem,N);

    double Kel[LEN+2][Nedf],Kch[LEN+2][Nchem],Wrad[LEN+2];
    double dTe[LEN+2],dTgas[LEN+2],dNel[LEN+2],dEzN[LEN+2],Icalc=0.0;
    double Di[N][LEN+2],Mui[N][LEN+2],Lam[LEN+2];
    double Del[LEN+2],Muel[LEN+2],Jel[LEN+2];

    gas_LenPrint(&Ni[0][0],N,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,Ez,Er,Fir,Jel,Icalc,&Ne[0][0],tic);
    gas_TimePrint(&Ni[0][1],N,Pgas[1],Tgas[1],Ngas[1],Nel[1],Te[1],Tv[1],Ez,0.0,Icalc,tic);

    for(i=1;i<=LEN;i++)
    {
        Wrad[i] = 0.0;
        Jel[i] = 0.0;
    }

	dot = 0, dot1 = 0, nt = 0;
	Ndot1 = int(Ndots/150);//20);///75);
	do
	{
		dot += 1;
		dot1 += 1;

        //if((nt==0) || (dot1==Ndot1))//(EQP==true)
            //Poisson_SORsolve(Fir,Er,&Ni[0][0],Npos,Nneg);

        if((nt==0) || (dot==Ndots) || (dTgas[1]>10.0) || (dTe[1]>0.1))///???????
            Trasport_coefs_calc(N,Npos+Nneg,&Ni[0][0],Del,Muel,&Di[0][0],&Mui[0][0],Lam,Pgas,Tgas,Te,EQP);

        TransportBoundary_mod(N,&Ni[0][0],&Xi[0][0],&Di[0][0],Ngas,Tgas,Pgas,Te,Tw);

        HeatTransport(Hgas,Ngas,Tgas,Lam,&Ni[0][0],&Xi[0][0],&Di[0][0],N,Ez,Jel,Wrad,dt);

        for(n=0;n<N;n++)
            SpecTransport(n,&Xi[n][0],&Ni[n][0],&Di[n][0],&Mui[n][0],Ngas,Tgas,Te,Er,EQP,dt);

        for(i=1;i<=LEN;i++)
        {
            if((nt==0) || (dot1==Ndot1) || (fabs(dEzN[i])>1.0))//|| || (dot1==Ndot1)
            {
                EEDF_calc(&Ne[i][0],&Ni[0][i],N,&Te[i],&dTe[i],Ez,Tgas[i],Nel[i],&Del[i],&Muel[i],&Jel[i],1.e-12,tic,dot);
                EEDF_const_calc(&Ne[i][0],N,&Kch[i][0],Nedf,Nel[i],tic);//&Kel[i][0]
            }

            if((nt==0) || (dTgas[i]>10.0) || (dTe[i]>0.1) || (dot==Ndots))//Ndots|| (dTe[i]>0.01)
                chem_const(&Kch[i][0],&Kel[i][0],Nedf,Nchem,N,Te[i],Tgas[i],tic);
            chem_runge_kutta4(&Ni[0][i],N,&Kch[i][0],Nchem,&Wrad[i],dt,tic,dot1);
        }

        //TransportBoundary(N,&Ni[0][0],&Xi[0][0],&Di[0][0],Ngas,Tgas,Te,Tw);

        //TeTransport(&Ni[0][0],Te,Lel,Vel,N+1,dt);

        for(i=1;i<=LEN;i++)
            gas_TP_calc(&Ni[0][i],&Xi[0][i],N,&Pgas[i],&Tgas[i],&dTgas[i],&Ngas[i],&Rogas[i],&Hgas[i],&Nel[i],&dNel[i],&Tv[i],Ez,&dEzN[i]);

        tic += dt;

        //Writing_data***********************************************
        if(dot==Ndots)
        {
            gas_LenPrint(&Ni[0][0],N,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,Ez,Er,Fir,Jel,Icalc,&Ne[0][0],tic);
            gas_TimePrint(&Ni[0][1],N,Pgas[1],Tgas[1],Ngas[1],Nel[1],Te[1],Tv[1],Ez,Jel[1],Icalc,tic);
            gas_SavePrint(&Ni[0][0],N,Pgas,Tgas,&Ne[0][0],Te,Ez,tic);
            dot = 0;
        }

        //Control_field_correction***********************************
        if((dot1==Ndot1) || (dEzN[1]>0.5 && dot1==int(Ndot1/10)))//(fabs(1.0-(Icalc/Iexp))>0.05//&& dot1==int(Ndot1/20
        {
            Ez = Field_correction(&Icalc,Iexp,Jel,Ez,Ngas[1],dEzN[1],Geom,&dt,dot);
            dot1 = 0;
        }

        nt += 1;

	}while(tic<tau);

	return 0;
}
double Field_correction(double *Icalc,double Iexp,double *Jel,double Eold,double Ngas,double dEzN,char *Geom,double *dt,int dot)
{
    //����� �� �����:
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
        Jsum *= 2*pi;//[���/�]
    }
    else if(!strcmp(Geom,"cartesian"))
    {
        for(int i=1;i<=LEN;i++)
            Jsum += Jel[i]*(l[i+1]-l[i]);
        Jsum *= Hght;//[���/�]
    }
    Jsum *= eKl/e; //[��/�=A]

    double th;
    double ENmin = 5.0;///[Td]

    Enew = Eold;
    if(fabs(1.0-(*Icalc/Jsum))<=0.01)///Stable_case
    {
        /*//th = Iexp/Jsum;///1.01;
        if(fabs(Iexp-Jsum)*1.0e3>0.5 && Iexp/Jsum<1.05)///
        {
            th = Iexp/Jsum;
            if(Iexp>Jsum)
                Enew = Eold*th;
            elseEnew = Eold/th;
            //printf("Stable_Ez-correction\n");
        }
        if(Iexp/Jsum>1.05)///(Jsum/Iexp>th) //&& (fabs(1.0-(*Icalc/Jsum))<0.001))
        {
            Enew = Eold*1.01;
            //printf("Stable_Ez-correction\n");
        }*/

        if(fabs(Iexp-Jsum)*1.0e3>0.1)
        {
            if(Iexp>Jsum)
                Enew = Eold*1.0025;
            if(Jsum>Iexp)
                Enew = Eold/1.0025;

            //printf("Stable_Ez-correction\n");
        }

    }
    if((fabs(1.0-(*Icalc/Jsum))>0.01) && Eold*Eabs*1e17/Ngas>ENmin)///Unstable_case_1(fabs(1.0-(*Icalc/Jsum))>0.025)
    {
        th = 1.01;
        if(Jsum/Iexp>1.1)
        {
            //Enew = Eold/th;
            Enew = Eold*Iexp/Jsum;

            //*dt *= 1/th;
            //printf("Unstable_Ez-correction_1\n");
        }
        else if(Jsum/Iexp>th)
        {
            Enew = Eold/th;
            //Enew = Eold*Iexp/Jsum;

            //*dt *= 1/th;
            //printf("Unstable_Ez-correction_1\n");
        }
        else if(Iexp/Jsum>th)
        {
            Enew = Eold*th;

            //printf("Unstable_Ez-correction_1\n");
        }
    }
    if(dEzN>0.5 && Eold*Eabs*1e17/Ngas>ENmin)///Unstable_case_2
    {
        th = 1.2;
        if(Jsum/Iexp>th)
        {
            //Enew = Eold/th;

            Enew = Eold*Iexp/Jsum;

            //*dt *= 1/th;
            //printf("Unstable_Ez-correction_2\n");
        }
    }

    //Enew = Eold*Iexp/Jsum;

    *Icalc = Jsum;///[A]

    if((Enew!=Eold) && dot==Ndots)
        printf("Ez/N = %.1lf[Td]\tEz = %.1f[B/cm]\n",Enew*Eabs*1e17/Ngas,Enew*Eabs);

    return Enew;
}
