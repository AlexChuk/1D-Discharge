# include "1D_MainFun.h"

int main(void)
{
    extern double   Ni[Nmax][LEN+2],Pgas[LEN+2],Tgas[LEN+2],Ngas[LEN+2],Rogas[LEN+2],Hgas[LEN+2],
                    Ne[LEN+2][NEmax],Nel[LEN+2],Te[LEN+2],Tv[LEN+2],
                    E[LEN+1],E_N[LEN+1],
                    dTgas[LEN+2],dNel[LEN+2];

    extern double tau,dt;//dte;
    extern int N;

    int Nedf,Nchem;
    int i,nt,dot,Nt;
    double tic=0.0;

    init_read();
	init_gasDBparser();
	init_data();

	mesh_calc(Len);
	mesh_GFcalc(Geom);

	gas_LenPrint(&Ni[0][0],N,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,E,&Ne[0][0],0.0);
    gas_TimePrint(&Ni[0][1],N,Pgas[1],Tgas[1],Ngas[1],Nel[1],Te[1],Tv[1],E[1],0.0);

    Nedf = EEDF_read_CS(N);
	Nchem = chem_make_react(Nedf);
	chem_read_react(Nchem,N);

    double Kel[LEN+2][Nedf],Kch[LEN+2][Nchem],dTe[LEN+2];

	Nt = int(tau/dt);
	//Nte = int(dt/dte);

	dot = 0;
	for(nt=0;nt<Nt;nt++)//Nt
	{
		dot += 1;

		//1DPoisson_SORsolve(Fi,&Ni[0][0]);

        //if((dot==Ndots) || (nt==0))
        {
            for(i=1;i<=LEN;i++)
            {
                //if(nt==0)
                {
                    EEDF_calc(&Ne[i][0],&Ni[0][i],N,&Te[i],&dTe[i],E[i],Tgas[i],Nel[i],tic,dot);
                    EEDF_const_calc(&Ne[i][0],N,&Kel[i][0],Nedf,Nel[i],tic);
                }

                //if(nt==0)//((dTgas[i]>10.0) || (nt==0))
                    chem_const(&Kch[i][0],&Kel[i][0],Nchem,N,Te[i],Tgas[i],tic);
                    chem_runge_kutta4(&Ni[0][i],N,&Kch[i][0],Nchem,dt,tic,dot);
            }

        }

		/*
		for(n=0;n<N;n++)
            Transport_SWEEPsolve(&Ni[n][0],N,Ngas,Tgas,dt/2.0,tic);

        HeatTransport_SWEEPsolve(&Ni[0][0],N,Ngas,Tgas,dt/2.0,tic);
        */

        //for(i=0;i<=LEN+1;i++)
            //gas_TP_calc(&Ni[0][i],N,Pgas,Tgas,dTgas,Ngas,Rogas,Hgas,Nel,dNel,Tv);

        //Writing_data***********************************************
        if(dot==Ndots)
        {
            gas_LenPrint(&Ni[0][0],N,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,E,&Ne[0][0],tic);
            gas_TimePrint(&Ni[0][1],N,Pgas[1],Tgas[1],Ngas[1],Nel[1],Te[1],Tv[1],E[1],tic);
            dot = 0;
        }

        tic += dt;

	}


	return 0;
}
