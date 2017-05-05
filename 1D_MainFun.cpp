# include "1D_MainFun.h"

int main(void)
{
    extern double   Ni[Nmax][LEN+2],Pgas[LEN+2],Tgas[LEN+2],Ngas[LEN+2],Rogas[LEN+2],Hgas[LEN+2],
                    Ne[LEN+2][NEmax],Nel[LEN+2],Te[LEN+2],Tv[LEN+2],
                    E[LEN+1],E_N[LEN+1];

    extern double tau,dt;//dte;
    extern int N;

    int Nedf,Nchem;
    int i,nt,dot,Nt,n;
    double tic=0.0;

    init_read();
	init_gasDBparser();
	init_data();

	mesh_calc(Len);
	//Transport_GFcalc(Geom);

	gas_LenPrint(&Ni[0][0],N,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,E,&Ne[0][0],0.0);
    gas_TimePrint(&Ni[0][1],N,Pgas[1],Tgas[1],Ngas[1],Nel[1],Te[1],Tv[1],E[1],0.0);

    Nedf = EEDF_read_CS(N);
	Nchem = chem_make_react(Nedf);
	chem_read_react(Nchem,N);

    double Kel[LEN+2][Nedf],Kch[LEN+2][Nchem];
    double dTe[LEN+2],dTgas[LEN+2],dNel[LEN+2];
    double Di[N][LEN+2],Mui[N][LEN+2];

	Nt = int(tau/dt);

	dot = 0;
	for(nt=0;nt<Nt;nt++)//Nt
	{
		dot += 1;

		//1DPoisson_SORsolve(Fi,&Ni[0][0]);

        for(i=1;i<=LEN;i++)
        {
            if((dot==Ndots) || (nt==0))
            {
                EEDF_calc(&Ne[i][0],&Ni[0][i],N,&Te[i],&dTe[i],E[i],Tgas[i],Nel[i],1.e-10,tic,dot);
                EEDF_const_calc(&Ne[i][0],N,&Kel[i][0],Nedf,Nel[i],tic);
            }

            if((dTgas[i]>10.0) || (dTe[i]>0.1) || (nt==0) || (dot==Ndots))
                chem_const(&Kch[i][0],&Kel[i][0],Nchem,N,Te[i],Tgas[i],tic);
            chem_runge_kutta4(&Ni[0][i],N,&Kch[i][0],Nchem,dt/2.0,tic,dot);//dt/2.0
        }

        /*for(n=0;n<N;n++)
        {
            if((dot==Ndots) || (nt==0))
                Trasport_coefs_calc(&Ni[n][0],&Di[n][0],&Mui[n][0]);
            Transport_SWEEPsolve(&Ni[n][0],n,&Di[n][0],&Mui[n][0],Ngas,Tgas,E,dt/2.0,tic);
        }*/

        /*
        HeatTransport_SWEEPsolve(&Ni[0][0],N,Ngas,Tgas,dt/2.0,tic);
        */

        for(i=1;i<=LEN;i++)
            gas_TP_calc(&Ni[0][i],N,&Pgas[i],&Tgas[i],&dTgas[i],&Ngas[i],&Rogas[i],&Hgas[i],&Nel[i],&dNel[i],&Tv[i]);

        tic += dt;

        //Writing_data***********************************************
        if(dot==Ndots)
        {
            gas_LenPrint(&Ni[0][0],N,Pgas,Tgas,Ngas,Hgas,Rogas,Nel,Te,Tv,E,&Ne[0][0],tic);
            gas_TimePrint(&Ni[0][1],N,Pgas[1],Tgas[1],Ngas[1],Nel[1],Te[1],Tv[1],E[1],tic);
            dot = 0;
        }

        //tic += dt;

	}


	return 0;
}
