# include "1D_MainFun.h"

//void GasData_print(int);

int main(void)
{
    extern double   Ni[Nmax][LEN+2],Roi[Nmax][LEN+2],Xi[Nmax][LEN+2],
                    Pgas,Tgas[LEN+2],Ngas[LEN+2],Rogas[LEN+2],Hgas[LEN+2],
                    Ne[LEN+2][NEmax],Nel[LEN+2],Te[LEN+2],Tv[LEN+2],
                    E[LEN+1],E_N[LEN+1];

    extern double tau,dt,tic;//dte;
    extern int N;

    int Nedf,Nchem;
    int nt,dot,Nt;

    init_read();
	init_gasDBparser();
	init_data();

	mesh_calc(Len);
	mesh_GFcalc(Geom);

    Nedf = EEDF_read_CS(N);
	Nchem = chem_make_react(Nedf);
	chem_read_react(Nchem);

    double Kel[LEN][Nedf],Kch[LEN][Nchem];

	Nt = int(tau/dt);
	//Nte = int(dt/dte);

	dot = 0;
	for(nt=0;nt<=Nt;nt++)
	{
		tic += dt;
		dot += 1;

		//1DPoisson_SORsolve(Fi,&Ni[0][0]);

        for(int i=1;i<=LEN;i++)
        {
            if((dot==Ndots) || (nt==0))
            {
                EEDF_calc(&Ne[i][0],&Ni[0][i],N,Te,E[i],Tgas[i],Nel[i],tic,dot);
                EEDF_const_calc(N,&Ne[i][0],&Kel[i][0],Nel[i]);
            }

            /*if((dTgas[i]>10.0) || (nt==0))
                chem_const(N,&Kch[i][0],&Kel[i][0],Te[i],Tgas[i],tic);
            chem_runge_kutta4(&Ni[0][i],N,dt,tic);*/
        }

		/*
		for(n=0;n<N;n++)
            1DTransport_SWEEPsolve(n,Ni);

        1DHeatTransport_SWEEPsolve(Ni);

		for(i=1;i<=LEN;i++)
            gas_TP_calc(&Ni[n][i],Tgas,Pgas,Hgas)
        */

        if(dot==Ndots)
			dot = 0;
	}

    //Writing_data***********************************************
	if(dot==Ndots)
		gas_print(&Ni[0][0],N,Nel,Tgas,Pgas,Hgas,Rogas,Te,Tv,tic);

	return 0;
}
