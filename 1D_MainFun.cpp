# include "1D_MainFun.h"

int main(void)
{
	init_read();

	mesh_read();
	mesh_calc();

	init_data();

	int nt,dot;
	Nt = int(tau/dt);
	Nte = int(dt/dte);

	int dot = 0, Ndots = 100;
	for(nt=0;nt<=Nt;nt++)
	{
		dot += 1;

		for(nte=0;nte<=Nte;nte++)
			1DPoisson_calc(nte);

		if(dot==Ndots)
			dot = 0;
	}

	return 0;
}
