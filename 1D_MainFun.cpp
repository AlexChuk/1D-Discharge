# include "1D_MainFun.h"
# include "Initials.h"
# include "1Dmesh.h"
# include "1DPoisson.h"
# include "1DTransport.h"

int main(void)
{
	init_read();

	mesh_calc(Len);

	init_data();

	int nt,dot;
	Nt = int(tau/dt);
	//Nte = int(dt/dte);

	dot = 0;
	for(nt=0;nt<=Nt;nt++)
	{
		dot += 1;

		Poisson_SORsolve(Fi,1.0);

		//Transport_SWEEPsolve(Ni);

        if(dot==Ndots)
			dot = 0;
	}

	return 0;
}
