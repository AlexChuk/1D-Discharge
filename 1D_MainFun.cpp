# include "1D_MainFun.h"
# include "Initials.h"
# include "1Dmesh.h"
# include "1DPoisson.h"
# include "1DTransport.h"

int main(void)
{
	init_read();

	1Dmesh_calc(Len);

	init_data();

	int nt,dot;
	Nt = int(tau/dt);
	//Nte = int(dt/dte);

	int dot = 0, Ndots = 100;
	for(nt=0;nt<=Nt;nt++)
	{
		dot += 1;

		1DPoisson_SORsolve(Fi[],1.0);

		1DTransport_SWEEPsolve(Ni);

        if(dot==Ndots)
			dot = 0;
	}

	if exist

	return 0;
}
