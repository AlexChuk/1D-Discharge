# include "1D_MainFun.h"
# include "1DTransport.h"

double GF_C[I+2],GF_L[I+2],GF_R[I+2];
double Di[I+2],Mui[I+2];
int Gf = 0;

void Trasport_coefs_calc()
{
    int i;
    for(i=0;i<I+2;i++)
    {
        Di[i] = 0.0;
        Mui[i] = 0.0;//e/me/Vm_av;
        //from EEDF_calc;
    }
}
int Trasport_GFcalc(char Geom)
{
    //—етка по длине:
	/*            left wall                                                               right wall
                  |                                                                       |
    Ni[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
                  |                                                                       |
                  |                                                                       |
    */

    double lC,lR,lL; //координаты центров €чеек
    double dlC,dlR,dlL; //разности граней [i]-й €чейки, между центрами €чеек слева ([i-1] и [i]) и справа ([i] и [i+1])

    int i;

    if(!strcmp(Geom,"radial"))
    {
        for(i=1;i<I+1;i++)
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
        for(i=1;i<I+1;i++)
        {
            lL = 0.5*(l[i-1]+l[i]);
            lC = 0.5*(l[i]+l[i+1]);
            lR = 0.5*(l[i+1]+l[i+2]);

            dlL = lC-lL;
            dlC = l[i+1]-l[i];
            dlR = lR-lC;

            GF_L[i] = 1.0/dlC*dlL);
            GF_R[i] = 1.0/dlC*dlR);
            GF_C[i] = -(GF_R[i]+GF_L[i]);
        }
    }
}
void Transport_SWEEPsolve(double Ni, int Gf)
{

    //—етка по длине:
	/*            left wall                                                               right wall
                  |                                                                       |
    Ni[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
                  |                                                                       |
                  |                                                                       |
    */

	if(Gf==0)
        Gf = Transport_GFcalc(char Geom);

    for(i=1;i<I+1;i++)
    {

    }


}
void Transport_boundary();
{

}
