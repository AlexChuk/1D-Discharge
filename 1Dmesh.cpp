# include "1D_MainFun.h"
# include "1Dmesh.h"

double l[LEN+3];
double GF_C[LEN+2],GF_L[LEN+2],GF_R[LEN+2];

void Mesh_calc(double Len)
{
    //—етка по длине:

	/*            left wall                                                               right wall
                  |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
                  |                                                                       |
                  |                                                                       |
    */


	int i;
	double dl;

	dl = Len/LEN;
	for(i=0;i<LEN+3;i++)
        l[i] = i*dl;
}
void Mesh_GFcalc(char *Geom)
{
    //—етка по длине:
	/*           left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    int i;
    double lC,lR,lL; //координаты центров €чеек
    double dlC,dlR,dlL; //разности граней [i]-й €чейки, между центрами €чеек слева ([i-1] и [i]) и справа ([i] и [i+1])

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

