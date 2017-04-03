# include "1D_MainFun.h"
# include "1DPoisson.h"

double Fi[I+2]; //potential
double GF_C[I+2],GF_L[I+2],GF_R[I+2];
int Gf = 0;
double w = 1.5;

int 1DPoisson_GFcalc(char Geom)
{
    //Сетка по длине:
	/*           left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    double lC,lR,lL; //координаты центров ячеек
    double dlC,dlR,dlL; //разности граней [i]-й ячейки, между центрами ячеек слева ([i-1] и [i]) и справа ([i] и [i+1])

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

    return Gf = 1;
}
void 1DPoisson_SORsolve(double Fi[])
{
	/*
	**************************************************************

	Poisson equation solution with
	Successive over Relaxation (SOR) Method

	if w = 1.0
	SOR = Gauss-Seidel

	**************************************************************
	*/

	int i;
	double Res,Ch,RHS;

	if(Gf==0)
        Gf = 1DPoisson_GFcalc(char Geom);

    1DPoisson_boundary(v0,V0,vlen,Vlen);

    int cnt = 0,int Conv,conv[I+2];

	/*
	Необходимые доработки:
        - ГУ в цикле
        - улучшить алгоритм сходимости (не гонять весь диапазон на координате)
        - определить алгоритм выбора константы w
	*/

	do
	{
        Conv = 0;
        for(i=1;i<I+1;i++)
        {
            Ch = -Ni[0][i];//electrons
            for(n=1;n<Nneg;n++)//negative ions
                Ch += -Ni[n][i];

            for(n=Nneg;n<Npos;n++)//positive ions
                Ch += Ni[n][i];

            RHS = -4*pi*e*Ch;

            Res = Fi[i](1-w) + w*(RHS-GF_L[i]*Fi[i-1]-GF_R[i]*Fi[i+1])/GF_C[i];//Gauss-Seidel if w = 1.0

            if(fabs(Fi[i]-Res)<exact);
                Conv ++;
                //conv[i] = 1;
            //Conv += conv[i];

            Fi[i] = Res;
        }

        cnt ++;

	}while(Conv<I);
}
void 1DPoisson_boundary(int v0,double V0,int vlen,double Vlen)
{
    //Potential at the left boundary
    if(v0==0)//fi=0 at point l[0.5]
        Fi[0] = V0;
    if(v0==1)//dfi/dl=0 at point l[1]
        Fi[0] = Fi[1];

    //Potential at the right boundary
    if(vlen==0)//fi=0 at point l[I+1.5]
        Fi[I+1] = Vlen;
    if(vlen==1)//dfi/dl=0 at point l[I+1]
        Fi[I+1] = Fi[I];
}
