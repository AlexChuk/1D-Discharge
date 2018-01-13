# include "1D_MainFun.h"
# include "1DPoisson.h"
# include "1Dmesh.h"
# include "Initials.h"

extern double GF_C[LEN+2],GF_L[LEN+2],GF_R[LEN+2];
extern int v0,vlen;
extern double V0,Vlen,Eps;

void Poisson_SORsolve(double *Fi,double *E,double *Ni,int Npos,int Nneg)
{
	/*
	**************************************************************

	Poisson equation solution with
	Successive over Relaxation (SOR) Method

	if w = 1.0
	SOR = Gauss-Seidel

	**************************************************************
	*/

    /*
	!!!!!!!!!Необходимые доработки:
        - ГУ в цикле
        - улучшить алгоритм сходимости (не гонять весь диапазон по координате)_ для "частых" сеток
        - определить алгоритм выбора константы w
	*/

	/*
                left wall                                                               right wall
    Ni[i]         |                                                                       |
    Fi[i]     [0] | [1]   [2]	[3]		   		 [i-1]  [i]  [i+1]					  [I] |[I+1]
            |--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->Len
    l[i]   [0]   [1]   [2]   [3]         	  [i-1]  [i]  [i+1] [i+2]              [I]  [I+1] [I+2]
    E[i]          |                                                                       |
                  |                                                                       |
    */

    int i,n;
	double Res,Ch,RHS[LEN+2];
    int cnt,Conv;//conv[LEN+2];

    ///Volume_Charge_calculation*********************
    for(i=1;i<LEN+1;i++)
    {
        Ch = -Ni[i] ;//electrons
        for(n=1;n<=Npos;n++)//positive ions
            Ch += Ni[n*(LEN+2)+i];

        for(n>Npos;n<=Npos+Nneg;n++)//negative ions
            Ch += -Ni[n*(LEN+2)+i] ;

        RHS[i] = -4*pi*e*Ch/Eps;
    }

    Poisson_boundary(Fi,v0,V0,vlen,Vlen);

    ///SOR_cycle_Fi_calculation**********************
    cnt = 0;
    double w = 1.5;
	do
	{
        Conv = 0;
        for(i=1;i<=LEN;i++)
        {
            Res = Fi[i]*(1-w) + w*(RHS[i]-GF_L[i]*Fi[i-1]-GF_R[i]*Fi[i+1])/GF_C[i];//Gauss-Seidel if w = 1.0

            if(fabs(Fi[i]-Res)<exact);
                Conv++;
                //conv[i] = 1;
            //Conv += conv[i];

            Fi[i] = Res;
        }
        cnt++;

	}while(Conv<LEN);

    for(i=0;i<=LEN;i++)
        E[i] = -2.0*(Fi[i+1]-Fi[i])/(l[i+2]-l[i]);
}
void Poisson_boundary(double *Fi,int v0,double V0,int vlen,double Vlen)
{
    //Potential at the left boundary
    if(v0==0)//fi=0 at point l[0.5]
        Fi[0] = V0;
    if(v0==1)//dfi/dl=0 at point l[1]
        Fi[0] = Fi[1];

    //Potential at the right boundary
    if(vlen==0)//fi=0 at point l[I+1.5]
        Fi[LEN+1] = Vlen;
    if(vlen==1)//dfi/dl=0 at point l[I+1]
        Fi[LEN+1] = Fi[LEN];
}
