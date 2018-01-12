# include "1D_MainFun.h"
# include "1DPoisson.h"
# include "1DTransport.h"

double w = 1.5;
extern double GF_C[LEN+2],GF_L[LEN+2],GF_R[LEN+2];

void Poisson_SORsolve(double *Fi,double *Ni,int Npos,int Nneg)
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
	!!!!!!!!!����������� ���������:
        - �� � �����
        - �������� �������� ���������� (�� ������ ���� �������� �� ����������)_ ��� "������" �����
        - ���������� �������� ������ ��������� w
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
	double Res,Ch,RHS;

	/*if(Gf==0)
        Gf = Poisson_GFcalc(char Geom);*/

    //Poisson_boundary(v0,V0,vlen,Vlen);

    int cnt = 0,Conv,conv[LEN+2];

	do
	{
        Conv = 0;
        for(i=1;i<LEN+1;i++)
        {
            Ch = -Ni[i] ;//electrons
            for(n=1;n<Npos;n++)//positive ions
                Ch += Ni[n*(LEN+2)+i] ;

            for(n=Npos;n<Nneg;n++)//negative ions
                Ch += -Ni[n*(LEN+2)+i] ;

            RHS = -4*pi*e*Ch;

            Res = Fi[i]*(1-w) + w*(RHS-GF_L[i]*Fi[i-1]-GF_R[i]*Fi[i+1])/GF_C[i];//Gauss-Seidel if w = 1.0

            if(fabs(Fi[i]-Res)<exact);
                Conv ++;
                //conv[i] = 1;
            //Conv += conv[i];

            Fi[i] = Res;
        }

        cnt ++;

	}while(Conv<LEN);
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
