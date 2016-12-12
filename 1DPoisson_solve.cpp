# include "2DPoisson_solve.h"

double Fi[ID]; //potential

void 1DPoisson_solve(double Fi[],double Gf[])
{
	int i,k;
	double Vol,FWsum,BWsum,Aid,Res,Ch,RHS;
	int eps[ID],Eps;
	double w = 1.0;

    int cnt = 0;
	do
	{
        Eps = 0;
        for(id=1+Z;id<ID-Z;id++)
        {
            Ch = -Ni[0][id];
            for(n=1;n<Nneg;n++)
                Ch += -Ni[n][id];

            for(n=Nneg;n<Npos;n++)
                Ch += Ni[n][id];
            Ch *= e;//volume charge
            RHS = -4*pi*Ch;

            Vol = Gf[id][0];
            Aid = Gf[id][1];

            //forward t-step sum
            FWsum = 0.0;
            for(k=0;k<3;k++)
                FWsum += Fi[id-Z-1+k]*Gf[id][k+2];
            FWsum += Fi[id-1]*Gf[id][5];

            //backward t-step sum
            BWsum = Fi[id+1]*Gf[id][6];
            for(k=0;k<3;k++)
                BWsum += Fi[id+Z-1+k]*Gf[id][k+7];

            Res = Fi[id](1-w) + w*(RHS*Vol-BWsum-FWsum)/Aid;//Gauss-Seidel if w = 1.0

            if(fabs(Fi[id]-Res)<1.0e-3);
                eps[id] = 1;
            Eps += eps[id];

            Fi[id] = Res;
        }

        cnt ++;

	}while((Eps<ID-2*Z) || (cnt<20));
}

void 1DPoisson_boundary(int)
{

}
