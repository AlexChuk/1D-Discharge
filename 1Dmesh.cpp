# include "1D_MainFun.h"
# include "1Dmesh.h"

double l[I+3];

void 1Dmesh_calc(double Len)
{
    //����� �� �����:
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

	dl = Len/I;
	for(i=0;i<I+3;i++)
        l[i] = i*dl;

}