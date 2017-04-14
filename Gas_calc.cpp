# include "1D_MainFun.h"
# include "Gas_calc.h"

double HCpSi[3][Nmax];

void gas_HCpSi_calc(double Tgas,int N)//нахождение внутренней энергии,энтальпии и теплоёмкостей по зад. температуре из аппроксиммаций Cp(T)
{
	double Hi,Cpi,Si,Mm,Cpt;
	int l,n,x;

	//T0 = Temp0;
	for(n=1;n<N;n++)
	{
		Mm = Mi[n];//[г]
		//Расчёт характеристик в точке T0
		if(Tgas>1000)
			x=0;
		else
			x=1;

		Cpi = 0;
		Hi = 0;
		Si = 0;
		for(l=0;l<5;l++)
		{
			Cpt = CXi[n][x][l]*pow(Tgas,l);

			Cpi += Cpt;//безразм
			Hi += Cpt/(l+1);//безразм

			if(l!=0)
				Si += Cpt/l;//безразм
		}
		Cpi = Cpi*kb/Mm;//[эрг/г*К]
		Hi = (Hi*kb*Tgas+CXi[n][x][5]*kb+CXi[n][x][7]*1.6e-12)/Mm;//[эрг/г]
		Si = (Si+CXi[n][x][0]*log(Tgas)+CXi[n][x][6])*kb/Mm;//[эрг/г*К]

		HCpSi[0][n] = Hi;//[эрг/г]
		HCpSi[1][n] = Cpi;//[эрг/г*К]
		HCpSi[2][n] = Si;//[эрг/г*К]
	}
}
void gas_TP_calc(double *Ni,int N,double *Pgas,double *Tgas,double *dTgas,double *Ngas,double *Rogas,double *Hgas,double *Nel,double *dNel,double *Tv)//расчёт температуры газа((H,P)-const)
{
	int n;
	double Xi[N],Nin,Roin,Tin,Pin,Hin;
	double Nout,Pout,Tout,Roout,Hout;
    double Hi,Cpi;
	double ftn,Ftn,Tnn,Tn,dT,Nn;

	dTgas = Tgas;
	dNel = Nel;

	//Присвоение значений из адресов переменных
	Tin = *Tgas;
	Pin = *Pgas;
	Hin = *Hgas;
	Nin = 0.0;
	Roin = 0.0;
	for(n=0;n<N;n++)
	{
		Nin += *Ni[n*LEN];
		Roin += *Ni[n*LEN]*Mi[n];
	}

    for(n=0;n<N;n++)
        Xi[n] = *Ni[n*LEN]/Nin;

    //Расчёт температуры_методом Ньютона
	Tn = Tin;//300;//Temp0;//Tin;//
	Nn = Nin;
	Tnn = Tin;
	do
	{
		if(!Tnn==0)
		{
            Tn = Tnn;
			Nn = Nin*Tin/Tn;//Pin/(kb*Tn);
		}

		Ftn = 0;
		ftn = 0;
		gas_HCpSi_calc(Tn);
		for(n=Npos+Nneg;n<N;n++)
		{
			Hi = HCpSi[0][n];//[эрг/г]
			Cpi = HCpSi[1][n];//[эрг/г*К]

			Ftn += Xi[n]*Nn*Mi[n]*Hi;//Roi1[n]*Hi;//Xi[n]*Hi;//
			ftn += Xi[n]*Nn*Mi[n]*Cpi;//Roi1[n]*Cpi;//Xi[n]*Cpi;//
		}

		Ftn = Ftn - Hin;

		Tnn = Tn - Ftn/ftn;

		dT = fabs(Tnn-Tn);

	}while(dT>0.01);
	Tout = Tnn;//Tin//300;

	//Isobaric process**************************************
	Pout = Pin;

	//Return to using variables:
	Nout = Pout/(kb*Tout);

	Roout = 0.0;
	for(n=0;n<N;n++)
		Roout += Nout*Xi[n]*Mi[n];

	gas_HCpSi_calc(Tout);
	Hout = 0;
	for(n=1;n<N;n++)
		Hout += Ni[n]*Mi[n]*HCpSi[0][n];

	//Return to using variables:
	*Pgas = Pout;
	*Tgas = Tout;
	*Ngas = Nout;
	*Rogas = Roout;
	*Hgas = Hout;
    *dTgas = fabs(*dTgas-Tout);

	for(n=0;n<N;n++)
        Ni[n*LEN] = Nout*Xi[n];

    *Nel = Nout*Xi[0];//концентрация электронов
    *dNel = fabs(*dNel-*Nel);

    //Tv-calculation
    for(n=0;n<N;n++)
    {
        if(!strcmp(Spec[n],"N2(0)"))
            break;
    }
    *Tv = fabs(HCpSi[0][n+1]-HCpSi[0][n])*Mi[n]/kb/log(Ni[n*LEN]/Ni[(n+1)*LEN]);

	//************************************************************
}
void gas_print(double *Ni,int N,double *Nel,double *Tgas,double Pgas,double *Hgas,double *Rogas,double *Te,double *Tv,double tic)//запись в файл
{
	int i,n;
	FILE *log;

	//запись параметров газа*******************************************************************
	log = fopen("Gas_data.txt", "a+");

	i = 1;
	printf("Time = %.2e[s]\n",tic);
	printf("Point - r[i] = %d\n",i);
	printf("P = %.1f[Torr]\tT = %.1f[K]\tXe = %.2e\t[%s] = %.2e\n",Pgas/p0,Tgas[i],Nel[i]/Ngas[i],Spec[5],Ni[5][i]);

	fprintf(log,"Time = %.2e[s]\n",tic);
	for(i=0;i<LEN+2;i++)
    {
        fprintf(log,"%.1f\t%.1f\t %.2e\t %.2e\t %.2e\t %.2e\t %.1f\t %.1f\t",Tgas[i],Pgas/p0,Hgas[i],Ngas[i],Rogas[i],Nel[i]/Ngas[i],Te[i]*eV_K,Tv);
        for(n=0;n<=11;n++)//N
            fprintf(log,"%.2e\t",Ni[n][i]);
        fprintf(log,"\n");
    }
    fprintf(log,"\n");

	fclose(log);

	//запись VDF********************************************************************************
	log = fopen("VDF_data.txt", "a+");

	fprintf(log,"Time = %.2e[s]\n",tic);
	for(n=11;n<N;n++)
		fprintf(log,"%.2e\t",Ni[n][i]/Ngas[i]);
	fprintf(log,"\n");

	fclose(log);
}
