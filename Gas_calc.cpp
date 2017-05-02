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
	double Xi[N],Nin,Roin,Tin,Pin,Hin,XMi,Xmi[N];
	double Nout,Pout,Tout,Roout,Hout;
    double Hi,Cpi;
	double ftn,Ftn,Tnn,Tn,dT,Nn;

	*dNel = *Nel;

	//Присвоение значений из адресов переменных
	Tin = *Tgas;
	Pin = *Pgas;
	Hin = *Hgas;

	Nin = 0.0;
	Roin = 0.0;
	for(n=0;n<N;n++)
		Nin += *(Ni+n*(LEN+2));//Nin += Ni[n*(LEN+2)];

	XMi = 0.0;
    for(n=0;n<N;n++)
    {
        Xi[n] = *(Ni+n*(LEN+2))/Nin;//Xi[n] = Ni[n*(LEN+2)]/Nin;

        Xmi[n] = Xi[n]*Mi[n];
        XMi += Xmi[n];
    }
    Roin = XMi*Nin;

    //Расчёт температуры_методом Ньютона
	Tn = Tin;//300;
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
		gas_HCpSi_calc(Tn,N);
		for(n=1;n<N;n++)//Npos+Nneg
		{
			Hi = HCpSi[0][n];//[эрг/г]
			Cpi = HCpSi[1][n];//[эрг/г*К]

			Ftn += Xmi[n]*Hi;//Roi1[n]*Hi;//Xi[n]*Hi;//
			ftn += Xmi[n]*Cpi;//Roi1[n]*Cpi;//Xi[n]*Cpi;//
		}

		Ftn = Ftn*Nn - Hin;

		Tnn = Tn - Ftn/(ftn*Nn);

		dT = fabs(Tnn-Tn);

	}while(dT>0.01);
	Tout = Tnn;//Tin//300;

	//Isobaric process**************************************
	Pout = Pin;

	//Return to using variables:
	Nout = Pout/(kb*Tout);
    Roout = XMi*Nout;

	gas_HCpSi_calc(Tout,N);
	Hout = 0;
	for(n=1;n<N;n++)
		Hout += Xmi[n]*HCpSi[0][n];//[эрг/cm^3]
    Hout *= Nout;

	//Return to using variables:
	*Pgas = Pout;
	*Tgas = Tout;
	*Ngas = Nout;
	*Rogas = Roout;
	*Hgas = Hout;
    *dTgas = fabs(Tin-Tout);

	for(n=0;n<N;n++)
        *(Ni+n*(LEN+2)) = Nout*Xi[n];//Ni[n*(LEN+2)] = Nout*Xi[n];

    //концентрация электронов
    *Nel = Nout*Xi[0];
    *dNel = fabs(*dNel-*Nel);

    //Tv-calculation
    /*for(n=0;n<N;n++)
    {
        if(!strcmp(Spec[n],"N2(0)"))
            break;
    }*/
    n=11;
    *Tv = fabs(HCpSi[0][n+1]-HCpSi[0][n])*Mi[n]/kb/log((*(Ni+n*(LEN+2)))/(*(Ni+(n+1)*(LEN+2))));//*Tv = fabs(HCpSi[0][n+1]-HCpSi[0][n])*Mi[n]/kb/log(Ni[n*(LEN+2)]/Ni[(n+1)*(LEN+2)]);

	//************************************************************
}
void gas_LenPrint(double *Ni,int N,double *Pgas,double *Tgas,double *Ngas,double *Hgas,double *Rogas,double *Nel,double *Te,double *Tv,double *E,double *Ne,double tic)//запись в файл
{
	int i,k,n;
	FILE *log;

	//Вывод на экран
	i = 1;
	printf("Time = %.2e[s]\n",tic);
	printf("Point - l[%d] = %.2lf\n",i,l[i]);
	printf("P = %.1lf[Torr]\tT = %.1lf[K]\tXe = %.2e\t[%s] = %.2e[cm-3]\n",Pgas[i]/p0,Tgas[i],Nel[i]/Ngas[i],Spec[5],Ni[5*(LEN+2)+i]);

    //запись параметров газа*******************************************************************
	log = fopen("Gas_data.txt", "a+");

	fprintf(log,"Time = %.2e[s]\n",tic);

    //Data_print
    fprintf(log,"Pgas,Torr\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.1lf\t",Pgas[i]/p0);
    fprintf(log,"\n");

    fprintf(log,"Tgas,K\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.1lf\t",Tgas[i]);
    fprintf(log,"\n");

    fprintf(log,"Tv,K\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.2lf\t",Tv[i]);
    fprintf(log,"\n");

    fprintf(log,"Ngas,cm-3\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.2e\t",Ngas[i]);
    fprintf(log,"\n");

    fprintf(log,"Rogas,g/cm^3\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.2e\t",Rogas[i]);
    fprintf(log,"\n");

    fprintf(log,"Hgas,erg/cm^3\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.2e\t",Hgas[i]);
    fprintf(log,"\n");

    //Ni[n]
    for(n=0;n<=21;n++)
    {
        fprintf(log,"%s,cm-3\t",Spec[n]);
        for(i=0;i<LEN+2;i++)
            fprintf(log,"%.2e\t",Ni[n*(LEN+2)+i]);//(LEN+2) - due to 2 additional virtual points
        fprintf(log,"\n");
    }

    //electron_data
    fprintf(log,"Xel\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.2e\t",Nel[i]/Ngas[i]);
    fprintf(log,"\n");

    fprintf(log,"E/N,Td\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.1lf\t",E[i]*1.0e17*E0/Ngas[i]);
    fprintf(log,"\n");

    fprintf(log,"Te,eV\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.2lf\t",Te[i]);
    fprintf(log,"\n");

    fprintf(log,"Te,K\t");
    for(i=0;i<LEN+2;i++)
        fprintf(log,"%.2lf\t",Te[i]*eV_K);
    fprintf(log,"\n");

    //VDF_print
    int I[] = {1,int(LEN/2),LEN};
    for(i=0;i<3;i++)
    {
        fprintf(log,"VDF(l=%.2lfcm)\t",l[I[i]]);
        for(n=11;n<N;n++)
            fprintf(log,"%.2e\t",Ni[n*(LEN+2)+I[i]]/Ni[11*(LEN+2)+I[i]]);
        fprintf(log,"\n");
    }
    fprintf(log,"\n");

    //EEDF_print
    fprintf(log,"E,eV\t");
	for(k=0;k<NEmax;k+=10)
        fprintf(log,"%.2lf\t",dEev*(k+0.5));
    fprintf(log,"\n");
    for(i=0;i<3;i++)
    {
        fprintf(log,"EEDF(l=%.2lfcm)\t",l[I[i]]);
        for(k=0; k<NEmax; k+=10)
            fprintf(log,"%.2e\t",Ne[NEmax*I[i]+k]/Nel[I[i]]);
        fprintf(log,"\n");
    }

    fprintf(log,"\n\n");

	fclose(log);
}
void gas_TimePrint(double *Ni,int N,double Pgas,double Tgas,double Ngas,double Nel,double Te,double Tv,double E,double tic)//запись в файл
{
    int i,k,n;
	FILE *log;

	if(tic==0.0)
    {
        log = fopen("Gas_TimeData.txt", "w");

        fprintf(log,"t,s\tPgas,Torr\tTgas,K\tTv,[K]\tTe,[K]\tE/N,Td\tNgas,[cm-3]\tXel\t");
        for(n=0;n<=11;n++)
            fprintf(log,"%s\t",Spec[n]);
        fprintf(log,"\n");
    }
    else
        log = fopen("Gas_TimeData.txt", "a+");


    fprintf(log,"%.2e\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.2e\t%.2e\t",tic,Pgas/p0,Tgas,Tv,Te*eV_K,E*1.0e17*E0/Ngas,Ngas,Nel/Ngas);
    for(n=0;n<=11;n++)
        fprintf(log,"%.2e\t",*(Ni+n*(LEN+2)));
    fprintf(log,"\n");

	fclose(log);
}

