# include "1D_MainFun.h"
# include "Initials.h"

double Ne[LEN+2][NEmax],Ni[Nmax][LEN+2],Mi[Nmax],LJi[Nmax][2],Pgas[LEN+2],Tgas[LEN+2],Ngas[LEN+2],Rogas[LEN+2],Hgas[LEN+2];
double Nel[LEN+2],Te[LEN+2],Tv[LEN+2];
double Gamma[Nmax][2];
double Ez,Er[LEN+1],Fi[LEN+2],Iexp;
double Tinit,Pinit,Xinit[Nmax],E_Ninit;
double dTgas,dTe,dNel;
double Len,Hght,Tw;
double tau,dt;
double Emax,dE,dEev;
int v0,vlen;
double V0,Vlen;

int N,NR,Nt,Nte,Ndots;
char Spec[Nmax][10],Spec_R[Nmax][10],GeomVect[10];
char Geom[20];
double CXi[Nmax][2][8];

void init_read()//считывание начальных данных
{
    char symb[] = "------------------------------------------------------------";

	FILE *init;
	init = fopen("Comps_N2.txt", "r");

	char Cmt[300];

    //считывание комментария
	fscanf(init,"%s",&Cmt);
	do
	{
		fscanf(init,"%s",&Cmt);
	}while(strcmp(Cmt,symb)!=0);

	//считывание начальных условий, названий компонент и их концентрации*************************
	fscanf(init,"%s",&Cmt);
	fscanf(init,"%d%s",&N,&Cmt);
	fscanf(init,"%lf%s",&tau,&Cmt);
	fscanf(init,"%lf%s",&dt,&Cmt);
	fscanf(init,"%d%s",&Ndots,&Cmt);
	fscanf(init,"%lf%s",&Pinit,&Cmt);
	fscanf(init,"%lf%s",&Tinit,&Cmt);
	fscanf(init,"%lf%s",&E_Ninit,&Cmt);
	fscanf(init,"%lf%s",&Iexp,&Cmt);
	fscanf(init,"%lf%s",&Emax,&Cmt);//считывание максимума энергетич шкалы
	fscanf(init,"%lf%s",&Len,&Cmt);
	fscanf(init,"%lf%s",&Hght,&Cmt);
	fscanf(init,"%lf%s",&Tw,&Cmt);

	fscanf(init,"%s%s",&Geom,&Cmt);
	fscanf(init,"%d%lf%s",&v0,&V0,&Cmt);
	fscanf(init,"%d%lf%s",&vlen,&Vlen,&Cmt);

	dEev = Emax/NEmax;
	dE = dEev*1.602e-12;//[eV]=1.602e-12[erg]

	int n,x,i,N1;
	fscanf(init,"%s",&Cmt);
	//fscanf(init,"%d%s%lf",&N1,&Spec[0],&Xi[0],&Gamma[0][0]);//электроны
	for(n=0;n<N;n++)
        fscanf(init,"%d%s%lf%lf%lf",&N1,&Spec[n],&Xinit[n],&Gamma[n][0],&Gamma[n][1]);

  	//Считывание_выбранных_компонент_для_вывода_скоростей_процессов*******************************
	fscanf(init,"%s",&Cmt);
	NR = 0;
	int err;
	fscanf(init,"%s",&Spec_R[NR]);
	while(strcmp(Spec_R[NR],symb)!=0)
	{
		err = 1;
		for(n=0;n<N;n++)
		{
			if(!strcmp(Spec[n],Spec_R[NR]))
			{
				NR += 1;
				err = 0;
				break;
			}
		}

		if(err == 1)
			printf("!!Warning:Unknown component - %s for ChemRates analysis!!Excluded!\n\n",Spec_R[NR]);

		fscanf(init,"%s",&Spec_R[NR]);
	}
	fclose(init);
}
void init_data()//задание начальных условий
{
	int i,n,k;

	//gas_HCpSi_calc(Tinit,N);

	//Газовые компоненты**********************************************
	for(i=0;i<=LEN+1;i++)
    {
        Pgas[i] = Pinit*p0; //Torr
        //Tgas[i] = Tinit;
        Tgas[i] = Tinit-i*i*(Tinit-Tw)/(LEN+1)/(LEN+1);

        gas_HCpSi_calc(Tgas[i],N);

        Ngas[i] = Pgas[i]/(kb*Tgas[i]);
        Hgas[i] = 0.0;
        Rogas[i] = 0.0;
        for(n=0;n<N;n++)
        {
            Ni[n][i] = Xinit[n]*Ngas[i];
            Rogas[i] += Ni[n][i]*Mi[n];
            if(n>=1)
                Hgas[i] += HCpSi[0][n]*Ni[n][i]*Mi[n];//[эрг/cm^3]
        }
        Nel[i] = Ni[0][i];
    }
	//****************************************************************

	//Потенциал электрического поля***********************************

	Iexp = Iexp;//[Amper]
    Ez = E_Ninit*Ngas[0]*1e-17;//E_N[Td] = E[B/cm]*1e17/Ngas[cm-3];
    Ez = Ez/Eabs;//E=E[B/cm]/E[abs]

	for(i=0;i<=LEN;i++)
    {
        //E[i] = E_Ninit*Ngas[i]*1e-17;//E_N[Td] = E[B/cm]*1e17/Ngas[cm-3];
        Er[i] = 0.0;//E_Ninit*Ngas[0]*1e-17;//E_N[Td] = E[B/cm]*1e17/Ngas[cm-3];
        Er[i] = Er[i]/Eabs;//E=E[B/cm]/E[abs]
    }

    /*for(i=0;i<=LEN+1;i++)
        Fi[i] = 0.0;*/

	//****************************************************************

	//Начальная EEDF**************************************************
    for(i=0;i<=LEN+1;i++)
    {
        Te[i] = 0.05;
        for(k=0;k<=NEmax-1;k++)
        {
            Ne[i][k] = 2*Nel[i]*pow((k+0.5)*dEev/(Te[i]*Te[i]*Te[i]*pi),0.5)*exp(-(k+0.5)*dEev/Te[i]);//Maxwell
            if(Ne[i][k]<1e-30)
                Ne[i][k] = 0.0;
        }

    }
	//****************************************************************

	//Writing_data****************************************************
	init_print();
	//****************************************************************
}
void init_print()//запись начальных данных
{
	int n,k;
	char symb[] = "------------------------------------------------------------";

	printf("%s\nInitial Data:\n\n",symb);
	printf("t = %.2lf[s]\n",0.0);
	printf("I = %.2lf[A]\tE/N = %.1f[Td]\tXe = %.2e\tTe = %.1f[eV]\n",Iexp,Ez/Ngas[1],Nel[1]/Ngas[1],Te[1]);
	printf("P = %.1f[Torr]\tT = %.1f[K]\t%s = %.2e\n",Pgas[1]/p0,Tgas[1],Spec[5],Ni[5][1]);
	printf("\n%s\n\n",symb);

	//File_rewriting**************************************************

    FILE *log;
    log = fopen("Gas_data.txt", "w");
	fclose(log);
}
void init_gasDBparser()//получение параметров газовых компонент из базы данных
{
    char symb[] = "------------------------------------------------------------";
    char Cmt[300];
    char str[10],strv[10],VCmt[100];
    double we,xewe,dH0,Enull;

    FILE *db;
	db = fopen("GasPropDB.txt", "r");

    int n,x,i,Nerr,err;
    int nvib = 0;
    err = 0;
    Nerr = 0;
    Mi[0] = me;//electrons
    for(n=1;n<N;n++)
    {
        //проверка элемента:колебательно-возбужденное или нет
        err = 0;
        i = 0;
        while((Spec[n][i]!='(') && (Spec[n][i]!='\0'))
            i++;

        if(Spec[n][i]=='(')
        {
            memset(str, 0, sizeof(str));
            int k = 0;
            while(Spec[n][k]!='(')
            {
                str[k] = Spec[n][k];
                k++;
            }
            nvib = 1;

            memset(VCmt, 0, sizeof(VCmt));
            sprintf(VCmt,"//%s_Spectrum_data:",str);
        }
        else
        {
            memset(str, 0, sizeof(str));
            strcat(str,Spec[n]);
        }

        //поиск_по_базе_данных**************************
        while(strcmp(Cmt,str)!=0)
        {
            fscanf(db,"%s",&Cmt);

            if(!strcmp(Cmt,"END"))
            {
                err++;
                break;
            }
        }

        if(err==0)
        {
            fscanf(db,"%lf%lf%lf",&Mi[n],&LJi[n][0],&LJi[n][1]);
            Mi[n] = Mi[n]*ma;
            for(x=0;x<2;x++)
            {
                for(i=0;i<8;i++)
                    fscanf(db,"%lf",&CXi[n][x][i]);//CXi[n][x][7]-доп.энтальпия образования компонент от основного состояния (не учт. в разложении)
            }

            if(nvib>0)
            {
                while(strcmp(Cmt,VCmt)!=0)
                    fscanf(db,"%s",&Cmt);
                fscanf(db,"%lf%lf",&we,&xewe);

                Enull = (0.5*we-0.5*0.5*xewe)/cm_eV;

                nvib -= 1;
                memset(strv, 0, sizeof(strv));
                sprintf(strv,"%s(%d)",str,nvib);
                while(!strcmp(strv,Spec[n+nvib]))
                {
                    dH0 = ((nvib+0.5)*we-(nvib+0.5)*(nvib+0.5)*xewe)/cm_eV-Enull;

                    Mi[n+nvib] = Mi[n];
                    LJi[n+nvib][0] = LJi[n][0];
                    LJi[n+nvib][1] = LJi[n][1];
                    for(x=0;x<2;x++)
                    {
                        for(i=0;i<8;i++)
                            CXi[n+nvib][x][i] = CXi[n][x][i];
                        if(nvib>0)
                            CXi[n+nvib][x][7] += dH0;
                    }
                    nvib ++;

                    sprintf(strv,"%s(%d)",str,nvib);
                }

                n += nvib-1;
                nvib = 0;
            }

            rewind(db);
        }
        else
        {
            printf("ERROR!!!_Element %s was found in DataBase!!!\n", Spec[n]);
            Nerr ++;

            rewind(db);
        }

    }

    fclose(db);

    //Запись_лога**************************************************************
    FILE *log;
    log = fopen("Log_GasDB.txt", "w");

    fprintf(log,"Parsed elements' parameters from LogGasDB.txt\n\n");
    for(n=1;n<N;n++)
        fprintf(log,"%d\t%s\t%.3lf\t%.2lf\t%.2lf\n",n,Spec[n],Mi[n]/ma,LJi[n][0],CXi[n][0][7]);

    fclose(log);
    //*************************************************************************

    //выдача отчета об ошибках:
    if(Nerr==0)
        printf("\nAll Components were parsed from DataBase\n");
    else
    {
        char a;
        printf("\nPlease, fill the DataBase!!Nerr=%d\n",Nerr);
        scanf("%c",&a);
        exit(1);
    }
}
