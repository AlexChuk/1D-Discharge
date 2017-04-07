# include "1D_MainFun.h"
# include "Initials.h"

double Ne[I+2][NEmax],Ni[Nmax][I+2],Mi[Nmax],LJi[Nmax][2],Roi[Nmax][I+2],Xi[Nmax][I+2],Pgas,Tgas[I+2],Ngas[I+2],Rogas[I+2],Hgas[I+2];
double Nel[I+2],Te[I+2],Tv[I+2];
double Gamma[Nmax][2];
double E[I+1],E_N[I+1];
double Tin,dTgas,dTe,dNel;
double Len,Tw,Lam;
double tau,dt,dte;
double Emax,dE,dEev;
int v0,vlen;
double V0,Vlen;

int N,NR,Nt,Nte,Ndots;
char Spec[Nmax][10],Spec_R[Nmax][10],GeomVect[10];
double CXi[Nmax][2][8];

void init_read()//���������� ��������� ������
{
    char symb[] = "------------------------------------------------------------";

	FILE *init;
	init = fopen("Comps_N2.txt", "r");

	char Cmt[300];

    //���������� �����������
	fscanf(init,"%s",&Cmt);
	do
	{
		fscanf(init,"%s",&Cmt);
	}while(strcmp(Cmt,symb)!=0);

	//���������� ��������� �������, �������� ��������� � �� ������������*************************
	fscanf(init,"%s",&Cmt);
	fscanf(init,"%d%s",&N,&Cmt);
	fscanf(init,"%lf%s",&tau,&Cmt);
	fscanf(init,"%lf%s",&dt,&Cmt);
	fscanf(init,"%lf%s",&dte,&Cmt);
	fscanf(init,"%d%s",&Ndots,&Cmt);
	fscanf(init,"%lf%s",&Pgas,&Cmt);
	fscanf(init,"%lf%s",&Tin,&Cmt);
	fscanf(init,"%lf%s",&E_N,&Cmt);
	fscanf(init,"%lf%s",&Emax,&Cmt);//���������� ��������� ��������� �����
	fscanf(init,"%lf%s",&Len,&Cmt);
	fscanf(init,"%lf%s",&Tw,&Cmt);

	fscanf(init,"%s%s",&Geom,&Cmt);
	fscanf(init,"%d%lf%s",&v0,&V0,&Cmt);
	fscanf(init,"%d%lf%s",&vlen,&Vlen,&Cmt);

	dEev = Emax/NEmax;
	dE = dEev*1.602e-12;//[eV]=1.602e-12[erg]

	int n,x,i,N1;
	fscanf(init,"%s",&Cmt);
	//fscanf(init,"%d%s%lf",&N1,&Spec[0],&Xi[0],&Gamma[0][0]);//���������
	for(n=0;n<N;n++)
        fscanf(init,"%d%s%lf%lf%lf",&N1,&Spec[n],&Xi[n],&Gamma[n][0],&Gamma[n][1]);

  	//����������_���������_���������_���_������_���������_���������*******************************
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
void init_data()//������� ��������� �������
{
	int i,n,k;

	gas_HCpSi_calc(Tgas);

	//������� ����������**********************************************
	for(i=0;i<=I+1;i++)
    {
        Pgas = Pgas*p0; //Torr
        Tgas[i] = Tin;
        Ni[0][i] = Xi[0][i];
        Ngas[i] = 0.0;
        Hgas[i] = 0.0;
        Rogas[i] = 0.0;
        for(n=1;n<N;n++)
        {
            Ni[n][i] = Xi[n]*Pgas[i]/(kb*Tgas[i]);
            Roi[n][i] = Ni[n][i]*Mi[n];

            Ngas[i] += Ni[n][i];
            Rogas[i] += Roi[n][i];
            Hgas[i] += HCpSi[0][n]*Roi[n][i];//[���/cm^3]
        }
    }
	//****************************************************************

	//��������� �������������� ����***********************************
	for(i=0;i<=I+1;i++)
        Fi[i] = 0.0;

	/*
	for(i=0;i<=I;i++)
    {
        E[i] = E_N*Ngas[i]*1e-17;//E_N[Td] = E[B/cm]*1e17/Ngas[cm-3];
        E[i] = E[i]/E0;//E[�/��] = E0*E[abs]
    }
    */
	//****************************************************************

	//��������� EEDF**************************************************
    for(i=0;i<=I+1;i++)
    {
        Te[i] = 2.0;
        Nel[i] = Xi[0];
        for(k=0;k<=NEmax-1;k++)
        {
            Ne[i][k] = 2*Nel*pow((k+0.5)*dEev/(Te*Te*Te*pi),0.5)*exp(-(k+0.5)*dEev/Te);//Maxwell
            if(Ne[i][k]<1e-30)
                Ne[i][k] = 0.0;
        }

    }
	//****************************************************************

	//Writing_data****************************************************
	init_print();
	//****************************************************************
}
void init_print()//������ ��������� ������
{
	int n,k;
	char symb[] = "------------------------------------------------------------";
	FILE *log;

	printf("%s\nInitial Data:\n\n",symb);

	//������ ����******************************************************************************

	printf("t = %.2lf[s]\n",0.0);

	printf("E/N = %.1f[Td]\tXe = %.2e\tTe = %.1f[eV]\n",E_N,Nel/Ngas,Te);

	log = fopen("EEDF_data.txt", "w");

    fprintf(log,"t,[s]\tTe,[eV]\t");
	for(k=0; k<NEmax; k++)
	{
	    fprintf(log,"%.2lf\t",dEev*(k+0.5));
		k += 10;
	}
	fprintf(log,"\n");

    fprintf(log,"t=%.1lf[s]\t %.2lf\t",0.0,Te);
	for(k=0; k<NEmax; k++)
	{
	    fprintf(log,"%.2e\t",Ne[k]/Nel);
		k += 10;
	}
	fprintf(log,"\n");

	/*
	fprintf(log,"t=0.0 s\nE,eV\t Ne(E),[cm-3/eV]\n");
	for(k=0; k<NEmax; k++)
	{
		fprintf(log,"%.2lf\t %.2e\n",dEev*(k+0.5),Ne[k]);
		k += 10;
	}
	fprintf(log,"\n\n");
	*/

	fclose(log);

	//������ ���������� ����*******************************************************************

	printf("P = %.1f[Torr]\tT = %.1f[K]\t[%s] = %.2e\n",Pgas/p0,Tgas,Spec[5],Ni[5]);

	log = fopen("Gas_data.txt", "w");

	fprintf(log,"t\t T,[K]\t P,[Torr]\t Hgas,[erg/cm3]\t Ngas,[cm-3]\t Ro,[g/cm3]\t Xe\t Te,[K]\t Tv,[K]\t X(N/N2)\t");
	for(n=0;n<N;n++)
		fprintf(log,"%s\t",Spec[n]);
	fprintf(log,"\n");

	fprintf(log,"%.2e\t%.1f\t%.1f\t %.2e\t %.2e\t %.2e\t %.2e\t %.1f\t %.1f\t %.2e\t",0.0,Tgas,Pgas/p0,Hgas,Ngas,Rogas,Nel/Ngas,Te*eV_K,Tv,Ni[8]/Ni[11]);
	for(n=0;n<N;n++)
		fprintf(log,"%.2e\t",Ni[n]);
	fprintf(log,"\n");

	fclose(log);

	printf("\n%s\n\n",symb);

	//������ VDF******************************************************************************
	log = fopen("VDF_data.txt", "w");

	fprintf(log,"t,[s]\t");
	for(n=11;n<N;n++)
		fprintf(log,"%d\t",n-11);
	fprintf(log,"\n");

	fclose(log);
}
void init_gasDBparser()//��������� ���������� ������� ��������� �� ���� ������
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
    for(n=1;n<N;n++)
    {
        //�������� ��������:������������-������������ ��� ���
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

        //�����_��_����_������**************************
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
                    fscanf(db,"%lf",&CXi[n][x][i]);//CXi[n][x][7]-���.��������� ����������� ��������� �� ��������� ��������� (�� ���. � ����������)
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

    //������_����**************************************************************
    FILE *log;
    log = fopen("Log_GasDB.txt", "w");

    fprintf(log,"Parsed elements' parameters from LogGasDB.txt\n\n");
    for(n=1;n<N;n++)
        fprintf(log,"%d\t%s\t%.3lf\t%.2lf\t%.2lf\n",n,Spec[n],Mi[n]/ma,LJi[n][0],CXi[n][0][7]);

    fclose(log);
    //*************************************************************************

    //������ ������ �� �������:
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
