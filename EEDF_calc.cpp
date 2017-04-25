# include "1D_MainFun.h"
# include "EEDF_calc.h"

extern double dte;

int CStype;
int CSref[Nmax][9][CSmax],CSR[CSmax][5];
double CS[NEmax][CSmax],Ith[CSmax];
double Ee,Vdr,Muel,Jel,Qel,QE;

int EEDF_read_CS(int N)//считывание сечений EEDF-процессов(возврат кол-ва реакций)
{
	FILE *cross;
	cross = fopen("Cross_N2_lxtype.txt", "r");

	FILE *log;
	log = fopen("Log_CS.txt", "w");
	fclose(log);

	log = fopen("Log_CS.txt", "a+");

	double X[50],Y[50],A[50],B[50];
	char CSstr[100],CSstr1[10],Cmt[100];

	CStype = 9;
	char KeyW[][15] =
	{
		"ELASTIC",
		"ROT-EXCITATION",
		"VIB-EXCITATION",
		"EXCITATION",
		"DISSOCIATION",
		"DEEXCITATION",
		"IONIZATION",
		"ATTACHMENT",
		"RECOMBINATION"
	};

	char symb[] = "------------------------------------------------------------";

	int i,k,K,j,J,err,n,n_t,Kw,add,d,s,Serr;

	//считывание комментария
	fscanf(cross,"%s",&Cmt);
	do
	{
		fscanf(cross,"%s",&Cmt);
	}while(strcmp(Cmt,symb)!=0);

	j = 0;
	Serr = 0;
	fscanf(cross,"%s",&CSstr);//считывание первой строки
	while(strcmp(CSstr,"END")!=0)
	{
		add = 0;

		err = 0;
		for(i=0;i<CStype;i++)//сравнение с ключевыми словами
		{
			if(!strcmp(CSstr,KeyW[i]))
			{
				Kw = i;
				CSR[j][0] = Kw;//CS-type

				err = 1;
				break;
			}
		}
		if(err==0)//нет совпадений
		{
			printf("!!!Can't find Keyword!!! %s\n", CSstr);
			fprintf(log,"!!!Can't find Keyword!!! %s\n", CSstr);

			Serr += 1;
			break;
		}

		//считывание второй строки
		fscanf(cross,"%s",&CSstr);
		err = 0;
		for(n=0;n<N;n++)
		{
			if(!strcmp(CSstr,Spec[n]))
			{
				n_t = n;//target particle
				CSR[j][2] = n_t+1;
				err += 1;
				break;
			}
		}
		if(err==0)
		{
			printf("!!!Unknown target particle - %s in #%d cross section!!!\n ",CSstr,j);
			fprintf(log,"!!!Unknown target particle - %s in #%d cross section!!!\n ",CSstr,j);

			Serr += 1;
			break;
		}

		//обработка 2-й строки
		if(Kw>=2)//w/o ELASTIC+ROT-EXCITATION
		{
			fscanf(cross,"%s%s",&Cmt,&CSstr);
			if(Kw==4)
				fscanf(cross,"%s",&CSstr1);

			if(!strcmp(Cmt,"="))//forvard or rev
				CSR[j][1] = -1;
			if(!strcmp(Cmt,"->"))
				CSR[j][1] = 1;

			err = 0;
			for(n=0;n<N;n++)
			{
				if(!strcmp(CSstr,Spec[n]))
				{
					CSR[j][3] = n+1;//product particle
					err += 1;

					if(err==2)
						break;
				}
				if((!strcmp(CSstr1,Spec[n]))&&(Kw==4))
				{
					CSR[j][4] = n+1;//product particle
					err += 1;

					if(err==2)
						break;
				}
			}

			//предупреждение
			if((Kw==4)&&(err<2))
			{
				printf("!!!Unknown products in CS: %s %s %s %s - !!!\n ",Spec[n_t],Cmt,CSstr,CSstr1);
				fprintf(log,"!!!Unknown products in CS: %s %s %s %s - !!!\n ",Spec[n_t],Cmt,CSstr,CSstr1);

				Serr += 1;
				break;
			}
			else if((err==0)&&(Kw!=4))
			{
				printf("!!!Unknown product particle in CS: %s %s %s - !!!\n ",Spec[n_t],Cmt,CSstr);
				fprintf(log,"!!!Unknown product particle in CS: %s %s %s - !!!\n ",Spec[n_t],Cmt,CSstr);

				Serr += 1;
				break;
			}
		}

		//обработка дополнения ко 2-й строке и считывание 3-ей строки
		fscanf(cross,"%s",CSstr);
		if(!strcmp(CSstr,"//add:"))//доп. реакции с таким же сечением "//add: N2(3) N2(5) //"
		{
			add = 0;
			fscanf(cross,"%s",CSstr);
			while(strcmp(CSstr,"//")!=0)
			{
				err = 0;
				for(n=0;n<N;n++)
				{
					if(!strcmp(CSstr,Spec[n]))
					{
						add += 1;
						CSR[j+add][0] = CSR[j][0];//type
						CSR[j+add][1] = CSR[j][1];//forvard or rev
						CSR[j+add][2] = CSR[j][2];//target particle
						CSR[j+add][3] = n+1;//product particle
						err += 1;
						break;
					}
				}

				if(err==0)
				{
					printf("!!!Unknown product particle in CS: %s -> add: %s - !!!\n ",Spec[n_t],CSstr);
					fprintf(log,"!!!Unknown product particle in CS: %s -> add: %s - !!!\n ",Spec[n_t],CSstr);

					Serr += 1;
					break;
				}

				fscanf(cross,"%s",CSstr);
			}

			fscanf(cross,"%lf%s",&Ith[j],Cmt);

			for(d=1;d<=add;d++)
			{
				n = CSR[j+d][3]-1;
				Ith[j+d] = Ith[j] + CXi[n][0][7];//написано только для случая постоянного реагента(таргета)!!!!!!!!!!внести правки!!!!!!
			}
		}
		else
		{
			Ith[j] = atof(CSstr);
			if(Kw==1)//rotational
				Ith[j] = Ith[j]/11605;//[eV]

			fscanf(cross,"%s",Cmt);
		}

		//формирование массива ссылок на номер сечения
		for(d=0;d<=add;d++)
		{
			s = CSref[n_t][Kw][0]+1;
			CSref[n_t][Kw][s] = j+1;//j+1??????????????????****************
			CSref[n_t][Kw][0] += 1;
		}

		//считывание 4-ей строки
		fscanf(cross,"%s",Cmt);

		//считывание и обработка сечения
		fscanf(cross,"%s",CSstr);
		if(!strcmp(CSstr,symb))
		{
			k = 0;
			//X[0] = 0.0;
			//Y[0] = 0.0;
			fscanf(cross,"%s",CSstr);
			while(strcmp(CSstr,symb)!=0)
			{
				X[k] = atof(CSstr);
				fscanf(cross,"%lf",&Y[k]);
				k += 1;
				fscanf(cross,"%s",CSstr);
			}
			K = k; //счётчик
			X[K] = 1500;//Emax;
			Y[K] = 0.0;//

			///Loging*****************************************************
			fprintf(log,"#CS=%d %s for particle %s\n",j,KeyW[Kw],Spec[n_t]);
			for(k=0;k<K;k++)
				fprintf(log,"%.2lf\t%.2e\n",X[k],Y[k]);
			fprintf(log,"\n");
			///Loging*****************************************************
		}
		else
		{
			printf("WARNING!Error while reading CS %s\t for particle %s\n",KeyW[Kw],Spec[n_t]);
			fprintf(log,"WARNING!Error while reading CS %s\t for particle %s\n",KeyW[Kw],Spec[n_t]);
			Serr += 1;
			break;
		}

		//линеаризация сечения
		for(k=0; k<K; k++)
		{
			A[k] = (Y[k+1]-Y[k])/(X[k+1]-X[k]);
			B[k] = Y[k]-A[k]*X[k];
		}
		//экстраполяция до лимита по энергии - сохранение наклона!!!!!!!!
		//A[K-1] = (A[K-2]+A[K-3]+A[K-4])/3;
		//B[K-1] = B[K-2];//(B[K-2]+B[K-3]+B[K-4])/3;

		//построение сечений по расчётной энергетической сетке
		k = 0;
		double Ei, CSj;
		for(i=0; i<NEmax; i++)
		{
			Ei = (i+0.5)*dEev;
			if(Ei<Ith[j])
				CSj = 0.0;
			else
			{
				if(Ei>X[k+1])//while(X[k+1]<Ei)
					k += 1;
				if(Ei>=X[K-1])
					k = K-1;

				CSj = A[k]*Ei+B[k];
				if(CSj<0)
				{
					//printf("WARNING!!!Cross section #%d is negative at point E=%.2lfeV\n",j,Ei);
					//fprintf(log,"WARNING!!!Cross section #%d is negative at point E=%.2lfeV\n",j,Ei);
					//Serr += 1;
					CSj = 0.0;
				}
			}

			CS[i][j] = CSj;//сечения,прочитанных из файла процессов, по точкам в центре ячеек сетки

			int I;
			for(d=1;d<=add;d++)
			{
				if(Ei<Ith[j+d])
					CS[i][j+d] = 0.0;
				else
				{
					I = i - int((Ith[j]-Ith[j+d])/dEev);//сдвиг порога по энергии
					if(I<0)
					{
						fprintf(log,"\nWARNING!!!Threshold shift is larger than Ith(target) - %s Cross section #%d \n",KeyW[Kw],j+d);
						break;
					}
					CS[I][j+d] = CS[i][j];

					if(i==NEmax-1)
					{
						while(I<=NEmax-1)
						{
							CS[I][j+d] = CS[i][j];//
							I += 1;
						}
					}
				}
			}

		}

		///Loging*****************************************************
		/*fprintf(log,"CS#%d at calc-energy net\n",j);
		for(i=0; i<NEmax; i++)
		{
			fprintf(log,"%.2lf\t%.2e\n",dEev*(i+0.5),CS[i][j]);
			i += 50;
		}
		fprintf(log,"\n");
		fclose(log);
		*/
		///Loging*****************************************************

		//пока не понял зачем
		//CS[j][0] = 0.0;
		//if(X[1]<dE1) {CXC[as][k][0]=Y[1];}

		j += 1+add;

		fscanf(cross,"%s",CSstr);
	}
	fclose(cross);

	//Printing diagnosic message*****************************************************
	printf("EEDF Cross Sections info:\n\n");
	if(Serr>0)
	{
		printf("ATTENTION!!!Number of errors = %d\n", Serr);
		fprintf(log,"ATTENTION!!!Number of errors = %d\n", Serr);
	}
	else
	{
		printf("Number of EEDF-processes = %d\n", j);
		fprintf(log,"Number of EEDF-processes = %d\n", j);
	}

	J = j;

	//Запись в Лог***********************************************************

	fprintf(log,"Log for CS-matrix\n\n");

	fprintf(log,"E,eV\t\t");
	for(j=0;j<J;j++)
		fprintf(log,"CSnum=%d\t\t",j+1);
	fprintf(log,"\n\n");

	for(k=0;k<NEmax;k++)
	{
		fprintf(log,"%.2lf\t",dEev*(k+0.5));
		for(j=0;j<J;j++)
			fprintf(log,"%.2e\t",CS[k][j]);
		fprintf(log,"\n");
		//k += 10;
	}
	fprintf(log,"\n\n");

	fprintf(log,"Log for CSR[][]-matrix\n\n");

	fprintf(log,"CSnum\t");
	for(i=0;i<5;i++)
		fprintf(log,"Ind=%d\t",i);
	fprintf(log,"\n\n");

	for(j=0;j<J;j++)
	{
		fprintf(log,"%d\t",j+1);
		for(i=0;i<5;i++)
			fprintf(log,"%d\t",CSR[j][i]);
		fprintf(log,"\n");
	}
	fprintf(log,"\n\n");

	int Nt = 3;
	fprintf(log,"Log for CSref[%s][][]-matrix\n\n",Spec[Nt]);

	fprintf(log,"CStype\t\t RSum for %s\t\t CS-num\n\n",Spec[Nt]);

	int jj;
	for(Kw=0;Kw<CStype;Kw++)
	{
		jj = CSref[Nt][Kw][0];//CSref[n_t][Kw][0];
		fprintf(log,"%s\t\t%d\t",KeyW[Kw],jj);
		if(jj>0)
		{
			for(j=1;j<=jj;j++)
				fprintf(log,"%d\t",CSref[Nt][Kw][j]);
		}
		fprintf(log,"\n");
	}
	fprintf(log,"\n\n");

	//**********добавление EEDF-реакций в файл кин.схемы**********************//

	int n0,n1,n2,n3;
	FILE *react;
	react = fopen("ReactionSet.txt", "w");

	fprintf(react,"%s\n",symb);
	fprintf(react,"\nEEDF_Reactions\n");
	fprintf(react,"\n%s\n",symb);
	jj = 0;
	for(j=0;j<J;j++)
	{
		n0 = CSR[j][0];//CS-type
		n1 = CSR[j][2]-1;//target
		n2 = CSR[j][3]-1;//product-1
		n3 = CSR[j][4]-1;//product-2

		if(n0==6)//ion
		{
			fprintf(react,"R%d\t e\t + %s\t ->\t  %s\t + e\t + e\t ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2]);
			jj += 1;
		}
		else if(n0>=7)//att+rec
		{
			fprintf(react,"R%d\t e\t + %s\t ->\t %s\t ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2]);
			jj += 1;
		}
		else
		{
			if(n3>0)
			{
				fprintf(react,"R%d\t e\t + %s\t ->\t e\t + %s\t + %s ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2],Spec[n3]);
				jj += 1;
			}
			else if(n2>0)
			{
				fprintf(react,"R%d\t e\t + %s\t ->\t e\t + %s ;\tEEDF\t//see_CS-set\n",jj+1,Spec[n1],Spec[n2]);
				jj += 1;
			}
		}
	}
	//fprintf(react,"\n");
	//fprintf(react,"// %s\n",symb);

	fclose(react);

	printf("%d of %d EEDF-processes were added to ReactionSet.txt\n\n",jj,J);
	fprintf(log,"%d of %d EEDF-processes were added to ReactionSet.txt\n\n",jj,J);

	fclose(log);

	return jj;
}
void EEDF_calc(double *Ne,double *Nni,int N,double *Te,double *dTe,double E,double Tgas,double Nel,double dte,double tic,int dot)//решение уравнения Больцмана
{
	int k,n,m,s,j,J,Jmax,nte;
    double Te0,Te1,Norm,E_kT;

	double A,B,C,F,den;
	double Vm[NEmax],Vmi[NEmax],D[NEmax],Vm_av;
	double Ur,Ul,Dr,Dl,Vr,Vl;
	double Ni[N];

	//Сетка по энергии:
	/*
	 Ne(k)[0]	[1]   [2]	[3]		   		 [k-1]  [k]  [k+1]						  [NEmax-1]
		|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|--x--|---------->E
	  k=0     1     2     3     4		   k-1    k    k+1                    			   NEmax
	*/

	for(n=0;n<N;n++)
        Ni[n] = Nni[n*(LEN+2)];

	nte = 0;
	Te1 = *Te;
	do
    //(nte=0;nte<20;nte++)
    {
        Te0 = Te1;

        //Difinition of D=A and Vm (Raizer)
        for(k=0;k<=NEmax-1;k++)
        {
            Vm[k] = 0.0;
            Vmi[k] = 0.0;
            double V = 0.0;
            for(n=1;n<N;n++)//по всем компонентам смеси (с кот. сталк. эл-ны), кроме электронов
            {
                //elastic collisions
                Jmax = CSref[n][0][0];//m=0

                if(Jmax>0)
                {
                    for(j=1;j<=Jmax;j++)
                    {
                        J = CSref[n][0][j]-1;

                        //CS[k][J] = 1e-15;//Test_for_Druvestein
                        //V = 2.0e11;//Test_for_Maxwell

                        V = CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ni[n];//sqrt(2/me)=4.69e13

                        Vm[k] += V;
                        Vmi[k] += V/Mi[n];
                    }
                }
            }

            D[k] = 2*e*e*E*E/(3*me*Vm[k]);// 2*e*e/(3*me*E0*E0)=1875.46*E^2/Vm - [B/см]^2*с
            Vmi[k] = 2*me*Vmi[k];//Sum(2*me*Vi/Mi);
        }

        //Цикл прогонки******************************************************
        //определение прогоночных коэффициентов
        double al[NEmax],bet[NEmax];
        double Qinel,Qex;

        for(k=0;k<=NEmax-1;k++)
        {
            //набор энергии от поля и упругие потери в соударениях
            if(k==0)
                Dl = D[k];
            else
                Dl = 0.5*(D[k]+D[k-1]);

            if(k==NEmax-1)
                Dr = D[k];
            else
                Dr = 0.5*(D[k+1]+D[k]);

            Vr = 0.5*(Vmi[k+1]+Vmi[k]);
            Vl = 0.5*(Vmi[k]+Vmi[k-1]);

            Ur = 0.5*Dr-Vr*(k+1)*dE;//point Ne[k]-right face
            Ul = 0.5*Dl-Vl*k*dE;//point Ne[k-1]-left face

            //without drift
            A = - Dl*k*dE/(dE*dE);//k-1(left)
            B = - Dr*(k+1)*dE/(dE*dE);//k+1(right)
            C = -(A + B) + 1/dte;//k(center)

            //with drift part
            if((Ur>=0)&&(Ul>=0))
            {
                A += -Ul/dE;
                C += Ur/dE;
            }
            if((Ur>=0)&&(Ul<0))
                C += (Ur-Ul)/dE;
            if((Ur<0)&&(Ul<0))
            {
                C += -Ul/dE;
                B += Ur/dE;
            }
            if((Ur<0)&&(Ul>=0))
            {
                A += -Ul/dE;
                B += Ur/dE;
            }

            //inelastic collisions*************************************************************************************
            Qinel = 0.0;
            s = 0;
            for(n=1;n<N;n++)//суммирование по всем компонентам смеси (с кот. сталк. эл-ны)
            {
                for(m=2;m<CStype;m++)//суммирование по всем типам процессов
                {
                    Jmax = CSref[n][m][0];//число пар (Ntarget,CStype)

                    if(Jmax>0)//только по ненулевым типам процессов
                    {
                        //excitation(vib,elec,diss)
                        if((m>=2)&&(m<=4))
                        {
                            for(j=1;j<=Jmax;j++)
                            {
                                J = CSref[n][m][j]-1;//номер сечения для пар (Ntarget,CStype)

                                s = k+int(Ith[J]/dEev);
                                if(s>=NEmax-1)
                                    s = NEmax-1;

                                Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k]+CS[s][J]*sqrt(2/me*dE*(s+0.5))*Ne[s];

                                Qinel += Qex*Ni[n];
                            }
                        }

                        //deexcitation(vib,elec)
                        if(m==5)
                        {
                            for(j=1;j<=Jmax;j++)
                            {
                                J = CSref[n][m][j]-1;

                                s = k+int(Ith[J]/dEev);
                                if(s>=NEmax-1)
                                    s = NEmax-1;

                                Qex = CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k]-CS[s][J]*sqrt(2/me*dE*(s+0.5))*Ne[s];
                                Qinel += Qex*Ni[n];
                            }
                        }

                        //attachment(recomb?)
                        if(m==6)
                        {
                            for(j=1;j<=Jmax;j++)
                            {
                                J = CSref[n][m][j]-1;

                                Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k];
                                Qinel += Qex*Ni[n];
                            }
                        }

                        //ionization_2-models
                        if(m==7)
                        {
                            for(j=1;j<=Jmax;j++)
                            {
                                J = CSref[n][m][j]-1;

                                /*//равнораспределение между 2-мя эл-нами********************************************
                                s = 2*k+int(Ith[J]/dEev);
                                if(s>=NEmax-1)
                                    s = NEmax-1;

                                Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k]+2*CS[s][J]*sqrt(2/me*dE*(s+0.5))*Ne[s];
                                */
                                //*********************************************************************************

                                //ионизирующий эл-н уносит остаток, новый уносит E0-Райзер+Хагеллар****************
                                E_kT = kb*Tgas;
                                s = k+int((Ith[J]+E_kT)/dEev);//E_kT=kT(Волошин)
                                if(s>=NEmax-1)
                                    s = NEmax-1;

                                Qex = 0.0;
                                if((k>=int(E_kT/dEev-0.5))&&((k<=int(E_kT/dEev+0.5))))
                                {
                                    int kI = int((Ith[J]+E_kT)/dEev);
                                    int ki = 0;
                                    for(ki=kI;ki<NEmax;ki++)
                                        Qex += CS[ki][J]*sqrt(2/me*dE*(ki+0.5))*Ne[ki];
                                }

                                Qex = -CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k]+CS[s][J]*sqrt(2/me*dE*(s+0.5))*Ne[s]+Qex;
                                //**********************************************************************************

                                Qinel += Qex*Ni[n];
                            }
                        }

                        //rotation??
                        //deexcitation??
                    }
                }
            }
            //*********************************************************************************************************

            //RH-part
            F = Ne[k]/dte + Qinel;

            if(k==0)
            {
                al[k+1] = -B/C;
                bet[k+1] = F/C;
            }
            else if(k==NEmax-1)
            {}
            else
            {
                den = 1.0/(A*al[k]+C);
                al[k+1] = -B*den;
                bet[k+1] = (F-A*bet[k])*den;
            }
        }

        //граничное условие**************************************************
        Ne[NEmax-1] = 0.0;
        //Ne[NEmax-1] = (F-A*bet[NEmax-1])/(A*al[NEmax-1]+C);

        //цикл обратной прогонки*********************************************
        for(k=NEmax-2;k>=0;k--)
        {
            Ne[k] = al[k+1]*Ne[k+1]+bet[k+1];
            if(Ne[k]<1e-30)
                Ne[k] = 0.0;
        }

        //Normalization***********************************************
        Norm = 0.0;
        for(k=0; k<NEmax; k++)
            Norm += Ne[k]*dEev;

        for(k=0; k<NEmax; k++)
            Ne[k] = Ne[k]*Nel/Norm;

        //************************************************************

        //Ee_Te-calculation*******************************************
        Ee = 0.0;
        for(k=0;k<NEmax;k++)
            Ee += Ne[k]*(k+0.5)*dEev*dEev;
        Ee = Ee/Nel;//[eV]
        Te1 = 2.0/3.0*Ee;//[eV]
        //************************************************************

        nte++;

    }while(fabs(Te1-Te0)>0.005);

    //Apply_to_used_matrices******************************************
    *Te = Te1;
    *dTe = Te1-Te0;

    //Vdr-calculation*************************************************

    Vm_av = 0.0;
    for(k=0;k<NEmax;k++)
        Vm_av += Vm[k]*Ne[k]*dEev;
    Vm_av *= 1.0/Nel;//
    Muel = e/me/Vm_av;
    //Vdr = Muel*E;
    //Jel = e*Nel*Vdr;//[СГС/cm2*s]
    //Qel = me/Mi[11]*Ee*Vm_av*Nel*1.602e-12;//[erg/cm3*s]

    //QE = E*Jel*3.14*pi*Rad*Rad; //[abs*СГС/s]

    //***************************************************************


    //Writing_data***************************************************
	//if(dot==Ndots)
		//EEDF_print(Ne,*Te,Nel,Norm,tic);

	//***************************************************************

}
void EEDF_const_calc(double *Ne,int N,double *Kel,int Nedf,double Nel,double tic)//Kel(EEDF)-calculation
{
	int I,j,k,J,Jmax,n,m;
    double Ki,K[Nedf];

    for(m=0;m<Nedf;m++)
    {
        Ki = 0.0;
        for(k=0;k<NEmax;k++)//int(Ith[m])
            Ki += CS[k][m+1]*sqrt(2/me*dE*(k+0.5))*Ne[k]*dEev;
        Kel[m] = Ki/Nel;
    }

    //Kel = K;

    /*
    I = 0;
    for(n=1;n<N;n++)//суммирование по всем компонентам смеси (с кот. сталк. эл-ны)
    {
        for(m=2;m<CStype;m++)//суммирование по всем типам процессов, кроме упругих столкновений и возб вращений
        {
            Jmax = CSref[n][m][0];//число пар (Ntarget,CStype)

            if(Jmax>0)//только по ненулевым типам процессов
            {
                for(j=1;j<=Jmax;j++)
                {
                    J = CSref[n][m][j]-1;//????
                    Ki = 0.0;
                    for(k=int(Ith[J]);k<NEmax;k++)
                        Ki += CS[k][J]*sqrt(2/me*dE*(k+0.5))*Ne[k]*dEev;
                    K[I] = Ki/Nel;
                    //Kel = K;

                    I++;
                }
            }
        }
    }
    */

    /*
    //Logging*******************************************
    FILE *log;
    if(tic==0)
        log = fopen("Kel_data.txt", "w");
    else
        log = fopen("Kel_data.txt", "a+");

    fprintf(log,"t=%.2e[s]\n",tic);
	for(m=0; m<Nedf; m++)
	    fprintf(log,"%d\t%.2e\n",m+1,Kel[m]);
	fprintf(log,"\n");
	fclose(log);
	//Logging*******************************************
	*/
}
void EEDF_print(double *Ne,double Te,double Nel,double Norm,double tic)//запись EEDF в файл
{
	int k;

    /*
	printf("\nt = %.2e[s]\n",tic);
	printf("Te = %.2lf[eV]\t Ne = %.2e[cm-3]\n",Te,Nel);
	*/

	FILE *log;
	log = fopen("EEDF_data.txt", "w");

    fprintf(log,"t=%.2e[s]\t Te=%.2lf[eV]\n",tic,Te);
	for(k=0; k<NEmax; k+=10)
        fprintf(log,"%.2lf\t%.2e\n",(k+0.5)*dEev,Ne[k]/Nel);
	fprintf(log,"\n\n");

	/*
	fprintf(log,"t=%.2e [s]\nNorm=%.2e [cm-3]\nTe=%.2lf [eV]\nE,[eV]\t Ne(E),[cm-3/eV]\n",dt*(nt+1),Norm,Te);//Nel
	for(k=0; k<NEmax; k+10)
		fprintf(log,"%.2lf\t%.2e\n",dEev*(k+0.5),Ne[k]/Nel);
	fprintf(log,"\n\n");
	*/
	fclose(log);
}
