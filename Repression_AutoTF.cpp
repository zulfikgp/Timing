// autorgeulation + target gene (repressed) + Decoy

//
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "time.h"

/* Memory allocation*/

double* dvmalloc(int n);
double* dvmalloc(int n) {
  double *vector;
  vector = (double *)malloc(sizeof(double)*n);
  for (int i=0;i<n;i++) {
    vector[i] = 0;
  }
  return vector;
}

int* ivmalloc(int n) {
  int *vector;
  vector = (int *)malloc(sizeof(int)*n);
  for (int i=0;i<n;i++) {
    vector[i] = 0;
  }
  return vector;
}

/****************************/
/* Find cumulative sum */

double* cumsum(double *A, int num);
double* cumsum(double *A, int num)
{
	double *CSum;
	CSum = dvmalloc(num);
	CSum[0] = A[0];
	for (int i=1;i<num;i++){
		CSum[i] = A[i] + CSum[i-1];
//		printf("%f %d\n",CSum[i],i);
	}
	return CSum;
}
/***********************/
/* Find index */

int FindInd(double *A, double flag);
int FindInd(double *A, double flag)
{
	int check = 0;
	int n=0, Ind;
	while(check==0){
		if (A[n]>=flag){
			check = 1;
			Ind = n;
//		printf("Psum %f Index %d\n",A[n],Ind);
		}
		n++;
		
	}
	return Ind;
}
/***********************/

/*Check if the site is available for initiation*/

int Initiate(int *A, int ind1, int LM);
int Initiate(int *A, int ind1, int LM)
{
	int check = 0;
	for (int i=ind1;i<ind1+LM;i++){
		if (A[i] == 1){
			check = 1;
		}
	}
	return check;
} 

/* Find a random floating number between min and max */

double unifrnd(double min, double max);
double unifrnd(double min, double max)
{
	double x1;
	/* x1 will be an element of [min,max] */
	x1=((double)rand()/RAND_MAX)*(max-min) + min;
	return x1;
}


/* FInd a random number from exponential distribution with mean K*/
/* x = -ln(1-y)/L, where y is uniform random number between 0 and 1*/
/* K = L */

double exprnd(double K);
double exprnd(double K)
{
	double x,y;
	y = unifrnd(0,0.999);
	x = -log(1-y)/K;
	return x;
}

/*** Gillespie for steady state**/
int* Gillespiess(int Decoy, double tend, int Nm1, int Nm2, int Np1, int Np2, int comp1, int comp2, int CDecoy, double ku1, double ku2, double kud, double gm0, double gp1, double gm2, double gp2, double fac);
int* Gillespiess(int Decoy, double tend, int Nm1, int Nm2, int Np1, int Np2, int comp1, int comp2, int CDecoy, double ku1, double ku2, double kud, double gm0, double gp1, double gm2, double gp2, double fac)
//int* Gillespie(int gene1, int gene2, int Decoy, double tend)
{	
	int gene1 = 1, gene2 = 1;
// Rates
	//double ku1 = 0.0069; // Unbinding first gene
	//double ku2 = 0.0069; // Unbinding first gene
	//double kud = 0.0024; // Binding/unbinding to decoy site

	double rp1 = 0.0003; //Degradation of protein
	double rp2 = 0.0003;

	double kb1 = 0.0027, kb2 = 0.0027, kbd = 0.0027; // Binding first gene
	//double gm2 = 0.11; // mRNA production defauls 0.33
	double gm1 = gm0*fac; // mRNA production	
	//double gp2 = 0.01; // Protein production
	double rm1 = 0.011, rm2 = 0.011; // Degradation of RNA

	double t = 0;

// Initiate the rates, To start all the rates are zero except the initiation of the first available site
	int NR = 17;
	double Rates[NR]; // Rate matrix
	Rates[0] = kb1*(gene1-comp1)*Np1;  // TF binding to gene1
	Rates[1] = kb2*(gene2-comp2)*Np1;  // TF binding to gene2
	Rates[2] = ku1*comp1; // Unbinding
	Rates[3] = ku2*comp2; // Unbinding
	Rates[4] = kbd*(Decoy-CDecoy)*Np1; // Binding of TF to Decoy
	Rates[5] = kud*CDecoy; // Unbinding

	Rates[6] = gm1*comp1 + gm0*(gene1-comp1); // mRNA production gene1
	Rates[7] = gm2*(gene2-comp2); // mRNA production gene2

	Rates[8] = gp1*Nm1; // protein production
	Rates[9] = gp2*Nm2; // protein production

	Rates[10] = rm1*Nm1; // mRNA degradation
	Rates[11] = rm2*Nm2; // mRNA degradation

	Rates[12] = rp1*Np1; // Protein degration
	Rates[13] = rp2*Np2; // Protein degration	

	Rates[14] = rp1*comp1; // Degradation of TFs bound to gene1
	Rates[15] = rp1*comp2; // gene2
	Rates[16] = rp1*CDecoy; // Degradation of TFs bound to decoy
	
	double *PSum, K, dt, flag;
	int Ind, Reaction;
	int *NOut;
	NOut = ivmalloc(7);

	while (t < tend){
		NOut[0] = Nm1;
		NOut[1] = Nm2;
		NOut[2] = Np1;
		NOut[3] = Np2;
		NOut[4] = comp1;
		NOut[5] = comp2;
		NOut[6] = CDecoy;

		PSum = cumsum(Rates,NR);
		K = PSum[NR-1];
		dt = exprnd(K);
		t = t+dt; // update time

/********* chose reaction **************/
		flag = unifrnd(0.00000001,K);
		Ind = FindInd(PSum,flag); // index of the reaction to execute
		free(PSum);
		Reaction = Ind; 
		switch (Reaction){
			case 0:
				comp1 = comp1 + 1; //Binding of protein1 to gene1
				Np1 = Np1 - 1;
				break;
				
			case 1:
				comp2 = comp2 + 1; //Binding of protein1 to gene1
				Np1 = Np1 - 1;
				break;
				
			case 2:
				comp1 = comp1 - 1; // unbinding
				Np1 = Np1 + 1;
				break;

			case 3:
				comp2 = comp2 - 1; // unbinding
				Np1 = Np1 + 1;
				break;
				
			case 4:
				CDecoy = CDecoy + 1; // Binding of protein to Decoy
				Np1 = Np1 - 1;
				break;
				
			case 5:
				CDecoy = CDecoy - 1; // Unbinding of protein to Decoy
				Np1 = Np1 + 1;
				break;
				
			case 6:
				Nm1 = Nm1 + 1; // mrna production
				break;
				
			case 7:
				Nm2 = Nm2 + 1; // mrna production
				break;
				
			case 8:
				Np1 = Np1 + 1;
				break;

			case 9:
				Np2 = Np2 + 1;
				break;
				
			case 10:
				Nm1 = Nm1 - 1; // mrna degradation
				break;
				
			case 11:
				Nm2 = Nm2 - 1; // mrna degradation
				break;
				
			case 12:
				Np1 = Np1 - 1;
				break;
				
			case 13:
				Np2 = Np2 - 1;
				break;
				
			case 14:
				comp1 = comp1 - 1;
				break;
				
			case 15:
				comp2 = comp2 - 1;
				break;
				
			default:
				CDecoy = CDecoy - 1;
		}
/*************** Recalculate rates ************************************/
		Rates[0] = kb1*(gene1-comp1)*Np1;  // TF binding to gene1
		Rates[1] = kb2*(gene2-comp2)*Np1;  // TF binding to gene2
		Rates[2] = ku1*comp1; // Unbinding
		Rates[3] = ku2*comp2; // Unbinding

		Rates[4] = kbd*(Decoy-CDecoy)*Np1; // Binding of TF to Decoy
		Rates[5] = kud*CDecoy; // Unbinding

		Rates[6] = gm1*comp1 + gm0*(gene1-comp1); // mRNA production gene1
		Rates[7] = gm2*(gene2-comp2); // mRNA production gene2

		Rates[8] = gp1*Nm1; // protein production
		Rates[9] = gp2*Nm2; // protein production

		Rates[10] = rm1*Nm1; // mRNA degradation
		Rates[11] = rm2*Nm2; // mRNA degradation

		Rates[12] = rp1*Np1; // Protein degration
		Rates[13] = rp2*Np2; // Protein degration	

		Rates[14] = rp1*comp1; // Degradation of TFs bound to gene1
		Rates[15] = rp1*comp2; // gene2
		Rates[16] = rp1*CDecoy; // Degradation of TFs bound to decoy
	}
	return NOut;
}

/****************************/

double* Gillespie(int Decoy, int Nm1, int Nm2, int Np1, int Np2, int comp1, int comp2, int CDecoy, double Tst, double Tsf, double Tss, double tend, double ku1, double ku2, double kud, double gm0, double gp1, double gm2, double gp2, double fac);
double* Gillespie(int Decoy, int Nm1, int Nm2, int Np1, int Np2, int comp1, int comp2, int CDecoy, double Tst, double Tsf, double Tss, double tend, double ku1, double ku2, double kud, double gm0, double gp1, double gm2, double gp2, double fac)
//int* Gillespie(int gene1, int gene2, int Decoy, double tend)
{	
// Gene1 TF complex, mRNA1, mRNA2, Protein1, Protein2, Decoy complex
//	int comp1 = 0, comp2 = 0, Nm1 = 0, Nm2 = 0, Np1 = 0, Np2 = 0, CDecoy = 0; 
	int gene1 = 1, gene2 = 1;
// Rates
	//double ku1 = 0.0909; // Unbinding first gene
	//double ku2 = 0.0909; // Unbinding first gene
	//double kud = 0.0024; // Binding/unbinding to decoy site

	double rp1 = 0.0003; //Degradation of protein
	double rp2 = 0.0003;

	double kb1 = 0.0027, kb2 = 0.0027, kbd = 0.0027; // Binding first gene
	//double gm2 = 0.11; // mRNA production defauls 0.33
	double gm1 = gm0*fac; // mRNA production	
	//double gp2 = 0.01; // Protein production
	double rm1 = 0.011, rm2 = 0.011; // Degradation of RNA

	double t = 0, t1 = 0, t2 = 0, t3 = 0;


// Initiate the rates, To start all the rates are zero except the initiation of the first available site
	int NR = 17;
	double Rates[NR]; // Rate matrix
	Rates[0] = kb1*(gene1-comp1)*Np1;  // TF binding to gene1
	Rates[1] = kb2*(gene2-comp2)*Np1;  // TF binding to gene2
	Rates[2] = ku1*comp1; // Unbinding
	Rates[3] = ku2*comp2; // Unbinding

	Rates[4] = kbd*(Decoy-CDecoy)*Np1; // Binding of TF to Decoy
	Rates[5] = kud*CDecoy; // Unbinding

	Rates[6] = gm1*comp1 + gm0*(gene1-comp1); // mRNA production gene1
	Rates[7] = gm2*(gene2-comp2); // mRNA production gene2

	Rates[8] = gp1*Nm1; // protein production
	Rates[9] = gp2*Nm2; // protein production

	Rates[10] = rm1*Nm1; // mRNA degradation
	Rates[11] = rm2*Nm2; // mRNA degradation

	Rates[12] = rp1*Np1; // Protein degration
	Rates[13] = rp2*Np2; // Protein degration	

	Rates[14] = rp1*comp1; // Degradation of TFs bound to gene1
	Rates[15] = rp1*comp2; // gene2
	Rates[16] = rp1*CDecoy; // Degradation of TFs bound to decoy
	
		
//Start Gillespie
	FILE *filet;	
	//filet = fopen("Time-series.txt","w");
	double *PSum, K, dt, flag;
	int Ind, Reaction;
	int check1 = 0, check2 = 0, check3 = 0, check = 0, Ttot = 0;

	//while (Np2< 0.5*Tss){
	while ( check == 0){
		PSum = cumsum(Rates,NR);
		K = PSum[NR-1];
		dt = exprnd(K);
//		printf("%f %f\n",K,dt);

		t = t+dt; // update time

/********* chose reaction **************/
		flag = unifrnd(0.00000001,K);
		Ind = FindInd(PSum,flag); // index of the reaction to execute
		free(PSum);
		Reaction = Ind; 
		switch (Reaction){
			case 0:
				comp1 = comp1 + 1; //Binding of protein1 to gene1
				Np1 = Np1 - 1;
				break;
				
			case 1:
				comp2 = comp2 + 1; //Binding of protein1 to gene1
				Np1 = Np1 - 1;
				break;
				
			case 2:
				comp1 = comp1 - 1; // unbinding
				Np1 = Np1 + 1;
				break;

			case 3:
				comp2 = comp2 - 1; // unbinding
				Np1 = Np1 + 1;
				break;
				
			case 4:
				CDecoy = CDecoy + 1; // Binding of protein to Decoy
				Np1 = Np1 - 1;
				break;
				
			case 5:
				CDecoy = CDecoy - 1; // Unbinding of protein to Decoy
				Np1 = Np1 + 1;
				break;
				
			case 6:
				Nm1 = Nm1 + 1; // mrna production
				break;
				
			case 7:
				Nm2 = Nm2 + 1; // mrna production
				break;
				
			case 8:
				Np1 = Np1 + 1;
				break;

			case 9:
				Np2 = Np2 + 1;
				break;
				
			case 10:
				Nm1 = Nm1 - 1; // mrna degradation
				break;
				
			case 11:
				Nm2 = Nm2 - 1; // mrna degradation
				break;
				
			case 12:
				Np1 = Np1 - 1;
				break;
				
			case 13:
				Np2 = Np2 - 1;
				break;
				
			case 14:
				comp1 = comp1 - 1;
				break;
				
			case 15:
				comp2 = comp2 - 1;
				break;
				
			default:
				CDecoy = CDecoy - 1;
		}
		Ttot = Np1 + comp1 + comp2 + CDecoy;
		if (Np2 >= 0.5*Tss & check1 == 0){
			check1 = 1;
			t3 = t;
			}
	
		if (Np1 >= 0.5*Tsf & check2 == 0){
			check2 = 1;
			t2 = t;
			}
		if (Ttot >= 0.5*Tst & check3 == 0){
			check3 = 1;
			t1 = t;
			}
				
		check = check1*check2*check3;
	
/*************** Recalculate rates ************************************/
		Rates[0] = kb1*(gene1-comp1)*Np1;  // TF binding to gene1
		Rates[1] = kb2*(gene2-comp2)*Np1;  // TF binding to gene2
		Rates[2] = ku1*comp1; // Unbinding
		Rates[3] = ku2*comp2; // Unbinding

		Rates[4] = kbd*(Decoy-CDecoy)*Np1; // Binding of TF to Decoy
		Rates[5] = kud*CDecoy; // Unbinding

		Rates[6] = gm1*comp1 + gm0*(gene1-comp1); // mRNA production gene1
		Rates[7] = gm2*(gene2-comp2); // mRNA production gene2

		Rates[8] = gp1*Nm1; // protein production
		Rates[9] = gp2*Nm2; // protein production

		Rates[10] = rm1*Nm1; // mRNA degradation
		Rates[11] = rm2*Nm2; // mRNA degradation

		Rates[12] = rp1*Np1; // Protein degration
		Rates[13] = rp2*Np2; // Protein degration	

		Rates[14] = rp1*comp1; // Degradation of TFs bound to gene1
		Rates[15] = rp1*comp2; // gene2
		Rates[16] = rp1*CDecoy; // Degradation of TFs bound to decoy

		}
	double *TOut;
	TOut = dvmalloc(3);
	TOut[0] = t1; // Total TF
	TOut[1] = t2; // Free TF
	TOut[2] = t3; // Target
	return TOut;
}

/********* MAIN CODE ***********/
int main()
{
	srand(time(NULL));
	time_t tstart, tfinish;
	tstart = time (NULL);
	FILE *file, *filer;

	int Decoy = 0;
	double Tss = 0, Tst = 0, Occ1, Occ2, Occ3, Tsf = 0;

	double ku1, ku2, kud;
	int Nit = 1*pow(10,5), Nit2 = 5*pow(10,4);
	double tend = 2*0.5*pow(10,5);
	double *TRes;
	int TF = 50, bufT;

	double bufr;
	filer = fopen("Infile.txt","r");
	fscanf(filer,"%lf",&bufr);
	ku1 = bufr;
	fscanf(filer,"%lf",&bufr);
	ku2 = bufr;
	fscanf(filer,"%lf",&bufr);
	kud = bufr;
	fscanf(filer,"%d",&bufT);
	TF = bufT;	
	//printf("%f %f %f\n",ku1,ku2,kud);
	fclose(filer);
	
	//ku2 = ku1, kud = ku1;
	char str1[200], buf1[10],buf2[10],buf3[10], buf4[10];
	sprintf(buf1,"%1.4f",ku1);
	sprintf(buf2,"%1.4f",ku2);
	sprintf(buf3,"%1.4f",kud);
	sprintf(buf4,"%d",TF);
	//strcpy(str1,"CLUSTER/VaryCons_Decoy_rp-0.0003_avTF-100_ku1-");
	

	strcpy(str1,"DATA/REPRESSION/Auto_Varyku2_AS-10.00_gm0-0.025_rp-0.0003_avTF-");
	strcat(str1,buf4);
	strcat(str1,"_ku1-");
	strcat(str1,buf1);
	strcat(str1,"_ku2-");
	strcat(str1,buf2);
	strcat(str1,"_kud-");
	strcat(str1,buf3);
	strcat(str1,".txt");
	file = fopen(str1,"w");

	double avgt1, avgt2, avgt3, avgt12,avgt22, avgt33, Var1, Var2, Var3;
	int Nin[7], *NOut, *Ntemp;
	double gm0 = 0.025, gp1 = 0.1; // TF production rates
	double gm2 = 0.05, gp2 = 0.1, gm1;
	double fac = 10.0; // fac > 1 (activator), <1 repressor
	double sig1, sig2, theta;
	double fccc[22] = {0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	double KUV[16] = {0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.15, 0.2, 0.3, 0.4,  0.6, 0.8, 1};
	for (int i=0;i<16;i++){	
		Decoy = 0;//ND[i];
//		ku1 = 0.002;
		ku2 = KUV[i] ;
//		fac = 20+i*20;//fccc[i];
		avgt1=0, avgt12=0, avgt2 = 0, avgt3 = 0, avgt22 = 0, avgt33=0, Tss = 0, Tst = 0, Tsf = 0, Occ1 = 0, Occ2 = 0, Occ3 = 0;

/*************** Tune transcription to get fixed TF **************/
		gm1 = fac*gm0;
		sig1 = 0.0027/(0.0003 + ku1);
		sig2 = 0.0027/(0.0003 + ku2);
		theta = 1/fac; // ratio of basal and activated rate of transcription
		//gm1 = TF*0.0003*0.011*(1+TF*sig1)*(1 + sig1/(1+TF*sig1) + sig2/(1+TF*sig2)) / (gp1*(theta + TF*sig1));
		gp1 = TF*0.0003*0.011*(1+TF*sig1)*(1 + sig1/(1+TF*sig1) + sig2/(1+TF*sig2)) / ((gm0 + TF*sig1*gm1));
/*****************************************************************/


		/*Find steady state values for TF and target*/
		Ntemp = Gillespiess(Decoy, pow(10,6), 0, 0, 0, 0, 0, 0, 0,ku1,ku2,kud,gm0,gp1,gm2,gp2, fac);
		Nin[0] = Ntemp[0], Nin[1] = Ntemp[1], Nin[2] = Ntemp[2], Nin[3] = Ntemp[3], Nin[4] = Ntemp[4], Nin[5] = Ntemp[5], Nin[6] = Ntemp[6];
		free(Ntemp);
		for (int ni=0;ni<Nit2;ni++){
			NOut = Gillespiess(Decoy, tend, Nin[0], Nin[1], Nin[2], Nin[3], Nin[4], Nin[5], Nin[6],ku1,ku2,kud,gm0,gp1,gm2,gp2,fac);
			Nin[0] = NOut[0], Nin[1] = NOut[1], Nin[2] = NOut[2], Nin[3] = NOut[3], Nin[4] = NOut[4], Nin[5] = NOut[5], Nin[6] = NOut[6];
			Tst = Tst + NOut[2] + NOut[4] + NOut[5] + NOut[6]; // total TF 
			Tss = Tss + NOut[3]; // Target
			Tsf = Tsf + NOut[2]; // free TF
			Occ1 = Occ1 + NOut[4];
			Occ2 = Occ2 + NOut[5];
			Occ3 = Occ3 + NOut[6];
			free(NOut);
		}
		Tst = Tst/Nit2;
		Tsf = Tsf/Nit2;
		Tss = Tss/Nit2; // average target expression
		Occ1 = Occ1/Nit2;
		Occ2 = Occ2/Nit2;
		Occ3 = Occ3/Nit2;
		//printf("%f %f\n",Tst,Tss);
		/************** */
	
		/*Response time*/
		for (int ni=0;ni<Nit;ni++){
			TRes = Gillespie(Decoy, 0, 0, 0, 0, 0, 0, 0, Tst, Tsf, Tss, tend, ku1, ku2, kud, gm0, gp1, gm2, gp2, fac);
			avgt1 = avgt1 + TRes[0]; //Total TF
			avgt2 = avgt2 + TRes[1]; // Free TF
			avgt3 = avgt3 + TRes[2]; // Target
			avgt12 = avgt12 + TRes[0]*TRes[0];
			avgt22 = avgt22 + TRes[1]*TRes[1];
			avgt33 = avgt33 + TRes[2]*TRes[2];
		}
		Var1 = avgt12/Nit-pow(avgt1/Nit,2);
		Var2 = avgt22/Nit-pow(avgt2/Nit,2);
		Var3 = avgt33/Nit-pow(avgt3/Nit,2);
	
		fprintf(file,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",fac,gm0,gp1,ku1,ku2,avgt1/Nit,avgt2/Nit,avgt3/Nit,Var1,Var2,Var3,Tst,Tsf,Tss,Occ1,Occ2,Occ3);
		printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",fac,gm0,gp1,ku1,ku2,avgt1/Nit,avgt2/Nit,avgt3/Nit,Var1,Var2,Var3,Tst,Tsf,Tss,Occ1,Occ2,Occ3);
		tfinish = time (NULL);
		long comp_time = (tfinish - tstart);
		printf("Computation time = %ld\n", comp_time);
	}
}
