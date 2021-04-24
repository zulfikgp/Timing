// All protein decays
// Target starts with saturation value
// TFs are negatively autoregulated
// 1 target + N Decoy
// at t=0 everything is set to zero


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
	int n=0, Ind = 0;
	while(check==0){
		if (A[n]>=flag){
//			if (A[n]==flag){printf("Ambiguous");}
			check = 1;
			Ind = n;
		//printf("Psum %f Index %d\n",A[n],Ind);
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

/************ Steady state ****************/

int* Gillespiess(double tend, double ku1, double ku2, double kud, int Decoy, double rp1, double gm1, double gp1, int Nm1, int Nm2, int Np1, int Np2, int comp1, int comp2, int CDecoy, double gm2, double gp2);
int* Gillespiess(double tend, double ku1, double ku2, double kud, int Decoy, double rp1, double gm1, double gp1, int Nm1, int Nm2, int Np1, int Np2, int comp1, int comp2, int CDecoy, double gm2, double gp2)
//int* Gillespie(int gene1, int gene2, int Decoy, double tend)
{

// Gene1 TF complex, mRNA1, mRNA2, Protein1, Protein2, Decoy complex
	int gene1 = 1, gene2 = 1; 
	
// Rates
/* Unbinding rates Oid = 1/420 = 0.0024, O1 = 1/144 = 0.0069, O2 = 1/11 = 0.0909, O3 = 1 / 0.47 = 2.1277 */ 
	double kb1 = 0.0027; // Binding first gene
	double kb2 = 0.0027; // Binding second gene
	double kbd = 0.0027; // Binding/unbinding to decoy site

	//double gm1 = gm, gm2 = gm; // mRNA production
	//double gp1 = gp, gp2 = gp; // Protein production
	double rp2 = rp1 ; //Degradation of protein
	double rm1 = 0.011, rm2 = 0.011; // Degradation of RNA

// Initiate the rates
	int NR = 17;
	double Rates[NR]; // Rate matrix
	Rates[0] = kb1*(gene1-comp1)*Np1;
	Rates[1] = kb2*(gene2-comp2)*Np1;
	Rates[2] = ku1*comp1;
	Rates[3] = ku2*comp2;

	Rates[4] = kbd*(Decoy-CDecoy)*Np1;
	Rates[5] = kud*CDecoy;

	Rates[6] = gm1*(gene1-comp1); // mRNA production
	Rates[7] = gm2*(gene2-comp2);

	Rates[8] = gp1*Nm1; // protein production
	Rates[9] = gp2*Nm2;

	Rates[10] = rm1*Nm1; // mRNA degradation
	Rates[11] = rm2*Nm2;

	Rates[12] = rp1*Np1; // Protein degration
	Rates[13] = rp2*Np2;
	
	Rates[14] = rp1*comp1; //Degradation of TFs bound to gene1
	Rates[15] = rp1*comp2; //Degradation of TFs bound to gene2
	Rates[16] = rp1*CDecoy; //Degradation of TFs bound to decoy

	//for (int i=0;i<17;i++){printf("The rates are %f  %f %f\n",Rates[i],gm1,gm2);}
		
//Start Gillespie
	double *PSum, K, dt, flag, t=0;
	int Ind, Reaction, cnt1 = 0, cnt2 = 0, Ttot = 0;
	int *NOut;
	NOut = ivmalloc(7);
	 
	while (t < tend){

		NOut[0] = Nm1; NOut[1] = Nm2; NOut[2] = Np1; NOut[3] = Np2;
		NOut[4] = comp1; NOut[5] = comp2; NOut[6] = CDecoy;

		PSum = cumsum(Rates,NR);
		K = PSum[NR-1];
		dt = exprnd(K);
		//printf("%f %f\n",K,dt);

		t = t+dt; // update time

/********* chose reaction **************/
		flag = unifrnd(0.00000001,K);
		//printf("Error due to %f %f\n",flag,K);
		Ind = FindInd(PSum,flag); // index of the reaction to execute
		free(PSum);
		Reaction = Ind; 
		switch (Reaction){
			case 0:
				comp1 = comp1 + 1; //Binding of protein1 to gene1
				Np1 = Np1 - 1;
				break;
				
			case 1:
				comp2 = comp2 + 1; //Binding of protein to gene2
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
				Nm2 = Nm2 + 1; // mRNA production target
				break;
				
			case 8:
				Np1 = Np1 + 1; // TF production
				break;
				
			case 9:
				Np2 = Np2 + 1; //Target production
				break;
				
			case 10:
				Nm1 = Nm1 - 1; // mrna degradation
				break;
				
			case 11:
				Nm2 = Nm2 - 1;
				break;
				
			case 12:
				Np1 = Np1 - 1; // Protein degration	
				//cnt1 = cnt1 + 1;
				break;
				
			case 13:
				Np2 = Np2 - 1;
				//cnt2 = cnt2 + 1;
				break;
				
			case 14:
				comp1 = comp1 - 1;
				//cnt1 = cnt1 + 1;
				break;
			
			case 15:
				comp2 = comp2 - 1;
				//cnt1 = cnt1 + 1;
				break;
				
			default:
				CDecoy = CDecoy - 1;
				//cnt1 = cnt1 + 1;

		}

/*************** Recalculate rates ************************************/
		Rates[0] = kb1*(gene1-comp1)*Np1;
		Rates[1] = kb2*(gene2-comp2)*Np1;
		Rates[2] = ku1*comp1;
		Rates[3] = ku2*comp2;

		Rates[4] = kbd*(Decoy-CDecoy)*Np1;
		Rates[5] = kud*CDecoy;

		Rates[6] = gm1*(gene1-comp1); // mRNA production
		Rates[7] = gm2*(gene2-comp2);
	
		Rates[8] = gp1*Nm1; // protein production
		Rates[9] = gp2*Nm2;

		Rates[10] = rm1*Nm1; // mRNA degradation
		Rates[11] = rm2*Nm2;

		Rates[12] = rp1*Np1; // Protein degration
		Rates[13] = rp2*Np2;
		
		Rates[14] = rp1*comp1; //Degradation of TFs bound to gene1
		Rates[15] = rp1*comp2; //Degradation of TFs bound to gene2
		Rates[16] = rp1*CDecoy; //Degradation of TFs bound to decoy

		//fprintf(filet,"%f %d %d\n",t,Np1+comp1+comp2+CDecoy,Np2);		
	}

	return NOut;
}


/************** Timing *************/

double* Gillespie(double ku1, double ku2, double kud, int Decoy, double Tst, double Tsf, double Tss, double rp1, double gm1, double gp1, double gm2, double gp2);
double* Gillespie(double ku1, double ku2, double kud, int Decoy, double Tst, double Tsf, double Tss, double rp1, double gm1, double gp1, double gm2, double gp2)
//int* Gillespie(int gene1, int gene2, int Decoy, double tend)
{

// Gene1 TF complex, mRNA1, mRNA2, Protein1, Protein2, Decoy complex
	int gene1 = 1, gene2 = 1, comp1 = 0, comp2 = 0, Nm1 = 0, Nm2 = 0, Np1 = 0, Np2 = 0, CDecoy = 0; 
	
// Rates
/* Unbinding rates Oid = 1/420 = 0.0024, O1 = 1/144 = 0.0069, O2 = 1/11 = 0.0909, O3 = 1 / 0.47 = 2.1277 */ 
	double kb1 = 0.0027; // Binding first gene
	double kb2 = 0.0027; // Binding second gene
	double kbd = 0.0027; // Binding/unbinding to decoy site

	//double gm1 = gm, gm2 = gm; // mRNA production
	//double gp1 = gp, gp2 = gp; // Protein production
	double rp2 = rp1 ; //Degradation of protein
	double rm1 = 0.011, rm2 = 0.011; // Degradation of RNA

	double t = 0, t1 = 0, t2 = 0, t3 = 0;
	int check1 = 0, check2 = 0, check = 0,check3 = 0;


// Initiate the rates
	int NR = 17;
	double Rates[NR]; // Rate matrix
	Rates[0] = kb1*(gene1-comp1)*Np1;
	Rates[1] = kb2*(gene2-comp2)*Np1;
	Rates[2] = ku1*comp1;
	Rates[3] = ku2*comp2;

	Rates[4] = kbd*(Decoy-CDecoy)*Np1;
	Rates[5] = kud*CDecoy;

	Rates[6] = gm1*(gene1-comp1); // mRNA production
	Rates[7] = gm2*(gene2-comp2);

	Rates[8] = gp1*Nm1; // protein production
	Rates[9] = gp2*Nm2;

	Rates[10] = rm1*Nm1; // mRNA degradation
	Rates[11] = rm2*Nm2;

	Rates[12] = rp1*Np1; // Protein degration
	Rates[13] = rp2*Np2;
	
	Rates[14] = rp1*comp1; //Degradation of TFs bound to gene1
	Rates[15] = rp1*comp2; //Degradation of TFs bound to gene2
	Rates[16] = rp1*CDecoy; //Degradation of TFs bound to decoy
	
//Start Gillespie
	double *PSum, K, dt, flag;
	int Ind, Reaction, cnt1 = 0, cnt2 = 0, Ttot = 0;;

	while (check == 0){

		PSum = cumsum(Rates,NR);
		K = PSum[NR-1];
		dt = exprnd(K);
//		printf("%f %f\n",K,dt);

		t = t+dt; // update time

/********* chose reaction **************/
		flag = unifrnd(0.00000001,K);
//		printf("%f %f\n",flag,K);
		Ind = FindInd(PSum,flag); // index of the reaction to execute
		free(PSum);
		Reaction = Ind; 
		switch (Reaction){
			case 0:
				comp1 = comp1 + 1; //Binding of protein1 to gene1
				Np1 = Np1 - 1;
				break;
				
			case 1:
				comp2 = comp2 + 1; //Binding of protein to gene2
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
				Nm2 = Nm2 + 1; // mRNA production target
				break;
				
			case 8:
				Np1 = Np1 + 1; // TF production
				break;
				
			case 9:
				Np2 = Np2 + 1; //Target production
				break;
				
			case 10:
				Nm1 = Nm1 - 1; // mrna degradation
				break;
				
			case 11:
				Nm2 = Nm2 - 1;
				break;
				
			case 12:
				Np1 = Np1 - 1; // Protein degration	
				//cnt1 = cnt1 + 1;
				break;
				
			case 13:
				Np2 = Np2 - 1;
				//cnt2 = cnt2 + 1;
				break;
				
			case 14:
				comp1 = comp1 - 1;
				//cnt1 = cnt1 + 1;
				break;
			
			case 15:
				comp2 = comp2 - 1;
				//cnt1 = cnt1 + 1;
				break;
				
			default:
				CDecoy = CDecoy - 1;
				//cnt1 = cnt1 + 1;

		}

/*************** Recalculate rates ************************************/
		Rates[0] = kb1*(gene1-comp1)*Np1;
		Rates[1] = kb2*(gene2-comp2)*Np1;
		Rates[2] = ku1*comp1;
		Rates[3] = ku2*comp2;

		Rates[4] = kbd*(Decoy-CDecoy)*Np1;
		Rates[5] = kud*CDecoy;

		Rates[6] = gm1*(gene1-comp1); // mRNA production
		Rates[7] = gm2*(gene2-comp2);
	
		Rates[8] = gp1*Nm1; // protein production
		Rates[9] = gp2*Nm2;

		Rates[10] = rm1*Nm1; // mRNA degradation
		Rates[11] = rm2*Nm2;

		Rates[12] = rp1*Np1; // Protein degration
		Rates[13] = rp2*Np2;
		
		Rates[14] = rp1*comp1; //Degradation of TFs bound to gene1
		Rates[15] = rp1*comp2; //Degradation of TFs bound to gene2
		Rates[16] = rp1*CDecoy; //Degradation of TFs bound to decoy

		//fprintf(filet,"%f %d %d\n",t,Np1+comp1+comp2+CDecoy,Np2);		

		Ttot = Np1+comp1+comp2+CDecoy;
		if (Ttot >= 0.5*Tst & check1 == 0){
			t1 = t;
			check1 = 1;
		}
		if (Np1 >= 0.5*Tsf  & check2 == 0){
			t2 = t;
			check2 = 1;
		}
		if (Np2 >= 0.5*Tss  & check3 == 0){
			t3 = t;
			check3 = 1;
		}
		check = check1*check2*check3;
	}
	//printf("%f %f\n", t1, t2);
	//fclose(filet);
	double *NOut;
	NOut = dvmalloc(3);
	 
	NOut[0] = t1;
	NOut[1] = t2;
	NOut[2] = t3;

	return NOut;
}


/**********************************************/
/**********************************************/

double* run_sim(double ku1, double ku2, double kud, double rpr, int Decoy, double gmt, double gpt, double gm2, double gp2);
double* run_sim(double ku1, double ku2, double kud, double rpr, int Decoy, double gmt, double gpt, double gm2, double gp2)

{
	srand(time(NULL));
	time_t tstart, tfinish;
	tstart = time (NULL);
	FILE *file;
	
	double DT, *DT2;
	int Nit2 = 4*pow(10,4), Nit = 1*pow(10,5);
	double avg1=0, avg2=0, avg3 = 0, avg12=0, avg22=0, avg33 = 0, Var1, Var2, Var3, Tst=0, Tss=0, Tsf = 0, Occ1 = 0, Occ2 = 0, Occ3 = 0;
	int Nin[7], *NOut, *Ntemp;
	double tend = pow(10,5);

	/*Find steady state values for TF, target*/


		Ntemp = Gillespiess(pow(10,6), ku1, ku2, kud, Decoy, rpr, gmt, gpt, 0, 0, 0, 0, 0, 0, 0,gm2,gp2);
		Nin[0] = Ntemp[0], Nin[1] = Ntemp[1], Nin[2] = Ntemp[2], Nin[3] = Ntemp[3];
		Nin[4] = Ntemp[4], Nin[5] = Ntemp[5], Nin[6] = Ntemp[6];
		free(Ntemp);
		for (int ni=0;ni<Nit2;ni++){
			NOut = Gillespiess(tend, ku1, ku2, kud, Decoy, rpr, gmt, gpt, Nin[0], Nin[1], Nin[2], Nin[3], Nin[4], Nin[5], Nin[6],gm2,gp2);
			Nin[0] = NOut[0], Nin[1] = NOut[1], Nin[2] = NOut[2], Nin[3] = NOut[3];
			Nin[4] = NOut[4], Nin[5] = NOut[5], Nin[6] = NOut[6];
			Tsf = Tsf + NOut[2]; // TF 
			Tss = Tss + NOut[3]; // Target
			Tst = Tst + NOut[2] + NOut[4] + NOut[5] + NOut[6];
			Occ1 = Occ1 + NOut[4]; // TF Occupanc
			Occ2 = Occ2 + NOut[5]; // Target occupancy
			Occ3 = Occ3 + NOut[6]; // Decoy occupancy
			free(NOut);
		}
		Tst = Tst/Nit2;
		Tsf = Tsf/Nit2;
		Tss = Tss/Nit2; // average target expression
		Occ1= Occ1/Nit2;
		Occ2= Occ2/Nit2;
		Occ3= Occ3/Nit2;
		printf("%f %f %f\n",Tst,Tsf,Tss);
	
	/*Find response time*/

	for (int ni=0;ni<Nit;ni++){
		DT2 = Gillespie(ku1, ku2, kud, Decoy, Tst, Tsf, Tss, rpr, gmt, gpt, gm2, gp2);
		//fprintf(file,"%f %f\n",DT2[0],DT2[1]);
		//printf("%f %f\n",DT2[0],DT2[1]);
		avg1 = avg1 + DT2[0];
		avg2 = avg2 + DT2[1];
		avg3 = avg3 + DT2[2];
		avg12 = avg12 + DT2[0]*DT2[0];
		avg22 = avg22 + DT2[1]*DT2[1];
		avg33 = avg33 + DT2[2]*DT2[2];
		free(DT2);
	}
	Var1 = avg12/Nit-pow(avg1/Nit,2);
	Var2 = avg22/Nit-pow(avg2/Nit,2);
	Var3 = avg33/Nit-pow(avg3/Nit,2);
	
	double *NT;
	NT = dvmalloc(12);
	NT[0] = avg1/Nit;
	NT[1] = avg2/Nit;
	NT[2] = avg3/Nit;
	NT[3] = Var1;
	NT[4] = Var2;
	NT[5] = Var3;
	NT[6] = Tst;
	NT[7] = Tsf;
	NT[8] = Tss;
	NT[9] = Occ1;
	NT[10] = Occ2;
	NT[11] = Occ3;
	//fclose(file);
	tfinish = time (NULL);
	long comp_time = (tfinish - tstart);
	printf("Computation time = %ld, \n", comp_time);
	
	return NT;
}

/********* MAIN CODE ***********/
int main()
{
	time_t tstart, tfinish;
	tstart = time (NULL);

	double ku1 = 0.01, ku2 = 0.0909, kud = 0.001;

	char str1[200], buf[10], buf2[10], buf3[10], buf4[10], buf5[10];


	int ND[15] = {0,10,20,30,40,50,60,70,80,90,100,150,200,250,300}, Decoy = 0;
	double kuv[20] = {0.0002,0.0004,0.0006,0.0008,0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,1};

	double rpr = 0.0003, *NOut; // basal rate and TF production rate
	double gmt = 0.011, gpt = 0.01, gm2 = 0.11, gp2 = 0.09; // gmt = transcription for TF mrna, gm2 for target
	double sig1, sig2, TF = 200;
	FILE *file;

	file = fopen("AutoRegulated/Auto_VaryTFAffinity_Ku2-0.0002_Kud-0.0010_TF-200-1.txt","w");
	
	for (int i=0;i<5;i++){
			ku2 = 0.0002;
			ku1 = 0.04+i*0.02;kuv[i];
			Decoy = 0;
			//gmt =0.01;
			sig1 = 0.0027/(ku1 + rpr);
			sig2 = 0.0027/(ku2 + rpr);
			gpt = TF*rpr*0.011*(1+TF*sig1)*(1 + sig1/(1+TF*sig1) + sig2/(1+TF*sig2))/gmt; // keep TF number fixed by gmt

			NOut = run_sim(ku1,ku2,kud,rpr,Decoy,gmt,gpt,gm2,gp2);
			fprintf(file,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",ku1,ku2,gmt,gpt,NOut[0],NOut[1],NOut[2],NOut[3],NOut[4],NOut[5],NOut[6],NOut[7],NOut[8],NOut[9],NOut[10],NOut[11]);
			printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",ku1,ku2,gmt,gpt,NOut[0],NOut[1],NOut[2],NOut[3],NOut[4],NOut[5],NOut[6],NOut[7],NOut[8],NOut[9],NOut[10],NOut[11]);
	}
	fclose(file);
	tfinish = time (NULL);
	long comp_time = (tfinish - tstart);
	printf("Computation time = %ld, \n", comp_time);

}
