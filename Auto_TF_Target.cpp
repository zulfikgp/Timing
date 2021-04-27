// Compute steady state protein expression and mean first passage time to reach a given threshold (Half maxima of steady stae)
// for an autoregulated TF gene and a target gene in presence of decoy sites
// compile and run the code, using something like: g++ -c Auto_TF_Ony.cpp, g++ -o exec_file Auto_TF_Ony.o, ./exec_file
// To edit the parameter values, go to the last section of the code
// For identical TF and target gene set gm01 = gm02, fac = fac2, gp1 = gp2, ku1 = ku2
// In the code, regulatory strength (fac, alpha in the paper) for the TF gene is varied
// Alternatively, one can vary binding affinities of TF gene

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
int* Gillespiess(int Decoy, double tend, int Nm1, int Nm2, int Np1, int Np2, int comp1, int comp2, int CDecoy, double ku1, double ku2, double kud, double gm0, double gm02, double gp1, double gp2, double gm1, double gm2);
int* Gillespiess(int Decoy, double tend, int Nm1, int Nm2, int Np1, int Np2, int comp1, int comp2, int CDecoy, double ku1, double ku2, double kud, double gm0, double gm02, double gp1, double gp2, double gm1, double gm2)
{	
	int gene1 = 1, gene2 = 1;
// Rates

	double rp1 = 0.0003; //Degradation of TF protein
	double rp2 = rp1; //Degradation of target protein

	double kb1 = 0.0027, kb2 = 0.0027, kbd = 0.0027; // Binding first gene
	//double gm2 = gm02*fac2; // mRNA production target gene
	//double gm1 = gm0*fac; // mRNA production	
	double rm1 = 0.011, rm2 = rm1; // Degradation of RNA

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
	Rates[7] = gm2*comp2 + gm02*(gene2-comp2); // mRNA production gene2

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
		Rates[0] = kb1*(gene1-comp1)*Np1;
		Rates[1] = kb2*(gene2-comp2)*Np1;
		Rates[2] = ku1*comp1; 
		Rates[3] = ku2*comp2; 

		Rates[4] = kbd*(Decoy-CDecoy)*Np1;
		Rates[5] = kud*CDecoy;

		Rates[6] = gm1*comp1 + gm0*(gene1-comp1); 
		Rates[7] = gm2*comp2 + gm02*(gene2-comp2);

		Rates[8] = gp1*Nm1;
		Rates[9] = gp2*Nm2;

		Rates[10] = rm1*Nm1;
		Rates[11] = rm2*Nm2;

		Rates[12] = rp1*Np1;
		Rates[13] = rp2*Np2;

		Rates[14] = rp1*comp1;
		Rates[15] = rp1*comp2;
		Rates[16] = rp1*CDecoy;
	}
	return NOut;
}

/****************************/

double* Gillespie(int Decoy, double Tst, double Tsf, double Tss, double tend, double ku1, double ku2, double kud, double gm0, double gm02, double gp1, double gp2, double gm1, double gm2);
double* Gillespie(int Decoy, double Tst, double Tsf, double Tss, double tend, double ku1, double ku2, double kud, double gm0, double gm02, double gp1, double gp2, double gm1, double gm2)
{	
// Gene1 TF complex, mRNA1, mRNA2, Protein1, Protein2, Decoy complex
	int comp1 = 0, comp2 = 0, Nm1 = 0, Nm2 = 0, Np1 = 0, Np2 = 0, CDecoy = 0; 
	int gene1 = 1, gene2 = 1;
// Rates

	double rp1 = 0.0003, rp2 = rp1;; //Degradation of TF protein and target protein
	double kb1 = 0.0027, kb2 = 0.0027, kbd = 0.0027; // TF Binding  TF-gene, target gene and 
//	double gm1 = gm0*fac, gm2 = gm02*fac2; // mRNA production TF, target
	double rm1 = 0.011, rm2 = rm1; // Degradation of RNA

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

	Rates[6] = gm1*comp1 + gm0*(gene1-comp1); // mRNA production TF
	Rates[7] = gm2*comp2 + gm02*(gene2-comp2);// mRNA production target

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
		Rates[0] = kb1*(gene1-comp1)*Np1;
		Rates[1] = kb2*(gene2-comp2)*Np1;
		Rates[2] = ku1*comp1; 
		Rates[3] = ku2*comp2; 

		Rates[4] = kbd*(Decoy-CDecoy)*Np1; 
		Rates[5] = kud*CDecoy;

		Rates[6] = gm1*comp1 + gm0*(gene1-comp1); 
		Rates[7] = gm2*comp2 + gm02*(gene2-comp2);

		Rates[8] = gp1*Nm1;
		Rates[9] = gp2*Nm2;

		Rates[10] = rm1*Nm1;
		Rates[11] = rm2*Nm2;

		Rates[12] = rp1*Np1;
		Rates[13] = rp2*Np2;

		Rates[14] = rp1*comp1;
		Rates[15] = rp1*comp2;
		Rates[16] = rp1*CDecoy;

		}
	double *TOut;
	TOut = dvmalloc(3);
	TOut[0] = t1; // Response time for Total TF
	TOut[1] = t2; // Response time for Free TF
	TOut[2] = t3; // Response time for Target
	return TOut;
}

/********* MAIN CODE ***********/

/* For identical TF and target gene set gm01 = gm02, fac = fac2, gp1 = gp2, ku1 = ku2 */
/* Also, comment the section where translation rate is tuned for fixed TF number*/


int main()
{
	srand(time(NULL));
	time_t tstart, tfinish;
	tstart = time (NULL);
	FILE *file, *filer;

	file = fopen("Auto_TF_target.txt","w"); // File to write the output

	double ku1 = 0.002, ku2 = 0.0002, kud = 0.0024; // Unbinding rates for TF gene, target gene and decoy.

	int Decoy = 0, Nin[7], *NOut, *Ntemp;;
	double Tss = 0, Tst = 0, Occ1, Occ2, Occ3, Tsf = 0;

	int Nit = 1*pow(10,5), Nit2 = 1*pow(10,5); // Number of iterations to compute steady state distribution and MFPT
	double tend = 1*0.5*pow(10,5); // Time step for sampling over time 
	double *TRes;
	int TF = 50; // # of TF proteins when a fixed TF is desired for varying parameters 


	double fac = 0.0, fac2 = 0.0; // fac > 1 (auto-activator), <1 auto-repressor, fac2 > 1(target activation) < 1 (target repression)
	double gm0 = 0.025, gm1 = fac*gm0, gm02, gm2; // transcription rates
	double gp1 = 0.1, gp2 = 0.1; // translation rates for TF and target

	/* For pure target repression (zero transcription when TF bound) set gm2 = 0 (fac2=0)*/
	/* For pure target activation (zero transcription when TF free) set gm02 = 0  and gm2 = value */
	gm02 = 0.025;
	gm2 = fac2*gm02;


	double sig1, sig2, theta;
	double fccc[22] = {0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	double KUV[16] = {0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.15, 0.2, 0.3, 0.4,  0.6, 0.8, 1};
	double avgt1, avgt2, avgt3, avgt12,avgt22, avgt33, Var1, Var2, Var3;

	for (int i=0;i<22;i++){	
		Decoy = 0;
//		ku2 = 0.0002;
//		ku1 = 0.002; // for identical sites set ku1 = ku2
//		ku2 = KUV[i];
		fac = fccc[i];
		gm1 = fac*gm0;
		avgt1=0, avgt12=0, avgt2 = 0, avgt3 = 0, avgt22 = 0, avgt33=0, Tss = 0, Tst = 0, Tsf = 0, Occ1 = 0, Occ2 = 0, Occ3 = 0;

/*************** Tune transcription or translation rate to get fixed number of TF based on mass action kinetics **************/
/********* Might have to use different rates in some cases as deterministic and stocahstic values differs  a lot ****/
/********* comment this section for identical TF and target *****/
		sig1 = 0.0027/(0.0003 + ku1);
		sig2 = 0.0027/(0.0003 + ku2);
		theta = 1/fac; // ratio of basal and activated rate of transcription
		//gm1 = TF*0.0003*0.011*(1+TF*sig1)*(1 + sig1/(1+TF*sig1) + sig2/(1+TF*sig2)) / (gp1*(theta + TF*sig1)); // tune transcription
		gp1 = TF*0.0003*0.011*(1+TF*sig1)*(1 + sig1/(1+TF*sig1) + sig2/(1+TF*sig2)) / ((gm0 + TF*sig1*gm1)); // Tune translation
/*****************************************************************/


		/*Find steady state values for TF and target expressions*/
		Ntemp = Gillespiess(Decoy, pow(10,6), 0, 0, 0, 0, 0, 0, 0,ku1,ku2,kud,gm0,gm02,gp1,gp2,gm1,gm2);
		Nin[0] = Ntemp[0], Nin[1] = Ntemp[1], Nin[2] = Ntemp[2], Nin[3] = Ntemp[3], Nin[4] = Ntemp[4], Nin[5] = Ntemp[5], Nin[6] = Ntemp[6];
		free(Ntemp);
		for (int ni=0;ni<Nit2;ni++){
			NOut = Gillespiess(Decoy, tend, Nin[0], Nin[1], Nin[2], Nin[3], Nin[4], Nin[5], Nin[6],ku1,ku2,kud,gm0,gm02,gp1,gp2,gm1,gm2);
			Nin[0] = NOut[0], Nin[1] = NOut[1], Nin[2] = NOut[2], Nin[3] = NOut[3], Nin[4] = NOut[4], Nin[5] = NOut[5], Nin[6] = NOut[6];
			Tst = Tst + NOut[2] + NOut[4] + NOut[5] + NOut[6]; // total TF 
			Tss = Tss + NOut[3]; // Target expression
			Tsf = Tsf + NOut[2]; // free TF expression
			Occ1 = Occ1 + NOut[4]; // TF gene occupancy
			Occ2 = Occ2 + NOut[5]; // Target gene occupancy
			Occ3 = Occ3 + NOut[6]; // Decoy occupancy
			free(NOut);
		}
		Tst = Tst/Nit2; // Average steady state total TF
		Tsf = Tsf/Nit2; // Average steady state free TF
		Tss = Tss/Nit2; // average target expression
		Occ1 = Occ1/Nit2; 
		Occ2 = Occ2/Nit2;
		Occ3 = Occ3/Nit2;
		//printf("%f %f\n",Tst,Tss);
		/************** */
	
		/*Find first passage time for TF and target genes*/
		for (int ni=0;ni<Nit;ni++){
			TRes = Gillespie(Decoy, Tst, Tsf, Tss, tend, ku1, ku2, kud, gm0, gm02, gp1, gp2, gm1, gm2);
			avgt1 = avgt1 + TRes[0]; //sum of FPT for Total TF
			avgt2 = avgt2 + TRes[1]; //sum of FPT for Free TF
			avgt3 = avgt3 + TRes[2]; //sum of FPT for Target
			avgt12 = avgt12 + TRes[0]*TRes[0];
			avgt22 = avgt22 + TRes[1]*TRes[1];
			avgt33 = avgt33 + TRes[2]*TRes[2];
		}
		avgt1 = avgt1/Nit; // MFPT for Total TF
		avgt2 = avgt2/Nit; // MFPT for free TF
		avgt3 = avgt3/Nit; // MFPT for Target
		Var1 = avgt12/Nit-pow(avgt1/Nit,2); // Variance in FPT total TF
		Var2 = avgt22/Nit-pow(avgt2/Nit,2); // Variance in FPT free TF
		Var3 = avgt33/Nit-pow(avgt3/Nit,2); // Variance in FPT target
	
		fprintf(file,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",fac,gm0,gp1,ku1,ku2,avgt1,avgt2,avgt3,Var1,Var2,Var3,Tst,Tsf,Tss,Occ1,Occ2,Occ3);
		//printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",fac,gm0,gp1,ku1,ku2,avgt1,avgt2,avgt3,Var1,Var2,Var3,Tst,Tsf,Tss,Occ1,Occ2,Occ3);
		tfinish = time (NULL);
		long comp_time = (tfinish - tstart);
		printf("Computation time = %ld\n", comp_time);
	}
	fclose(file);
}
