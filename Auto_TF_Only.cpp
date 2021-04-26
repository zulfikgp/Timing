// Compute steady state protein expression and mean first passage time to reach a given threshold (Half maxima of steady stae)
// for an autoregulated TF gene in presence of decoy sites
// compile and run the code, using something like: g++ -c Auto_TF_Ony.cpp, g++ -o exec_file Auto_TF_Ony.o, ./exec_file
// To edit the parameter values, go to the last section of the code

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

/************ Function to compute Steady state protein expression****************/

int* Gillespiess(double tend, double ku1, double kud, int Decoy, double rp1, double gm0, double gp1, int Nm1, int Np1, int comp1, int CDecoy, double fac);
int* Gillespiess(double tend, double ku1, double kud, int Decoy, double rp1, double gm0, double gp1, int Nm1, int Np1, int comp1, int CDecoy, double fac)
{

// Nm1 = # of mrna, Np1 = number of free TF, comp1 = bound TF gene, CDecoy = # of bound decoy sites, Decoy = # of total decoy sites

// Rates
	int gene1 = 1;
	double gm1 = gm0*fac; // TF bound Transcription 
	double kb1 = 0.0027; // Binding first gene
	double kbd = 0.0027; // Binding/unbinding to decoy site
	double rm1 = 0.011; // Degradation of RNA

// Initiate the rates
	int NR = 10;
	double Rates[NR]; // Rate matrix
	Rates[0] = kb1*(gene1-comp1)*Np1; // TF binding to TF gene
	Rates[1] = ku1*comp1; // TF unbinding from TF gene

	Rates[2] = kbd*(Decoy-CDecoy)*Np1; // TF binding to decoy
	Rates[3] = kud*CDecoy; // TF unbinding from decoy

	Rates[4] = gm0*(gene1-comp1) + gm1*comp1; // mRNA production

	Rates[5] = gp1*Nm1; // protein production

	Rates[6] = rm1*Nm1; // mRNA degradation

	Rates[7] = rp1*Np1; // Protein degration
	
	Rates[8] = rp1*comp1; //Degradation of TFs bound to gene1
	Rates[9] = rp1*CDecoy; //Degradation of TFs bound to decoy

		
//Start Gillespie
	double *PSum, K, dt, flag, t=0;
	int Ind, Reaction, cnt1 = 0, cnt2 = 0, Ttot = 0;
	int *NOut;
	NOut = ivmalloc(4);
	 
	while (t < tend){

		NOut[0] = Nm1; NOut[1] = Np1; NOut[2] = comp1; NOut[3] = CDecoy;

		PSum = cumsum(Rates,NR); // Cumulative sum of the rates
		K = PSum[NR-1];
		dt = exprnd(K); 

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
				comp1 = comp1 - 1; // unbinding
				Np1 = Np1 + 1;
				break;
				
			case 2:
				CDecoy = CDecoy + 1; // Binding of protein to Decoy
				Np1 = Np1 - 1;
				break;

			case 3:
				CDecoy = CDecoy - 1; // Unbinding of protein from Decoy
				Np1 = Np1 + 1;
				break;
				
			case 4:
				Nm1 = Nm1 + 1; // mrna production
				break;
								
			case 5:
				Np1 = Np1 + 1; // TF production
				break;
				
			case 6:
				Nm1 = Nm1 - 1; // mrna degradation
				break;
				
			case 7:
				Np1 = Np1 - 1; // Protein degration	
				//cnt1 = cnt1 + 1;
				break;
				
			case 8:
				comp1 = comp1 - 1; // Degradation of TFs bound to gene1
				//cnt1 = cnt1 + 1;
				break;
				
			default:
				CDecoy = CDecoy - 1; // Degradation of TFs bound to decoy
				//cnt1 = cnt1 + 1;

		}

/*************** Recalculate rates ************************************/
		Rates[0] = kb1*(gene1-comp1)*Np1;
		Rates[1] = ku1*comp1;
		Rates[2] = kbd*(Decoy-CDecoy)*Np1;
		Rates[3] = kud*CDecoy;
		Rates[4] = gm0*(gene1-comp1) + gm1*comp1;
		Rates[5] = gp1*Nm1;
		Rates[6] = rm1*Nm1;
		Rates[7] = rp1*Np1;
		Rates[8] = rp1*comp1;
		Rates[9] = rp1*CDecoy;
	}

	return NOut;
}


/************** Function to Find time to first passage time to reach certain threshold *************/
/* at t=0, the expression is zero, i.e., Nm1 = 0, Np1 = 0, comp1 = 0, CDecpy = 0 */


double* Gillespie(double ku1, double kud, int Decoy, double Tst, double Tsf, double rp1, double gm0, double gp1, double fac, double Yth);
double* Gillespie(double ku1, double kud, int Decoy, double Tst, double Tsf, double rp1, double gm0, double gp1, double fac, double Yth)
{
	int gene1 = 1, comp1 = 0, Nm1 = 0, Np1 = 0, CDecoy = 0; 
	double gm1 = gm0*fac;
	
// Rates
	double kb1 = 0.0027;
	double kbd = 0.0027;
	double rm1 = 0.011; // Degradation of RNA

	double t = 0, t1 = 0, t2 = 0, t3 = 0;
	int check1 = 0, check2 = 0, check = 0,check3 = 0;


// Initiate the rates
	int NR = 10;
	double Rates[NR]; // Rate matrix

	Rates[0] = kb1*(gene1-comp1)*Np1;
	Rates[1] = ku1*comp1;

	Rates[2] = kbd*(Decoy-CDecoy)*Np1;
	Rates[3] = kud*CDecoy;

	Rates[4] = gm0*(gene1-comp1) + gm1*comp1;

	Rates[5] = gp1*Nm1;

	Rates[6] = rm1*Nm1;

	Rates[7] = rp1*Np1;
	
	Rates[8] = rp1*comp1;
	Rates[9] = rp1*CDecoy;
	
//Run SSA until the desired threshold of protein expression is achieved

	double *PSum, K, dt, flag;
	int Ind, Reaction, cnt1 = 0, cnt2 = 0, Ttot = 0;;

	while (check == 0){

		PSum = cumsum(Rates,NR);
		K = PSum[NR-1];
		dt = exprnd(K);

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
				comp1 = comp1 - 1; // unbinding
				Np1 = Np1 + 1;
				break;
				
			case 2:
				CDecoy = CDecoy + 1; // Binding of protein to Decoy
				Np1 = Np1 - 1;
				break;

			case 3:
				CDecoy = CDecoy - 1; // Unbinding of protein to Decoy
				Np1 = Np1 + 1;
				break;
				
			case 4:
				Nm1 = Nm1 + 1; // mrna production
				break;
				
			case 5:
				Np1 = Np1 + 1; // TF production
				break;
				
			case 6:
				Nm1 = Nm1 - 1; // mrna degradation
				break;
				
			case 7:
				Np1 = Np1 - 1; // Protein degration	
				//cnt1 = cnt1 + 1;
				break;
				
			case 8:
				comp1 = comp1 - 1;
				//cnt1 = cnt1 + 1;
				break;
				
			default:
				CDecoy = CDecoy - 1;
				//cnt1 = cnt1 + 1;

		}

/*************** Recalculate rates ************************************/

		Rates[0] = kb1*(gene1-comp1)*Np1;
		Rates[1] = ku1*comp1;
		Rates[2] = kbd*(Decoy-CDecoy)*Np1;
		Rates[3] = kud*CDecoy;
		Rates[4] = gm0*(gene1-comp1) + gm1*comp1; // mRNA production
		Rates[5] = gp1*Nm1; // protein production
		Rates[6] = rm1*Nm1; // mRNA degradation
		Rates[7] = rp1*Np1; // Protein degration
		Rates[8] = rp1*comp1; //Degradation of TFs bound to gene1
		Rates[9] = rp1*CDecoy; //Degradation of TFs bound to decoy

		Ttot = Np1+comp1+CDecoy; // Total TF is the sum of free, bound to TF gene and decoys

// Check f the desired threshold is achieved //

		if (Ttot >= Yth*Tst & check1 == 0){
			t1 = t;
			check1 = 1;
		}
		if (Np1 >= Yth*Tsf  & check2 == 0){
			t2 = t;
			check2 = 1;
		}
		check = check1*check2; // check = 1, means both total TF and free TF threshold is reached
	}
	
	double *NOut;
	NOut = dvmalloc(2);
	 
	NOut[0] = t1;
	NOut[1] = t2;

	return NOut;
}


/********* function to run the simulation **********************************/

double* run_sim(double ku1, double kud, double rpr, int Decoy, double gm0, double gpt, double fac, double Yth);
double* run_sim(double ku1, double kud, double rpr, int Decoy, double gm0, double gpt, double fac, double Yth)

{
	srand(time(NULL));
	time_t tstart, tfinish;
	tstart = time (NULL);
	FILE *file;
	
	double DT, *DT2;
	int Nit2 = 1*pow(10,4), Nit = 1*pow(10,4); // Number of iterations to compute steady state distribution and MFPT
	double avg1=0, avg2=0, avg12=0, avg22=0, Var1, Var2, Tst=0, Tsf = 0, Occ1 = 0, Occ2 = 0;
	int Nin[4], *NOut, *Ntemp;
	double tend = 2*0.5*pow(10,5); // Time step for sampling over time 

	//Find steady state values for TF

		//file = fopen("SSDistribution_ku-0.0010_gm0-0.005_gp-0.10_fac-10.txt","w"); // File to store data from each iteration

		Ntemp = Gillespiess(pow(10,6), ku1, kud, Decoy, rpr, gm0, gpt, 0, 0, 0, 0, fac); // Run simulation to get rid of transient
		Nin[0] = Ntemp[0], Nin[1] = Ntemp[1], Nin[2] = Ntemp[2], Nin[3] = Ntemp[3]; // Store the end point data for input
		free(Ntemp);
		for (int ni=0;ni<Nit2;ni++){
			NOut = Gillespiess(tend, ku1, kud, Decoy, rpr, gm0, gpt, Nin[0], Nin[1], Nin[2], Nin[3],fac);
			Nin[0] = NOut[0], Nin[1] = NOut[1], Nin[2] = NOut[2], Nin[3] = NOut[3];
			Tsf = Tsf + NOut[1]; // free TF
			Tst = Tst + NOut[1] + NOut[3] + NOut[2]; // Total TF
			Occ1 = Occ1 + NOut[2]; // TF Occupancy
			Occ2 = Occ2 + NOut[3]; // Decoy occupancy
			//fprintf(file,"%d\n",NOut[1] + NOut[3] + NOut[2]);
			free(NOut);
		}
		//fclose(file);
		Tst = Tst/Nit2; // Average steady state TF
		Tsf = Tsf/Nit2; // Average steady state free TF
		Occ1= Occ1/Nit2; // Average steady state TF gene occupancy
		Occ2= Occ2/Nit2;// Average steady state decoy occupancy
		//printf("%f %f\n",Tst,Tsf);
	
	/*Find Mean First passage time*/

	// file = fopen("Distribution_ku-0.0001_gm0-0.005_gp-0.10_fac-10.txt","w"); // File to store first passage time for each iteration

	for (int ni=0;ni<Nit;ni++){
		DT2 = Gillespie(ku1, kud, Decoy, Tst, Tsf, rpr, gm0, gpt, fac, Yth);
		//fprintf(file,"%f %f\n",DT2[0],DT2[1]);
		avg1 = avg1 + DT2[0];
		avg2 = avg2 + DT2[1];
		avg12 = avg12 + DT2[0]*DT2[0];
		avg22 = avg22 + DT2[1]*DT2[1];
		free(DT2);
	}
	avg1 = avg1/Nit; // MFPT for total TF
	avg2 = avg2/Nit; // MFPT for free TF
	Var1 = avg12/Nit-pow(avg1/Nit,2); // Variance in FPT for total TF
	Var2 = avg22/Nit-pow(avg2/Nit,2); // Variance in FPT for free TF
	
	double *NT;
	NT = dvmalloc(8);
	NT[0] = avg1;
	NT[1] = avg2;
	NT[2] = Var1;
	NT[3] = Var2;
	NT[4] = Tst;
	NT[5] = Tsf;
	NT[6] = Occ1;
	NT[7] = Occ2;

	tfinish = time (NULL);
	long comp_time = (tfinish - tstart);
	//printf("Computation time = %ld, \n", comp_time);
	
	return NT;
}

/********* MAIN CODE ***********/
/* Run the main script as a standard c++ code for */

int main()
{
	time_t tstart, tfinish;
	tstart = time (NULL);

	double ku1 = 0.002, kud = 0.001; // Unbinding rate for TF gene and decoy sites
	double rpr = 0.0003; // rpr = ln(2)/tau, where tau is cell division time 
	double fac = 10; // Regulatory strength, fac < 1 (auto-repression), fac > 1 (auto-activation), fac = 1 (constitutive), 
	double Yth = 0.5; // Threshold for protein expression
	double gm0 = 0.0025, gp = 0.025; // gm0 = basal transcription for TF mrna, gpt = translation rate for TF protein
	double *NOut;

	FILE *file;
	file = fopen("VaryAffinity_gm0-0.05_gp-0.1_fac-0.1.txt","w");
//	file = fopen("VaryDecoy_Ku1-0.005_gm-0.025_gp-0.025_fac-0.1.txt","w");

	int Decoy = 0; // Number of decoy sites
	double kuv[17] = {0.0001,0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,1}; // Set of TF unbinding rate
	
	for (int i=0;i<17;i++){
			ku1 = kuv[i];
			//Decoy = i*100;
			//gm0 =0.05+i*0.0;
			//gpt = 0.1;
			NOut = run_sim(ku1,kud,rpr,Decoy,gm0,gp,fac,Yth);
			fprintf(file,"%d %f %f %f %f %f %f %f %f %f %f %f %f\n",Decoy,fac,ku1,gm0,gp,NOut[0],NOut[1],NOut[2],NOut[3],NOut[4],NOut[5],NOut[6],NOut[7]);
			printf("%d %f %f %f %f %f %f %f %f %f %f %f %f\n",Decoy,fac,ku1,gm0,gp,NOut[0],NOut[1],NOut[2],NOut[3],NOut[4],NOut[5],NOut[6],NOut[7]);
	}
	fclose(file);
	tfinish = time (NULL);
	long comp_time = (tfinish - tstart);
	printf("Computation time = %ld, \n", comp_time);

}
