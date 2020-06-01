//*******************************************************************
//* 3D (Re)Construction from Radial S2(r) using Orthogonal Sampling *
//*******************************************************************


//ORTHOGONAL SAMPLING is used to generate large systems
//Both lattice-point and gray-scale method has limitations on system size
//Sometimes, small systems can not accurately reproduce the correlations...



//modified: 05/13/2012
//using a more efficient sampling method for computing S_2 for each line
// i.e., instead of actually keeping track of line segments moving through 
// the system, we focus on the black pixels as particles and directly
// compute the pair distances and then normalize them 
// This is like the orthgoanl version of the lattice-point method



using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>


#define MAXX 192
#define MAXY 250000

#define Nt MAXX/2


double f1;//volume fraction of black pixels
double f2;
long int NP;

int indexi; int indexj; int indexm; int indexn;
int indexk; int indexl;
//int cindexi; int cindexj; int cindexm; int  cindexn;
//in case to resume the configuration...


//the following is to treat pix as particles
int pix_position[MAXX]; //the position of black pixels
int pix_counter; //the number of pix on each line/column/height...
int temp_npL = 0; 
int temp_npC = 0;
int temp_npH = 0;


//the 3D configuration...
int config[MAXX][MAXX][MAXX];
int best_config[MAXX][MAXX][MAXX];

//need to index to pin the 1D array
long int lineS2[MAXX][MAXX][Nt];
long int columeS2[MAXX][MAXX][Nt];
long int heightS2[MAXX][MAXX][Nt];

double obj[Nt];


//**************************************
//*** The cooling schedule

double alpha = 0.85;
//double beta = 0.85;

int TN = 150; // # of times lower down the temperature
double T = 0.00001;

//int Nconf = 1;
int Nevl = 120000;//times of evlovtion, on average, each pixel is moved 0.5 times... far form equilirum at each stage... but just for test...

int flag_iconfig; //indicating whether read in a configuration (1) or starting from random (0)

//**************************************

void read_parameter()
{
  cout<<"Reading parameters for annealing reconstruction from standard input"<<endl;

  cout<<"Init config flag flag_iconfig = "; cin>>flag_iconfig; cout<<endl;

  cout<<"Number of black pixels: NP = "; cin>>NP; cout<<endl;

  //now computing the volume fraction...
  f1 = (double)NP/(double)(MAXX*MAXX*MAXX);
  f2 = 1.0 - f1;

  //now for the cooling schedule
  cout<<"Starint temp T0 = "; cin>>T; cout<<endl;
  cout<<"Decreasing ratio: alpha = "; cin>>alpha; cout<<endl;
  cout<<"Number of decreasing T stages: TN = "; cin>>TN; cout<<endl;
  cout<<"Number of pxiel move per stage: Nevl = "; cin>>Nevl; cout<<endl;


  
}

void get_obj()
{
  
  FILE* fp;

  if((fp = fopen("sobj.txt","r"))==NULL)
    {
      printf("Can not open file sobj.txt! Abort!\n");
      exit(1);
    }
  
  int ind;
  float value;

  //printf("here\n");

  fscanf(fp, "%d", &ind);
  fscanf(fp, "%f", &value);

  //printf("here\n");

  double fr = f1;

  if((1.0-value)<0.001)
    {
      obj[0] = fr*(1-fr)*value + fr*fr;
      for(int i=1; i<Nt; i++)
	{
	  fscanf(fp, "%d", &ind);
	  fscanf(fp, "%f", &value);

	   obj[i] = fr*(1-fr)*value + fr*fr;
	}

      //printf("here\n");
    }
  else
    {
      obj[0] = value;
      for(int i=1; i<Nt; i++)
	{
	  fscanf(fp, "%d", &ind);
	  fscanf(fp, "%f", &value);

	  
	  obj[i] = value;
	}

      //printf("here\n");
    }

  fclose(fp);
 

  //printf("f = %f\n", fr);
  //printf("obj[0] = %f \n", obj[0]);
  //printf("obj[1] = %f \n", obj[1]);

  fp = fopen("obj.txt","w");
   for(int r=0; r<Nt; r++)
     fprintf(fp,"%d \t %f\n", r, obj[r]);
  fclose(fp);  
}


void init_configII()
{
  srand(time(NULL));

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      for(int k=0; k<MAXX; k++)
      {
        config[i][j][k] = 0;
      }


  for(int i=0; i<NP; i++)
    {
      int m = rand()%MAXX;
      int n = rand()%MAXX;
      int l = rand()%MAXX;

      while(config[m][n][l]==1)
	{
	  m = rand()%MAXX;
	  n = rand()%MAXX;
	  l = rand()%MAXX;
	}

      config[m][n][l] = 1;
    }

  //****************************
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      for(int k=0; k<MAXX; k++)
         best_config[i][j][k] = 0;

}



void read_config()
{
  FILE* fp;

  if((fp = fopen("Mconfig.txt","r"))==NULL)
    {
      printf("Can not open file Mconfig.txt! Abort!\n");
      exit(1);
    }


  int xt;
  int yt; 
  int zt;

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      for(int k=0; k<MAXX; k++)
      {
        config[i][j][k] = 0;
      }

  for(int i=0; i<NP; i++)
    {
      fscanf(fp, "%d", &xt);
      fscanf(fp, "%d", &yt);
      fscanf(fp, "%d", &zt);

      config[xt][yt][zt] = 1;
    }

  fclose(fp);

}



int abs(int a)
{
  if(a>=0) return a;
  else return -a;
}


void sampleS2line(int index1, int index2)
{
  
  for(int r=0; r<Nt; r++)
    {
      lineS2[index1][index2][r] = 0;
    }

  //serach the line for pixel positions
  pix_counter = 0;

  for(int i=0; i<MAXX; i++)
    {
      if(config[i][index1][index2] == 1)
	{
	  pix_position[pix_counter] = i;
	  pix_counter++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position[i]-pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	lineS2[index1][index2][temp_dist]++;

	//if(temp_dist == 0) temp_np++;
      }

  //temp_npL += lineS2[index1][index2][0];
  

  //cout<<"s2Line_Np = "<<lineS2[index1][index2][0]<<endl;
}



void sampleS2colume(int index1, int index2)
{

  for(int r=0; r<Nt; r++)
    {
      columeS2[index1][index2][r] = 0;
    }
  
  //serach the line for pixel positions
  pix_counter = 0;

  for(int i=0; i<MAXX; i++)
    {
      if(config[index1][i][index2]==1)
	{
	  pix_position[pix_counter] = i;
	  pix_counter++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position[i]-pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	columeS2[index1][index2][temp_dist]++;

	//if(temp_dist == 0) temp_np++;
	
      }

  // temp_npC += columeS2[index1][index2][0];
}


void sampleS2height(int index1, int index2)
{

  for(int r=0; r<Nt; r++)
    {
      heightS2[index1][index2][r] = 0;
    }
  
  //serach the line for pixel positions
  pix_counter = 0;
  
  for(int i=0; i<MAXX; i++)
    {
      if(config[index1][index2][i]==1)
	{
	  pix_position[pix_counter] = i;
	  pix_counter++;
	}
    }
  
  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position[i]-pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	heightS2[index1][index2][temp_dist]++;

	
      }
  
  //temp_npH += heightS2[index1][index2][0];
  
}



void change_config()
{
  int i, j, k, m, n, l;
  int lim = 0;
  
 
  do{   i = rand() % MAXX;
        m = rand() % MAXX;

        j = rand() % MAXX; 
        n = rand() % MAXX;

        k = rand() % MAXX; 
        l = rand() % MAXX;

        lim++;

     } while(config[i][j][k] == config[m][n][l] && lim < 1000);


  int temp;
  temp = config[i][j][k];
  config[i][j][k] = config[m][n][l];
  config[m][n][l] = temp;

  indexi = i;
  indexj = j;
  indexk = k;
  indexm = m;
  indexn = n;
  indexl = l;

}

/*
//this is just for testing 
void change_config_fake()
{
  int i, j, k, m, n, l;
  int lim = 0;
  
 
  do{   i = rand() % MAXX;
        m = rand() % MAXX;

        j = rand() % MAXX; 
        n = rand() % MAXX;

        k = rand() % MAXX; 
        l = rand() % MAXX;

        lim++;

     } while(config[i][j][k] == config[m][n][l] && lim < 1000);


  //int temp;
  //temp = config[i][j][k];
  //config[i][j][k] = config[m][n][l];
  //config[m][n][l] = temp;

  indexi = i;
  indexj = j;
  indexk = k;
  indexm = m;
  indexn = n;
  indexl = l;

}
*/

void resume_config()
{
  int temp;
  //first we resume the config
  temp = config[indexi][indexj][indexk];
  config[indexi][indexj][indexk] = config[indexm][indexn][indexl];
  config[indexm][indexn][indexl] = temp;

}


double energy(double S[Nt])
{
  double E=0;

  for(int i=1; i<Nt; i++)
    {
      E = E + (S[i] - obj[i])*(S[i] - obj[i]);
    }
  
  
  return E;
}


double d_energy(double S2[Nt], double ST2[Nt])
{
  double d_E = 0;
  d_E = energy(ST2) - energy(S2);
  return d_E;
}

double PE(double dE, double T) // the probability that should compare ...
{
  if(dE > 0) return exp(-dE/T);
  else return 1;
}


main()
{
  double S1 = 0;

  double S2[Nt];
  double ST2[Nt];

  long int  SS2[Nt];

  long int SLT[2][Nt];
  long int SCT[2][Nt];
  long int SHT[2][Nt];//save the current values for changed lines, columes and heights...


  double energyb = MAXY;
  double energyt = 0;


  for(int i=0; i<Nt; i++)
    {
      S2[i] = 0;
      ST2[i] = 0;
      SS2[i] = 0;
    }


  get_obj();
 
  //read_config(); 
  //init_configII();//initialize configuration. volume fraction preserved...

  read_parameter();
  
  if(flag_iconfig == 0)
    init_configII();//initialize configuration. volume fraction preserved...
  else
    read_config();



  //now we sample S2 for the first time...

  
  cout<<"initial sampling S2...."<<endl;
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      {
	sampleS2line(i,j);
	sampleS2colume(i,j);
	sampleS2height(i,j);
	
        temp_npL += lineS2[i][j][0];
	temp_npC += columeS2[i][j][0];
	temp_npH += heightS2[i][j][0];

	for(int r=0; r<Nt; r++)
	  {
	    SS2[r] += lineS2[i][j][r];
	    SS2[r] += columeS2[i][j][r];
	    SS2[r] += heightS2[i][j][r];
	  }


	
     } 

  //cout<<"temp_npL = "<<temp_npL<<endl;
  //cout<<"temp_npC = "<<temp_npC<<endl;
  //cout<<"temp_npH = "<<temp_npH<<endl;

  for(int r=0; r<Nt; r++)
    {
      S2[r] = (double) SS2[r]/(double)(3*MAXX*MAXX*MAXX);

      printf("%d \t  %f \n", r, S2[r]);
    }
  


  /*
  for(int r=0; r<Nt; r++)
    {
      temp_npL = 0;
      temp_npC = 0;
      temp_npH = 0;

      for(int i=0; i<MAXX; i++)
	for(int j=0; j<MAXX; j++)
	  {
	    //SS2[r] += lineS2[i][j][r];
	    //SS2[r] += columeS2[i][j][r];
	    //SS2[r] += heightS2[i][j][r];
	    
	    temp_npL += lineS2[i][j][r];
	    temp_npC += columeS2[i][j][r];
	    temp_npH += heightS2[i][j][r];
	    
	  }

      cout<<"temp_npL = "<<temp_npL<<endl;
      cout<<"temp_npC = "<<temp_npC<<endl;
      cout<<"temp_npH = "<<temp_npH<<endl;

      
      SS2[r] = temp_npL + temp_npC + temp_npH;
      
      S2[r] = (double) SS2[r]/(double)(3*MAXX*MAXX*MAXX);
      
      printf("%d \t  %f \n", r, S2[r]);

    }
  */

  // cout<<"temp_npL = "<<temp_npL<<endl;
  // cout<<"temp_npC = "<<temp_npC<<endl;
  // cout<<"temp_npH = "<<temp_npH<<endl;
  
   cout<<"*****************************************"<<endl;

   FILE* fp = fopen("TS2.txt", "w");
   for(int r=0; r<Nt; r++)
     {
       fprintf(fp, "%d \t  %lf \n", r, S2[r]);
     }
   fclose(fp);
   
   

   //simulated annealing procedure to evlove the system
   //*****************************************************************
   //*****************************************************************
   
   int N_acc = 0; //acceptance rate...

   //these are for updating SS2 efficiently...
   int TEMP[Nt];
   int LCHV[Nt];

   cout<<"Staring the simulated annealing reconstruction process..."<<endl;
   for(int q=0; q<TN; q++)
     {
       T = alpha*T;

       N_acc = 0;

       cout<<"Stage "<<q+1<<" with T = "<<T<<endl;
       
       for(int i=0; i< Nevl; i++)
	 {
	   change_config();
	   //sample S2 for the new configuration, using time saving methods
	   
	   
	   //first save the values of lines and columes that will be changed
	   for(int r=0; r<Nt; r++)
	     {
	       SLT[0][r] = lineS2[indexj][indexk][r];
	       SLT[1][r] = lineS2[indexn][indexl][r];
	       
	       
	       SCT[0][r] = columeS2[indexi][indexk][r];
	       SCT[1][r] = columeS2[indexm][indexl][r];
	       
	       
	       SHT[0][r] = heightS2[indexi][indexj][r]; 
	       SHT[1][r] = heightS2[indexm][indexn][r];
	     }
	   
	   
	   sampleS2line(indexj, indexk);
	   sampleS2line(indexn, indexl);	     
	   //only sample the lines that changed...new results are stored in lineS2[lindexi,lindexm][*] 
	   
	   sampleS2colume(indexi, indexk);
	   sampleS2colume(indexm, indexl);
	   //only sample the columes that changed...new results are stored in columeS2[cindexi,cindexm][*]
	   
	   sampleS2height(indexi, indexj);
	   sampleS2height(indexm, indexn);
	   
	   
	   //Now we compute the S2 for the new configuration...
	   for(int r=0; r<Nt; r++)
	     {       
	       //the following method only consider the changes.. 
	       //************************************************
	       TEMP[r] = 0; //old value 
	       LCHV[r] = 0; //current value
	       
	       for(int vv=0; vv<2; vv++)
		 {
		   TEMP[r] += (SLT[vv][r]+SCT[vv][r]+SHT[vv][r]);
		 }
	       
	       LCHV[r] = (lineS2[indexj][indexk][r]+lineS2[indexn][indexl][r])+
		 (columeS2[indexi][indexk][r]+columeS2[indexm][indexl][r])+
		 (heightS2[indexi][indexj][r]+heightS2[indexm][indexn][r]);

	       /*
	       if(r!=0)
		 {
		   if(TEMP != LCHV)
		     {
		       cout<<"TEMP = "<<TEMP<<endl;
		       cout<<"LCHV = "<<LCHV<<endl;
		     }
		 }
	       */
	       
	       SS2[r] = (int)(SS2[r] - TEMP[r] + LCHV[r]);
	       
	       ST2[r] = (double)SS2[r]/(double)(3*MAXX*MAXX*MAXX);

	       //cout<<"ST2[r] = "<<ST2[r]<<endl;
	       
	       //************************************************
	     }
	   
	   //Monte Carlo steps...

	   //cout<<"dE = "<<d_energy(S2, ST2)<<endl;
	   
	   double P = double (rand() % MAXY)/(double) MAXY;
	   
	   if( P > PE(d_energy(S2, ST2), T))
	     {
	       resume_config();
	       //this just resumes the 'configuration', still need to resume S2...
	       
	       
	       for(int r=0; r<Nt; r++)
		 {
		   /*
		   TEMP = 0;
		   LCHV = 0;
		   
		   for(int vv=0; vv<2; vv++)
		     {
		       TEMP += (SLT[vv][r]+SCT[vv][r]+SHT[vv][r]);
		     }
		   
		   LCHV = (lineS2[indexj][indexk][r]+lineS2[indexn][indexl][r])+
		     (columeS2[indexi][indexk][r]+columeS2[indexm][indexl][r])+
		     (heightS2[indexi][indexj][r]+heightS2[indexm][indexn][r]);
		   */		   

		   
		   SS2[r] = (int)(SS2[r] + TEMP[r] - LCHV[r]);
		  
		   
		   //S2 does not change, do not need to resume that...
		   //now we resume lineS2, columeS2 and heightS2;
		   lineS2[indexj][indexk][r] = SLT[0][r];
		   lineS2[indexn][indexl][r] = SLT[1][r];
		   
		   columeS2[indexi][indexk][r] = SCT[0][r];
		   columeS2[indexm][indexl][r] = SCT[1][r];
		   
		   heightS2[indexi][indexj][r] = SHT[0][r];
		   heightS2[indexm][indexn][r] = SHT[1][r];
		      
		 }
	     } 
	   
	   else 
	     {
	       for(int r=0; r<Nt; r++)
		 {
		   S2[r] = ST2[r];
		 }

	       N_acc ++;
	     }
	   

	   //compare and record the best energy and configuration...
	   energyt = energy(S2);
	   
	   
	   if(energyt < energyb)
	     {
	       energyb = energyt;
	       
	       /*
	       for(int i=0; i<MAXX; i++)
		 for(int j=0; j<MAXX; j++)
		   for(int k=0; k<MAXX; k++)
		     {
		       best_config[i][j][k] = config[i][j][k];
		     }
	       */
	     }
	   
	   //printf("%f   %d change has finished... \n", energyt, i+1 );
	   
	 }
       printf("%d th change of temperature has finished... \n",q+1 );
       cout<<"The acceptance rate: "<<(double)N_acc/(double)Nevl<<endl;
       cout<<"The energy E = "<<energyb<<endl;
       
       
       fp = fopen("E.txt","a");
       fprintf(fp, "%e\n", energyb);
       fclose(fp);
       
       
       
       printf("*************************************************\n");
       
     }
   
   //*****************************************************************
   //*****************************************************************
   //this is the end of simulated annealing
   
   fp = fopen("S2.txt", "w");
   fprintf(fp, "%d \t %f \n", 0, f1);
   for(int r=1; r<Nt; r++)
     {
       fprintf(fp, "%d \t %f \n", r, S2[r]);
     }
   fclose(fp);

   fp = fopen("Fconfig.txt","w");
   for(int i=0; i<MAXX; i++)
     for(int j=0; j<MAXX; j++)
       for(int k=0; k<MAXX; k++)
	 {
	   if(config[i][j][k] == 1)
	     {
	       fprintf(fp, "%d \t %d \t %d \n", i, j, k );
	     }
	 }
   fclose(fp);
   
   //this is the end of the codes...
   
}
