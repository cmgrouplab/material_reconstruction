//*******************************************************************
//* 3D (Re)Construction from Radial S2(r)+L(r) using Orthogonal Sampling *
//*******************************************************************


//ORTHOGONAL SAMPLING is used to generate large systems
//Both lattice-point and gray-scale method has limitations on system size
//Sometimes, small systems can not accurately reproduce the correlations...


//started: May 17th, 2009
//bug found on May 17th, 2009, seems to work fine now...
//It is important to note for L-sampling function used in the reconstruction, N2H, N2V, N2D need to be cleared to zero each time the sampling functions are evoked


//modified: May 11, 2012
//make a new function to read in parameters, then no need to re-compile for
//  each volume fraction with a fixed system size;
//need to recompile for different system sizes...

//modified: 05/13/2012
//using a more efficient sampling method for computing S_2 for each line
// i.e., instead of actually keeping track of line segments moving through
// the system, we focus on the black pixels as particles and directly
// compute the pair distances and then normalize them
// This is like the orthgoanl version of the lattice-point method

//modified: 05/14/2012
//L was not sampled correctly before
//missing important events when periodic boundary is used!!!

//modified 03/16/2018
//generalized to treat each direction separately


using namespace std;


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>


#define MAXX 50
#define MAXY 50
#define MAXZ 50
#define MAXS 1500000

#define NtX MAXX/2
#define NtY MAXY/2
#define NtZ MAXZ/2

double f1;//volume fraction of black pixels
double f2;
int NP;

int indexi; int indexj; int indexm; int indexn;
int indexk; int indexl;
//int cindexi; int cindexj; int cindexm; int  cindexn;
//in case to resume the configuration...

//the following is to treat pix as particles
int pix_position_X[MAXX]; //the position of black pixels
int pix_counter_X; //the number of pix on each line/column/height...
int pix_position_Y[MAXY]; //the position of black pixels
int pix_counter_Y;
int pix_position_Z[MAXZ]; //the position of black pixels
int pix_counter_Z;


//the 3D configuration...
int config[MAXX][MAXY][MAXZ];
int best_config[MAXX][MAXY][MAXZ];

//need to index to pin the 1D array
long int lineS2[MAXY][MAXZ][NtX]; // for S2(r)
long int columeS2[MAXX][MAXZ][NtY];
long int heightS2[MAXX][MAXY][NtZ];

long int N2H[MAXY][MAXZ][NtX];
long int N2V[MAXX][MAXZ][NtY]; // for L(r)
long int N2D[MAXX][MAXY][NtZ];

double objS2X[NtX];
double objS2Y[NtY];
double objS2Z[NtZ];

double objLX[NtX];
double objLY[NtY];
double objLZ[NtZ];

//**************************************
//*** The cooling schedule

double alpha = 0.90;
//double beta = 0.85;

int TN = 200; // # of times lower down the temperature
double T = 0.0000001;

//int Nconf = 1;
int Nevl = 10000000;//times of evlovtion, on average, each pixel is moved 0.5 times... far form equilirum at each stage... but just for test...

int flag_iconfig; //indicating whether read in a configuration (1) or starting from random (0)

//**************************************

void read_parameter()
{
  cout<<"Reading parameters for annealing reconstruction from standard input"<<endl;

  cout<<"Init config flag flag_iconfig = "; cin>>flag_iconfig; cout<<endl;

  cout<<"Number of black pixels: NP = "; cin>>NP; cout<<endl;

  //now computing the volume fraction...
  f1 = (double)NP/(double)(MAXX*MAXY*MAXZ);
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
  //**************************************************************
  //For L ...

  ifstream fin;


  if((fp = fopen("sobjS2.txt","r"))==NULL)
    {
      printf("Can not open file sobjS2.txt! Abort!\n");
      exit(1);
    }
  fclose(fp);

  fin.open("sobjS2.txt");


  int ind;
  float value;

  //read in B-B correlation along X, Y, Z directions...
  for(int i=0; i<NtX; i++)
    {
      //fscanf(fp, "%d", &ind);
      //fscanf(fp, "%f", &value);

      fin >> ind;
      fin >> value;

      if(i!=ind)
	{
	  cout<<"index is not right when reading objS2X, abort!"<<endl;
	  exit(1);
	}

      objS2X[i] = value;
    }
  for(int i=0; i<NtY; i++)
    {
      //fscanf(fp, "%d", &ind);
      //fscanf(fp, "%f", &value);

      fin >> ind;
      fin >> value;

      if(i!=ind)
	{
	  cout<<"index is not right when reading objS2Y, abort!"<<endl;
	  exit(1);
	}


      objS2Y[i] = value;
    }
  for(int i=0; i<NtZ; i++)
    {
      //fscanf(fp, "%d", &ind);
      //fscanf(fp, "%f", &value);

      fin >> ind;
      fin >> value;

      if(i!=ind)
	{
	  cout<<"index is not right when reading objS2Z, abort!"<<endl;
	  exit(1);
	}


      objS2Z[i] = value;
    }

  //fclose(fp);
  fin.close();

  //print out the read-in ...
  fp = fopen("objS2.txt","w");
  for(int r=0; r<NtX; r++)
    fprintf(fp,"%d \t %f\n", r, objS2X[r]);
  fprintf(fp, "#######################\n");
  for(int r=0; r<NtY; r++)
    fprintf(fp,"%d \t %f\n", r, objS2Y[r]);
  fprintf(fp, "#######################\n");
  for(int r=0; r<NtZ; r++)
    fprintf(fp,"%d \t %f\n", r, objS2Z[r]);
  fprintf(fp, "#######################\n");
  fclose(fp);

  //*********************************************************
  //For L ...

   if((fp = fopen("sobjL.txt","r"))==NULL)
    {
      printf("Can not open file sobjL.txt! Abort!\n");
      exit(1);
    }
  fclose(fp);

  fin.open("sobjL.txt");


//  int ind;
 // float value;

  //read in B-B correlation along X, Y, Z directions...
  for(int i=0; i<NtX; i++)
    {
      //fscanf(fp, "%d", &ind);
      //fscanf(fp, "%f", &value);

      fin >> ind;
      fin >> value;

      if(i!=ind)
	{
	  cout<<"index is not right when reading objLX, abort!"<<endl;
	  exit(1);
	}

      objLX[i] = value;
    }
  for(int i=0; i<NtY; i++)
    {
      //fscanf(fp, "%d", &ind);
      //fscanf(fp, "%f", &value);

      fin >> ind;
      fin >> value;

      if(i!=ind)
	{
	  cout<<"index is not right when reading objLY, abort!"<<endl;
	  exit(1);
	}


      objLY[i] = value;
    }
  for(int i=0; i<NtZ; i++)
    {
      //fscanf(fp, "%d", &ind);
      //fscanf(fp, "%f", &value);

      fin >> ind;
      fin >> value;

      if(i!=ind)
	{
	  cout<<"index is not right when reading objLZ, abort!"<<endl;
	  exit(1);
	}


      objLZ[i] = value;
    }

  //fclose(fp);
  fin.close();


  //print out the read-in ...

   fp = fopen("objL.txt","w");
  for(int r=0; r<NtX; r++)
    fprintf(fp,"%d \t %f\n", r, objLX[r]);
  fprintf(fp, "#######################\n");
  for(int r=0; r<NtY; r++)
    fprintf(fp,"%d \t %f\n", r, objLY[r]);
  fprintf(fp, "#######################\n");
  for(int r=0; r<NtZ; r++)
    fprintf(fp,"%d \t %f\n", r, objLZ[r]);
  fprintf(fp, "#######################\n");
  fclose(fp);

}


void read_config()
{
  FILE* fp;

  if((fp = fopen("Mconfig.txt","r")) == NULL)
    {
      printf("Can not open target Mconfig file! Abort!\n");
      exit(1);
    }

  int xt;
  int yt;
  int zt;

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
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



void init_configII()
{
  srand(time(NULL));

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
      {
        config[i][j][k] = 0;
      }


  for(int i=0; i<NP; i++)
    {
      int m = rand()%MAXX;
      int n = rand()%MAXY;
      int l = rand()%MAXZ;

      while(config[m][n][l]==1)
	{
	  m = rand()%MAXX;
	  n = rand()%MAXY;
	  l = rand()%MAXZ;
	}

      config[m][n][l] = 1;
    }

  //****************************
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
         best_config[i][j][k] = 0;

}

int abs(int a)
{
  if(a>=0) return a;
  else return -a;
}


void sampleS2line(int index1, int index2)
{

  for(int r=0; r<NtX; r++)
    {
      lineS2[index1][index2][r] = 0;
    }

  //serach the line for pixel positions
  pix_counter_X = 0;

  for(int i=0; i<MAXX; i++)
    {
      if(config[i][index1][index2] == 1)
	{
	  pix_position_X[pix_counter_X] = i;
	  pix_counter_X++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter_X; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position_X[i]-pix_position_X[j]);

	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;

	lineS2[index1][index2][temp_dist]++;

	//if(temp_dist == 0) temp_np++;
      }

  //temp_npL += lineS2[index1][index2][0];


  //cout<<"s2Line_Np = "<<lineS2[index1][index2][0]<<endl;
}



void sampleS2colume(int index1, int index2)
{

  for(int r=0; r<NtY; r++)
    {
      columeS2[index1][index2][r] = 0;
    }

  //serach the line for pixel positions
  pix_counter_Y = 0;

  for(int i=0; i<MAXY; i++)
    {
      if(config[index1][i][index2]==1)
	{
	  pix_position_Y[pix_counter_Y] = i;
	  pix_counter_Y++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter_Y; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position_Y[i]-pix_position_Y[j]);

	if(temp_dist>=MAXY/2) temp_dist = MAXY-temp_dist;

	columeS2[index1][index2][temp_dist]++;

	//if(temp_dist == 0) temp_np++;

      }

  // temp_npC += columeS2[index1][index2][0];
}


void sampleS2height(int index1, int index2)
{

  for(int r=0; r<NtZ; r++)
    {
      heightS2[index1][index2][r] = 0;
    }

  //serach the line for pixel positions
  pix_counter_Z = 0;

  for(int i=0; i<MAXZ; i++)
    {
      if(config[index1][index2][i]==1)
	{
	  pix_position_Z[pix_counter_Z] = i;
	  pix_counter_Z++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter_Z; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position_Z[i]-pix_position_Z[j]);

	if(temp_dist>=MAXZ/2) temp_dist = MAXZ-temp_dist;

	heightS2[index1][index2][temp_dist]++;


      }

  //temp_npH += heightS2[index1][index2][0];

}


void sample_horizontal(int lindex, int cindex)
{

  for(int r=0; r<NtY; r++)
    N2H[lindex][cindex][r] = 0;


  int ener[MAXY];
  int flag_empty = 0;

  for(int i=0; i<MAXY; i++)
    {
      if(config[lindex][i][cindex]==0) ener[i]=-1;
      else
	{
	  int en = 0;
	  int neb1 = i - 1;
	  if(neb1<0) neb1 = neb1 + MAXY;
	  if(config[lindex][neb1][cindex]==1) en++;
	  int neb2 = i+1;
	  if(neb2>=MAXY) neb2 = neb2 - MAXY;
	  if(config[lindex][neb2][cindex]==1) en++;

	  ener[i] = en;
	}
    }

  int position[MAXY];
  for(int i=0; i<MAXY; i++)
    {
      position[i] = -1;
    }

  int ctp = 0;
  for(int i=0; i<MAXY; i++)
    {
      if(ener[i]==1)
	{
	  position[ctp] = i;
	  ctp++;
	}
      else if(ener[i]==0)
	{
	  N2H[lindex][cindex][0]++;
	}
    }

  if(config[lindex][0][cindex]==1&&config[lindex][MAXY-1][cindex]==1)//when the cord is at the bd
    {
      if(ctp>2)
	{
	  for(int i=1; i<ctp-1; i=i+2)
	    {
	      int len = position[i+1]-position[i]+1;
	      for(int r=0; r<=len; r++)
		{
		  if(r<NtY) N2H[lindex][cindex][r] = N2H[lindex][cindex][r]+(len-r);
		}
	    }

	  int len = (position[0]+1)+(MAXY - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtY) N2H[lindex][cindex][r] = N2H[lindex][cindex][r]+(len-r);
	    }
	}
      else if(ctp == 2)
	{
	  int len = (position[0]+1) + (MAXY - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtY) N2H[lindex][cindex][r] = N2H[lindex][cindex][r]+(len-r);
	    }
	}
      else if(ctp == 0 && flag_empty !=0) //this the rare event that the entrie row is full!!!
	{
	  int len = MAXY;
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtY) N2H[lindex][cindex][r] = N2H[lindex][cindex][r]+(len-r);
	    }
	}

    }
  else
    {
      for(int i=0; i<ctp; i=i+2)
	{
	  int len = position[i+1]-position[i]+1;
	  for(int r= 0; r<=len; r++)
	    {
	      if(r<NtY) N2H[lindex][cindex][r] = N2H[lindex][cindex][r]+(len-r);
	    }
	}
    }


  //cout<<"N2H["<<lindex<<" "<<cindex<<"] = "<<N2H[lindex][cindex][0]<<endl;

}

void sample_vertical(int lindex, int cindex)
{
  for(int r=0; r<NtX; r++)
    N2V[lindex][cindex][r] = 0;

  int flag_empty = 0; //gurantee this is not a null row

  int ener[MAXX];
  for(int i=0; i<MAXX; i++)
    {
      if(config[i][lindex][cindex]==0) ener[i]=-1;
      else
	{
	  int en = 0;
	  int neb1 = i - 1;
	  if(neb1<0) neb1 = neb1 + MAXX;
	  if(config[neb1][lindex][cindex]==1) en++;
	  int neb2 = i+1;
	  if(neb2>=MAXX) neb2 = neb2 - MAXX;
	  if(config[neb2][lindex][cindex]==1) en++;

	  ener[i] = en;

	  flag_empty++;
	}
    }

  int position[MAXX];
  for(int i=0; i<MAXX; i++)
    {
      position[i] = -1;
    }

  int ctp = 0;
  for(int i=0; i<MAXX; i++)
    {
      if(ener[i]==1)
	{
	  position[ctp] = i;
	  ctp++;
	}
      else if(ener[i]==0)
	{
	  N2V[lindex][cindex][0]++;
	}
    }

  if(config[0][lindex][cindex]==1&&config[MAXX-1][lindex][cindex]==1)//when the cord is at the bd
    {
      if(ctp>2)
	{
	  for(int i=1; i<ctp-1; i=i+2)
	    {
	      int len = position[i+1]-position[i]+1;
	      for(int r=0; r<=len; r++)
		{
		  if(r<NtX) N2V[lindex][cindex][r] = N2V[lindex][cindex][r]+(len-r);
		}
	    }

	  int len = (position[0]+1) + (MAXX - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtX) N2V[lindex][cindex][r] = N2V[lindex][cindex][r]+(len-r);
	    }
	}
      else if(ctp == 2)
	{

	  //cout<<"ctp = 2"<<endl;
	  int len = (position[0]+1) + (MAXX - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtX) N2V[lindex][cindex][r] = N2V[lindex][cindex][r]+(len-r);
	    }
	}
      else if(ctp == 0 && flag_empty !=0) //this the rare event that the entrie row is full!!!
	{
	  int len = MAXY;
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtX) N2V[lindex][cindex][r] = N2V[lindex][cindex][r]+(len-r);
	    }
	}
    }
  else
    {
      for(int i=0; i<ctp; i=i+2)
	{
	  int len = position[i+1]-position[i]+1;
	  for(int r= 0; r<=len; r++)
	    {
	      if(r<NtX) N2V[lindex][cindex][r] = N2V[lindex][cindex][r]+(len-r);
	    }
	}
    }

  //cout<<"N2V["<<lindex<<" "<<cindex<<"] = "<<N2V[lindex][cindex][0]<<endl;
  //cout<<N2V[lindex][cindex][0]<<endl;

}


void sample_perpendicular(int lindex, int cindex)
{
  for(int r=0; r<NtZ; r++)
    N2D[lindex][cindex][r] = 0;

  int flag_empty = 0;
  int ener[MAXZ];
  for(int i=0; i<MAXZ; i++)
    {
      if(config[lindex][cindex][i]==0) ener[i]=-1;
      else
	{
	  int en = 0;
	  int neb1 = i - 1;
	  if(neb1<0) neb1 = neb1 + MAXZ;
	  if(config[lindex][cindex][neb1]==1) en++;
	  int neb2 = i+1;
	  if(neb2>=MAXZ) neb2 = neb2 - MAXZ;
	  if(config[lindex][cindex][neb2]==1) en++;

	  ener[i] = en;
	}
    }

  int position[MAXZ];
  for(int i=0; i<MAXZ; i++)
    {
      position[i] = -1;
    }

  int ctp = 0;
  for(int i=0; i<MAXZ; i++)
    {
      if(ener[i]==1)
	{
	  position[ctp] = i;
	  ctp++;
	}
      else if(ener[i]==0)
	{
	  N2D[lindex][cindex][0]++;
	}
    }

  if(config[lindex][cindex][0]==1 && config[lindex][cindex][MAXZ-1]==1)//when the cord is at the bd
    {
      if(ctp>2)
	{
	  for(int i=1; i<ctp-1; i=i+2)
	    {
	      int len = position[i+1]-position[i]+1;
	      for(int r=0; r<=len; r++)
		{
		  if(r<NtZ) N2D[lindex][cindex][r] = N2D[lindex][cindex][r]+(len-r);
		}
	    }

	  int len = (position[0]+1) + (MAXZ - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtZ) N2D[lindex][cindex][r] = N2D[lindex][cindex][r]+(len-r);
	    }
	}

      else if(ctp == 2)
	{
	  int len = (position[0]+1) + (MAXX - position[ctp-1]);
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtZ) N2D[lindex][cindex][r] = N2D[lindex][cindex][r]+(len-r);
	    }
	}
      else if(ctp == 0 && flag_empty !=0) //this the rare event that the entrie row is full!!!
	{
	  int len = MAXZ;
	  for(int r=0; r<=len; r++)
	    {
	      if(r<NtZ) N2D[lindex][cindex][r] = N2D[lindex][cindex][r]+(len-r);
	    }
	}

    }
  else
    {
      for(int i=0; i<ctp; i=i+2)
	{
	  int len = position[i+1]-position[i]+1;
	  for(int r= 0; r<=len; r++)
	    {
	      if(r<NtZ) N2D[lindex][cindex][r] = N2D[lindex][cindex][r]+(len-r);
	    }
	}
    }



}



void change_config()
{
  int i, j, k, m, n, l;
  int lim = 0;


  do{   i = rand() % MAXX;
        m = rand() % MAXX;

        j = rand() % MAXY;
        n = rand() % MAXY;

        k = rand() % MAXZ;
        l = rand() % MAXZ;

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

void change_config_fake()
{
  int i, j, k, m, n, l;
  int lim = 0;


  do{   i = rand() % MAXX;
        m = rand() % MAXX;

        j = rand() % MAXY;
        n = rand() % MAXY;

        k = rand() % MAXZ;
        l = rand() % MAXZ;

        lim++;

     } while(config[i][j][k] == config[m][n][l] && lim < 1000);


  indexi = i;
  indexj = j;
  indexk = k;
  indexm = m;
  indexn = n;
  indexl = l;

}


void resume_config()
{
  int temp;
  //first we resume the config
  temp = config[indexi][indexj][indexk];
  config[indexi][indexj][indexk] = config[indexm][indexn][indexl];
  config[indexm][indexn][indexl] = temp;

}


double energy(double SX[NtX],double SY[NtY],double SZ[NtZ], double LX[NtX], double LY[NtY], double LZ[NtZ] )
{
  double E=0;

  for(int i=1; i<NtX; i++)
    {
      E = E + (SX[i] - objS2X[i])*(SX[i] - objS2X[i])+ (LX[i]-objLX[i])*(LX[i]-objLX[i]);
    }
  for(int i=1; i<NtY; i++)
    {
      E = E + (SY[i] - objS2Y[i])*(SY[i] - objS2Y[i])+ (LY[i]-objLY[i])*(LY[i]-objLY[i]);
    }
  for(int i=1; i<NtZ; i++)
    {
      E = E + (SZ[i] - objS2Z[i])*(SZ[i] - objS2Z[i])+ (LZ[i]-objLZ[i])*(LZ[i]-objLZ[i]);
    }

  return E;
}


double d_energy(double S2X[NtX],double S2Y[NtY],double S2Z[NtZ], double ST2X[NtX], double ST2Y[NtY],
                double ST2Z[NtZ],double LX[NtX], double LY[NtY], double LZ[NtZ], double LTX[NtX], double LTY[NtY],  double LTZ[NtZ])
//double d_energy(double SBX[NtX], double SBY[NtY], double SBZ[NtZ], double TSBX[NtX], double TSBY[NtY], double TSBZ[NtZ])
{
  double d_E = 0;
  d_E = energy(ST2X,ST2Y,ST2Z,LTX,LTY,LTZ) - energy(S2X,S2Y,S2Z,LX,LY,LZ);
  //d_E = energy(TSBX, TSBY, TSBZ) - energy(SBX, SBY, SBZ);
  return d_E;
}

double PE(double dE, double T) // the probability that should compare ...
{
  if(dE > 0) return exp(-dE/T);
  else return 1;
}


main()
{
  FILE* fp;
  //double S1 = 0;

  double S2X[NtX]; //for S2...
  double S2Y[NtY];
  double S2Z[NtZ];

  double ST2X[NtX];
  double ST2Y[NtY];
  double ST2Z[NtZ];

  int  SS2X[NtX];
  int  SS2Y[NtY];
  int  SS2Z[NtZ];

  double LX[NtX]; // for L...
  double LY[NtY];
  double LZ[NtZ];

  double LTX[NtX];
  double LTY[NtY];
  double LTZ[NtZ];

  int LSX[NtX];
  int LSY[NtY];
  int LSZ[NtZ];

  int SLT[2][NtX];
  int SCT[2][NtY];
  int SHT[2][NtZ];//save the current values for changed lines, columes and heights...

  int LHT[2][NtX];
  int LVT[2][NtY];
  int LDT[2][NtZ];


  double energyb = MAXS;
  double energyt = 0;


  for(int i=0; i<NtX; i++)
    {
      S2X[i] = 0;
      ST2X[i] = 0;
      SS2X[i] = 0;

      LX[i] = 0;
      LTX[i] = 0;
      LSX[i] =0;
    }

 for(int i=0; i<NtY; i++)
    {
      S2Y[i] = 0;
      ST2Y[i] = 0;
      SS2Y[i] = 0;

      LY[i] = 0;
      LTY[i] = 0;
      LSY[i] =0;
    }

 for(int i=0; i<NtZ; i++)
    {
      S2Z[i] = 0;
      ST2Z[i] = 0;
      SS2Z[i] = 0;

      LZ[i] = 0;
      LTZ[i] = 0;
      LSZ[i] =0;
    }


  for(int i=0; i<MAXY; i++)
    for(int j=0; j<MAXZ; j++)
      for(int k=0; k<NtX; k++)
      {
	lineS2[i][j][k] = 0;
	N2H[i][j][k] = 0;
      }

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXZ; j++)
      for(int k=0; k<NtY; k++)
      {
	columeS2[i][j][k] = 0;
	N2V[i][j][k] = 0;
      }

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<NtZ; k++)
      {
	heightS2[i][j][k] = 0;
	N2D[i][j][k] =0;
      }

  get_obj();

  //read_config();

  read_parameter();

  if(flag_iconfig == 0)
    init_configII();//initialize configuration. volume fraction preserved...
  else
    read_config();

  //now we sample S2 and L for the first time...


  cout<<"initial sampling S2 and L...."<<endl;
  for(int i=0; i<MAXY; i++)
    for(int j=0; j<MAXZ; j++)
      {
	     sampleS2line(i,j);
	     sample_horizontal(i,j);
      }
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXZ; j++)
      {
	    sampleS2colume(i,j);
        sample_vertical(i,j);
      }
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      {
	   sampleS2height(i,j);
       sample_perpendicular(i, j);
      }
  for(int r=0; r<NtX; r++)
    {
      for(int i=0; i<MAXY; i++)
	    for(int j=0; j<MAXZ; j++)
	  {
	    //SS2[r] += (lineS2[i][j][r] + columeS2[i][j][r] + heightS2[i][j][r]);
	    SS2X[r] += lineS2[i][j][r];
	    LSX[r] += N2H[i][j][r];
      }

        S2X[r] = (double) SS2X[r]/(double)(MAXX*MAXY*MAXZ);
        LX[r] = (double) LSX[r]/(double)(MAXX*MAXY*MAXZ);

  //    printf("%d \t  %f \n", r, S2[r]);

    }
   for(int r=0; r<NtY; r++)
    {
      for(int i=0; i<MAXX; i++)
	    for(int j=0; j<MAXZ; j++)
	  {
	    //SS2[r] += (lineS2[i][j][r] + columeS2[i][j][r] + heightS2[i][j][r]);
	    SS2Y[r] += columeS2[i][j][r];
	    LSY[r] += N2V[i][j][r];
      }

        S2Y[r] = (double) SS2Y[r]/(double)(MAXX*MAXY*MAXZ);
        LY[r] = (double) LSY[r]/(double)(MAXX*MAXY*MAXZ);

    //  printf("%d \t  %f \n", r, S2[r]);

    }
      for(int r=0; r<NtZ; r++)
    {
      for(int i=0; i<MAXX; i++)
	    for(int j=0; j<MAXY; j++)
	  {
	    //SS2[r] += (lineS2[i][j][r] + columeS2[i][j][r] + heightS2[i][j][r]);
	    SS2Z[r] += heightS2[i][j][r];
        LSZ[r] += N2D[i][j][r];
      }
        S2Z[r] = (double) SS2Z[r]/(double)(MAXX*MAXY*MAXZ);
        LZ[r] = (double) LSZ[r]/(double)(MAXX*MAXY*MAXZ);
  //    printf("%d \t  %f \n", r, S2[r]);

    }
  cout<<"*****************************************"<<endl;

  fp = fopen("TS2.txt", "w");
  for(int r=0; r<NtX; r++)
   {
      fprintf(fp, "%d \t  %f \n", r, S2X[r]);
    }
    fprintf(fp, "#######################\n");
  for(int r=0; r<NtY; r++)
   {
      fprintf(fp, "%d \t  %f \n", r, S2Y[r]);
    }
    fprintf(fp, "#######################\n");
  for(int r=0; r<NtZ; r++)
   {
      fprintf(fp, "%d \t  %f \n", r, S2Z[r]);
    }
    fprintf(fp, "#######################\n");
  fclose(fp);

  fp = fopen("TL.txt", "w");
  for(int r=0; r<NtX; r++)
    {
      fprintf(fp, "%d \t  %f \n", r, LX[r]);
    }
  for(int r=0; r<NtY; r++)
    {
      fprintf(fp, "%d \t  %f \n", r, LY[r]);
    }
  for(int r=0; r<NtZ; r++)
    {
      fprintf(fp, "%d \t  %f \n", r, LZ[r]);
    }
  fclose(fp);



  //simulated annealing procedure to evlove the system
  //*****************************************************************
  //*****************************************************************

  int Nacc; //for computing acceptance rate...

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //these are for efficiently updating the functions
  int TEMP_S2X[NtX];
  int TEMP_S2Y[NtY];
  int TEMP_S2Z[NtZ];

  int LCHV_S2X[NtX];
  int LCHV_S2Y[NtY];
  int LCHV_S2Z[NtZ];

  int TEMP_LX[NtX];
  int TEMP_LY[NtY];
  int TEMP_LZ[NtZ];

  int LCHV_LX[NtX];
  int LCHV_LY[NtY];
  int LCHV_LZ[NtZ];


  cout<<"Staring the simulated annealing reconstruction process..."<<endl;
  for(int q=0; q<TN; q++)
    {
      T = alpha*T;

      Nacc = 0;

      cout<<"Stage "<<q+1<<" with T = "<<T<<endl;

      for(int i=0; i< Nevl; i++)
	{
	  change_config();
	  //sample S2 for the new configuration, using time saving methods


	  //first save the values of lines and columes that will be changed


	    for(int r=0; r<NtX; r++)
	    {
	      SLT[0][r] = lineS2[indexj][indexk][r];
	      SLT[1][r] = lineS2[indexn][indexl][r];
          LHT[0][r] = N2H[indexj][indexk][r];
	      LHT[1][r] = N2H[indexn][indexl][r];
	    }
        for(int r=0; r<NtY; r++)
        {
	      SCT[0][r] = columeS2[indexi][indexk][r];
	      SCT[1][r] = columeS2[indexm][indexl][r];
	      LVT[0][r] = N2V[indexi][indexk][r];
	      LVT[1][r] = N2V[indexm][indexl][r];
	    }
        for(int r=0; r<NtZ; r++)
        {
          SHT[0][r] = heightS2[indexi][indexj][r];
	      SHT[1][r] = heightS2[indexm][indexn][r];
          LDT[0][r] = N2D[indexi][indexj][r];
	      LDT[1][r] = N2D[indexm][indexn][r];
        }

	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  sampleS2line(indexj, indexk);
	  sampleS2line(indexn, indexl);
	  //only sample the lines that changed...new results are stored in lineS2[lindexi,lindexm][*]

	  sampleS2colume(indexi, indexk);
	  sampleS2colume(indexm, indexl);
	  //only sample the columes that changed...new results are stored in columeS2[cindexi,cindexm][*]

	  sampleS2height(indexi, indexj);
	  sampleS2height(indexm, indexn);

	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  sample_horizontal(indexj, indexk);
	  sample_horizontal(indexn, indexl);

	  sample_vertical(indexi, indexk);
	  sample_vertical(indexm, indexl);

	  sample_perpendicular(indexi, indexj);
	  sample_perpendicular(indexm, indexn);


	  //Now we compute the S2 and L for the new configuration...
	  for(int r=0; r<NtX; r++)
	    {
	      //the following method only consider the changes..
	      //************************************************
	      TEMP_S2X[r] = 0; //old value
	      LCHV_S2X[r] = 0; //current value

	      TEMP_LX[r] = 0; //old value
	      LCHV_LX[r] = 0; //current value


	      for(int vv=0; vv<2; vv++)
		{
		  TEMP_S2X[r] += SLT[vv][r];
		  TEMP_LX[r]  += LHT[vv][r];
		}

	       LCHV_S2X[r] = lineS2[indexj][indexk][r]+lineS2[indexn][indexl][r];
	       LCHV_LX[r] = N2H[indexj][indexk][r]+N2H[indexn][indexl][r];


	      if(r!=0)
		{
		  //if(TEMP_S2 != LCHV_S2) cout<<"TEMP_S2, LCHV_S2 = "<<TEMP_S2<<" "<<LCHV_S2<<endl;
		  //if(TEMP_L[r] != LCHV_L[r]) cout<<"TEMP_L, LCHV_L = "<<TEMP_L[r]<<" "<<LCHV_L[r]<<endl;
		}


	      SS2X[r] =(int)( SS2X[r] - TEMP_S2X[r] + LCHV_S2X[r]);

	      LSX[r] = (int)(LSX[r] - TEMP_LX[r] + LCHV_LX[r]);

	      ST2X[r] = (double)SS2X[r]/(double)(MAXX*MAXY*MAXZ);

	      LTX[r] = (double)LSX[r]/(double)(MAXX*MAXY*MAXZ);

	      //************************************************
	    }
        for(int r=0; r<NtY; r++)
	    {
	      //the following method only consider the changes..
	      //************************************************
	      TEMP_S2Y[r] = 0; //old value
	      LCHV_S2Y[r] = 0; //current value

	      TEMP_LY[r] = 0; //old value
	      LCHV_LY[r] = 0; //current value


	      for(int vv=0; vv<2; vv++)
		{
		  TEMP_S2Y[r] += SCT[vv][r];
		  TEMP_LY[r] +=  LVT[vv][r];
		}

           LCHV_S2Y[r] = columeS2[indexi][indexk][r]+columeS2[indexm][indexl][r];
           LCHV_LY[r] = N2V[indexi][indexk][r]+N2V[indexm][indexl][r];

	      if(r!=0)
		{
		  //if(TEMP_S2 != LCHV_S2) cout<<"TEMP_S2, LCHV_S2 = "<<TEMP_S2<<" "<<LCHV_S2<<endl;
		  //if(TEMP_L[r] != LCHV_L[r]) cout<<"TEMP_L, LCHV_L = "<<TEMP_L[r]<<" "<<LCHV_L[r]<<endl;
		}



	      SS2Y[r] =(int)( SS2Y[r] - TEMP_S2Y[r] + LCHV_S2Y[r]);

	      LSY[r] = (int)(LSY[r] - TEMP_LY[r] + LCHV_LY[r]);

	      ST2Y[r] = (double)SS2Y[r]/(double)(MAXX*MAXY*MAXZ);

	      LTY[r] = (double)LSY[r]/(double)(MAXX*MAXY*MAXZ);


	      //************************************************
	    }

	    for(int r=0; r<NtZ; r++)
	    {
	      //the following method only consider the changes..
	      //************************************************
	      TEMP_S2Z[r] = 0; //old value
	      LCHV_S2Z[r] = 0; //current value

	      TEMP_LZ[r] = 0; //old value
	      LCHV_LZ[r] = 0; //current value


	      for(int vv=0; vv<2; vv++)
		{
		  TEMP_S2Z[r] += SHT[vv][r];
		  TEMP_LZ[r] +=  LDT[vv][r];
		}

           LCHV_S2Z[r] =heightS2[indexi][indexj][r]+heightS2[indexm][indexn][r];
	       LCHV_LZ[r] = N2D[indexi][indexj][r]+N2D[indexm][indexn][r];


	      if(r!=0)
		{
		  //if(TEMP_S2 != LCHV_S2) cout<<"TEMP_S2, LCHV_S2 = "<<TEMP_S2<<" "<<LCHV_S2<<endl;
		  //if(TEMP_L[r] != LCHV_L[r]) cout<<"TEMP_L, LCHV_L = "<<TEMP_L[r]<<" "<<LCHV_L[r]<<endl;
		}



	      SS2Z[r] =(int)( SS2Z[r] - TEMP_S2Z[r] + LCHV_S2Z[r]);

	      LSZ[r] = (int)(LSZ[r] - TEMP_LZ[r] + LCHV_LZ[r]);

	      ST2Z[r] = (double)SS2Z[r]/(double)(MAXX*MAXY*MAXZ);

	      LTZ[r] = (double)LSZ[r]/(double)(MAXX*MAXY*MAXZ);

	      //************************************************
	    }
	  //Monte Carlo steps...

	  double P = double (rand() % MAXS)/(double) MAXS;

	  if( P > PE(d_energy(S2X,S2Y,S2Z,ST2X,ST2Y,ST2Z,LX,LY,LZ,LTX,LTY,LTZ), T))
      	    {
	      resume_config();
	      //this just resumes the 'configuration', still need to resume S2...



	      for(int r=0; r<NtX; r++)
		{

		  SS2X[r] = (int)(SS2X[r] + TEMP_S2X[r] - LCHV_S2X[r]);
		  LSX[r] = (int)(LSX[r] + TEMP_LX[r] - LCHV_LX[r]);

		  //S2 does not change, do not need to resume that...
		  //now we resume lineS2, columeS2 and heightS2;
		  lineS2[indexj][indexk][r] = SLT[0][r];
		  lineS2[indexn][indexl][r] = SLT[1][r];

		  N2H[indexj][indexk][r] = LHT[0][r];
		  N2H[indexn][indexl][r] = LHT[1][r];
		}
	      for(int r=0; r<NtY; r++)
		{

		  SS2Y[r] = (int)(SS2Y[r] + TEMP_S2Y[r] - LCHV_S2Y[r]);
		  LSY[r] = (int)(LSY[r] + TEMP_LY[r] - LCHV_LY[r]);

		  columeS2[indexi][indexk][r] = SCT[0][r];
		  columeS2[indexm][indexl][r] = SCT[1][r];

		  N2V[indexi][indexk][r] = LVT[0][r];
		  N2V[indexm][indexl][r] = LVT[1][r];
		}
         for(int r=0; r<NtZ; r++)
		{

		  SS2Z[r] = (int)(SS2Z[r] + TEMP_S2Z[r] - LCHV_S2Z[r]);
		  LSZ[r] = (int)(LSZ[r] + TEMP_LZ[r] - LCHV_LZ[r]);

		  heightS2[indexi][indexj][r] = SHT[0][r];
		  heightS2[indexm][indexn][r] = SHT[1][r];

		  N2D[indexi][indexj][r] = LDT[0][r];
		  N2D[indexm][indexn][r] = LDT[1][r];
		}
	        }

         else
	    {
	      for(int r=0; r<NtX; r++)
		{
		  S2X[r] = ST2X[r];
		  LX[r] = LTX[r];
		}
		   for(int r=0; r<NtY; r++)
		{
		  S2Y[r] = ST2Y[r];
		  LY[r] = LTY[r];
		}
		   for(int r=0; r<NtZ; r++)
		{
		  S2Z[r] = ST2Z[r];
		  LZ[r] = LTZ[r];
		}

	      Nacc++;
	    }


	  //compare and record the best energy and configuration...
	  energyt = energy(S2X,S2Y,S2Z,LX,LY,LZ);

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
      cout<<"The acceptance rate: "<<(double)Nacc/(double)Nevl<<endl;
      cout<<"The energy E = "<<energyb<<endl;


      fp = fopen("E.txt","a");
      fprintf(fp, "%e\n", energyb);
      fclose(fp);


      printf("*************************************************\n");

    }

  //*****************************************************************
  //*****************************************************************
  //this is the end of simulated annealing

  fp = fopen("Fconfig.txt","w");
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
	{
	  if(config[i][j][k] == 1)
	    {
	      fprintf(fp, "%d \t %d \t %d \n", i, j, k );
	    }
	}
  fclose(fp);

   fp = fopen("S2.txt", "w");
   fprintf(fp, "%d \t %f \n", 0, f1);
   for(int r=1; r<NtX; r++)
     {
       fprintf(fp, "%d \t %f \n", r, S2X[r]);
     }
   fprintf(fp, "####################\n\n");
   fprintf(fp, "%d \t %f \n", 0, f1);
   for(int r=1; r<NtY; r++)
     {
       fprintf(fp, "%d \t %f \n", r, S2Y[r]);
     }
   fprintf(fp, "####################\n\n");
   fprintf(fp, "%d \t %f \n", 0, f1);
   for(int r=1; r<NtZ; r++)
     {
       fprintf(fp, "%d \t %f \n", r, S2Z[r]);
     }
   fprintf(fp, "####################\n\n");
   fclose(fp);

   fp = fopen("L.txt", "w");
   fprintf(fp, "%d \t %f \n", 0, f1);
   for(int r=1; r<NtX; r++)
     {
       fprintf(fp, "%d \t %f \n", r, LX[r]);
     }
   fprintf(fp, "####################\n\n");
   fprintf(fp, "%d \t %f \n", 0, f1);
   for(int r=1; r<NtY; r++)
     {
       fprintf(fp, "%d \t %f \n", r, LY[r]);
     }
   fprintf(fp, "####################\n\n");
   fprintf(fp, "%d \t %f \n", 0, f1);
   for(int r=1; r<NtZ; r++)
     {
       fprintf(fp, "%d \t %f \n", r, LZ[r]);
     }
   fprintf(fp, "####################\n\n");
   fclose(fp);

  //this is the end of the codes...

}
