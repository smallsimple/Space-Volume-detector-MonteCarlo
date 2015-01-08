#include <stdio.h>
#include <string.h>
#include <cmath>
#include <stdlib.h>
#include <sys/times.h>
#include <limits.h>
#define N2 100
#define N3 1000
#define N4 10000
#define N5 100000
#define N6 1000000
#define N7 10000000
#define two31minus1 2147483647

FILE *fp_metaball, *fp_write;
char metaball[30];
char writename[50];
int nMetaballs;
double *rho, *lambda, constantC;
hrtime_t start, end;
float ElapsedTime;

void readfine();
double EvalGaussian(double x, double y, double z);
unsigned int bracket(double x, double y, double z);
void MCArea(long int ,int );


typedef struct {
  double cx, cy, cz;
} Center;

Center *centers;

int main(int argc, char **argv)
{
    if (argc != 3)
    {
      fprintf(stderr, "usage: MonteCarlo <inFile .mball> <outFile .txt> \n");
      return 1;
    }
    strcpy(metaball, argv[1]);
    strcpy(writename,argv[2]);
    readfine();
    //printf("%f\t %f\t %d\t %d\t %f\n", rho[0], lambda[0], nMetaballs, N2, constantC);

    MCArea(N2, 1);
    MCArea(N3, 1);
    MCArea(N4, 1);
    MCArea(N5, 1);
    MCArea(N6, 1);
    MCArea(N7, 1);
	fprintf(fp_write, "\n");
    MCArea(N2, 2);
    MCArea(N3, 2);
    MCArea(N4, 2);
    MCArea(N5, 2);
    MCArea(N6, 2);
    MCArea(N7, 2);
	fprintf(fp_write, "\n");
    MCArea(N2, 3);
    MCArea(N3, 3);
    MCArea(N4, 3);
    MCArea(N5, 3);
    MCArea(N6, 3);
    MCArea(N7, 3);

    //printf("hello\n");
    fclose(fp_metaball);
    fclose(fp_write);
}

void MCArea(long int number, int option)
{
    double xr, yr, zr, volume;
    int nRoots;
    long int in=0, number2;
    start = gethrtime();

	for (long int j = 0; j < number; j++)
    {
        if (option == 3)
        {
            xr = 2*drand48()-1;
            yr = 2*drand48()-1;
            zr = 2*drand48()-1;
        }
        if (option == 2)
        {
            xr = 2*(double)random()/(double)(two31minus1) -1.0;
            yr = 2*(double)random()/(double)(two31minus1) -1.0;
            zr = 2*(double)random()/(double)(two31minus1) -1.0;
        }
        if (option == 1)
        {
            xr = 2*(double)rand()/(double)RAND_MAX -1.0;
            yr = 2*(double)rand()/(double)RAND_MAX -1.0;
            zr = 2*(double)rand()/(double)RAND_MAX -1.0;
			//printf("%f\t%f\t%f\n", xr, yr, zr);
        }
        nRoots = bracket(xr, yr, zr);

        if (nRoots % 2 == 1)
        in++;

		//printf("%d\t%d\n", in, j);
    }

    volume = (double)in * 8 / (double)number;
    end = gethrtime();
    ElapsedTime = (float)(end - start) * 1.0e-09;
    if(option == 1)
    {
        fprintf(fp_write, "rand() mode: N = %10d \t volume = %f    Time = %.9f seconds\n", number, volume, ElapsedTime);
        printf("rand() mode: N = %10d \t volume = %f    Time = %.9f seconds\n", number, volume, ElapsedTime);
    }

    if(option == 2)
    {
        fprintf(fp_write, "random() mode: N = %10d \t volume = %f    Time = %.9f seconds\n", number, volume, ElapsedTime);
        printf("random() mode: N = %10d \t volume = %f    Time = %.9f seconds\n", number, volume, ElapsedTime);
    }
    if(option == 3)
    {
        fprintf(fp_write, "drand48() mode: N = %10d \t volume = %f    Time = %.9f seconds\n", number, volume, ElapsedTime);
        printf("drand48() mode: N = %10d \t volume = %f    Time = %.9f seconds\n", number, volume, ElapsedTime);
    }
}


void readfine()
{
    fp_write = fopen(writename, "w");
      if ((fp_metaball = fopen(metaball, "r")) == 0)
    {
      fprintf(fp_write, "ERROR: Opening file '%s' for read.\n", metaball);
      exit(0);
    }

  fscanf(fp_metaball, "%d\n", &(nMetaballs));

  if (nMetaballs >= 1)
    {

      if ((centers = (Center *)malloc(sizeof(Center) * nMetaballs)) == 0)
	{
	  fprintf(fp_write, "ERROR: Allocate centers array.\n", metaball);
	  exit(0);
	}

      for (unsigned int i = 0; i < nMetaballs; i++)
	{
	  float cx, cy, cz;
	  fscanf(fp_metaball, "%f %f %f\n", &cx, &cy, &cz);

	  centers[i].cx = cx;
	  centers[i].cy = cy;
	  centers[i].cz = cz;
	}

    if ((rho = ( double *)malloc(sizeof(double) * nMetaballs)) == 0)
    {
        fprintf(fp_write, "ERROR: Allocate rho array.\n", metaball);
        exit(0);
    }
    for (unsigned int i = 0; i < nMetaballs; i++)
    {
        float r;
        fscanf (fp_metaball, "%f\n", &r);
        rho[i] = r;
    }

    if ((lambda = (double *)malloc(sizeof(double) * nMetaballs)) == 0)
    {
        fprintf(fp_write, "ERROR: Allocate lambda array. \n", metaball);
        exit(0);
    }
    for (unsigned int i = 0; i< nMetaballs; i++)
    {
        float l;
        fscanf (fp_metaball, "%f\n", &l);
        lambda[i] =l;
    }

    if (1)
    {
        float c;
        fscanf (fp_metaball, "%f\n", &c);
        constantC = c;
    }
    ////
    }
  else
    {
      fprintf(stderr, "WARNING: No metaballs in file.\n");
      nMetaballs = 0;
      centers = 0;
    }
}

double EvalGaussian(double x, double y, double z)
{
    double sum=0;
    for(int m=0; m<nMetaballs; m++)
    {
        sum += lambda[m]*exp(-rho[m]*((x - centers[m].cx)*(x - centers[m].cx) + (y - centers[m].cy)*(y - centers[m].cy)+(z-centers[m].cz)*(z-centers[m].cz)));
    }

    return sum - constantC;
}

unsigned int bracket(double x, double y, double z)
{
    double *zroot;
    int count = 0;
    int n= (int)((z+1)/.1);
    n = n + 1;
    zroot=(double *)malloc(sizeof(double)*(n));


    for (int i= 0; i< n-1; i++)
    {zroot[i]=-1 + i*.1; }

    zroot[n-1] = z;

    for (int i= 0; i< n -1; i++)
    {

        if (EvalGaussian(x, y, zroot[i])*EvalGaussian(x,y,zroot[i+1])<=0)
        {
            count ++;
        }
    }
	free(zroot);
    return count;
}
