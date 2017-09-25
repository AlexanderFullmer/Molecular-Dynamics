#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define VAdd(v1, v2, v3){\
    (v1).x = (v2).x + (v3).x;\
    (v1).y = (v2).y + (v3).y;\
}                                       //Vector Addition
#define VVAdd(v1, v2)   VAdd(v1, v1, v2)  //v1 = v1 + v2
#define VSub(v1, v2, v3){\
    (v1).x = (v2).x - (v3).x;\
    (v1).y = (v2).y - (v3).y;\
}                                       //Vector Subtraction
#define VMul(v1, v2, v3){\
  (v1).x = (v2).x * (v3).x;\
  (v1).y = (v2).y * (v3).y;\
}                                       //Multiplication
#define VDiv(v1, v2, v3){\
  (v1).x = (v2).x / (v3).x;\
  (v1).y = (v2).y / (v3).y;\
}                                       //Division
#define VDot(v1, v2)\
    (v1).x * (v2).x + (v1).y * (v2).y   //Dot Product
#define VLenSq(v)        VDot(v, v)     //Vector Squared
#define VLen(v) sqrt(VDot (v, v))       //Vector length
#define VSAdd(v1, v2, s3, v3){\
  (v1).x = (v2).x + (s3) * (v3).x;\
  (v1).y = (v2).y + (s3) * (v3).y;\
}                                         //Vector Addition with Constant
#define VVSAdd(v1, s2, v2)    VSAdd(v1, v1, s2, v2) //(v1) = (v1) + (s2) * (v2)
#define VCSum(v)\
    ((v).x + (v).y)
#define VScale(v, s){\
  (v).x *= s;\
  (v).y *= s;\
}                                       //Rescale a Vector
#define VSet(v, sx, sy){\
  (v).x = sx;\
  (v).y = sy;\
}                                       //Set Vector Components
#define VSetAll(v, s)         VSet(v,s,s)     //Set All Components to Const.
#define VZero(v)              VSetAll(v, 0)  //Set All Components to 0
#define VSCopy(v2, s1, v1){\
  (v2).x = (s1) * (v1).x;\
  (v2).y = (s1) * (v1).y;\
}                                 //Multiplication of Vector with Const.
#define VProd(v)  ((v).x * (v).y)       //Multiplication of x and y Components
#define VWrap(v, t){\
  if (v.t >= 0.5 * region.t) v.t -= region.t;\
  else if (v.t < -0.5 * region.t) v.t += region.t;\
}                                           //Periodic Wraparound
#define VWrapAll(v){\
    VWrap(v, x);\
    VWrap(v, y);\
}
#define Sqr(x)              ((x) * (x))
#define Cube(x)             ((x) * (x) * (x))
#define DO_MOL              for (n = 0; n < nMol; n++)
#define NDIM  2
#define AllocMem(a, n, t) a = (t *) malloc ((n) * sizeof (t))
#define PropZero(v){\
  v.sum = 0.;\
  v.sum2 = 0.;\
}
#define PropAccum(v){\
  v.sum += v.val;\
  v.sum2 += Sqr (v.val);\
}
#define Min(x1, x2)\
  ({typeof (x1) _x1 = (x1);\
    typeof(x2) _x2 = (x2);\
    _x1 < _x2 ? _x1 : _x2;})
#define Max(x1,x2)\
({typeof (x1) _x1 = (x1);\
    typeof(x2) _x2 = (x2);\
    _x1 > _x2 ? _x1 : _x2;})
#define PropAvg(v, n){\
  v.sum /= n;\
  v.sum2 = sqrt (Max (v.sum2 /n - Sqr(v.sum), 0.));\
}
#define PropEst(v)\
  v.sum, v.sum2
#define NameI(x) {#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x) {#x, &x, N_R, sizeof (x) / sizeof (real)}
#define NP_I ((int*) (nameList[k].vPtr) + j)
#define NP_R ((real*) (nameList[k].vPtr) + j)
#define IADD 453806245
#define IMUL 314159269
#define MASK 2147483647
#define SCALE 0.4656612873e-9

typedef double real;
typedef enum {N_I, N_R} VType;

typedef struct{
    real x, y;
} VecR;                     //2D Vector

typedef struct{             //Molecule with Position, Velocity, & Accel Vectors
  VecR r, rv, ra;
} Mol;


typedef struct{
  int x, y;
} VecI;


typedef struct{
  real val, sum, sum2;
} Prop;

typedef struct{
    char* vName;
    void* vPtr;
    VType vType;
    int vLen, vStatus;
} NameList;

Mol *mol;
VecR region, vSum;
Prop kinEnergy, pressure, totEnergy;
real rCut, timeNow, uSum, velMag,
     virSum, vvSum, *histVel, rangeVel, hFunction;
int  moreCycles, nMol, stepCount, countVel, limitVel,
     sizeHistVel, stepVel, randSeedP = 17;
FILE* fp;

real deltaT = 0.005;
real density = 0.8;
VecI initUcell = {20, 20};
int  stepAvg = 100;
int  stepEquil = 0;
int  stepLimit = 10000;
real  temperature = 1.;

int GetNameList();
void PrintNameList();
void SetParams ();
void SetupJob();
void AllocArrays ();
void InitCoords ();
real RandR();
void VRand(VecR*);
void InitVels ();
void InitAccels ();
void SingleStep();
void LeapfrogStep (int);
void ApplyBoundaryCond ();
void ComputeForces();
void EvalProps();
void AccumProps (int);
void PrintSummary (FILE*);
void EvalVelDist ();
void PrintVelDist (FILE*);







int main(int argc, char** argv)
{
  SetParams();
  SetupJob();
  moreCycles = 1;
  fp = fopen("MD.txt","w");
  while (moreCycles){
    SingleStep();
    if (stepCount >= stepLimit) moreCycles = 0;
  }
  fclose(fp);
}



void SetParams ()
{
  rCut = pow (2., 1./6.);
  VSCopy (region, 1. / sqrt (density), initUcell);
  nMol = VProd (initUcell);
  velMag = sqrt (NDIM * (1. - 1. / nMol) * temperature);
}


void SetupJob()
{
  AllocArrays();
  stepCount = 0;
  InitCoords ();
  InitVels();
  InitAccels();
  AccumProps(0);
}


void AllocArrays ()
{
  AllocMem (mol, nMol, Mol);
}


void InitCoords ()
{
  VecR c, gap;
  int n, nx, ny;

  VDiv (gap, region, initUcell);
  n = 0;
  for (ny = 0; ny < initUcell.y; ny++){
    for (nx = 0; nx < initUcell.x; nx++){
      VSet (c, nx + 0.5, ny + 0.5);
      VMul (c, c, gap);
      VVSAdd (c, -0.5, region);
      mol[n].r = c;
      ++n;
    }
  }
}

real RandR(){
    randSeedP = (randSeedP*IMUL + IADD) & MASK;
    return (randSeedP * SCALE);
}

void VRand(VecR* p){
    real s;

    s = 2. * M_PI * RandR();
    p -> x = cos(s);
    p -> y = sin(s);
}


void InitVels ()
{
  int n;

  VZero (vSum);
  DO_MOL {
    VRand (&mol[n].rv);
    VScale (mol[n].rv, velMag);
    VVAdd (vSum, mol[n].rv);
  }
  DO_MOL VVSAdd (mol[n].rv, -1. / nMol, vSum);
}


void InitAccels ()
{
  int n;

  DO_MOL VZero (mol[n].ra);
}


void SingleStep()
{
  ++stepCount;
  timeNow = stepCount * deltaT;
  LeapfrogStep(1);
  ApplyBoundaryCond();
  ComputeForces();
  LeapfrogStep(2);
  EvalProps();
  AccumProps(1);
  if(stepCount % stepAvg == 0){
    AccumProps(2);
    PrintSummary(fp);
    AccumProps(0);
  }
}


void LeapfrogStep (int part)
{
  int n;

  if (part == 1){
    DO_MOL{
      VVSAdd (mol[n].rv, 0.5 * deltaT, mol[n].ra);
      VVSAdd (mol[n].r, deltaT, mol[n].rv);
    }
  }
  else{
    DO_MOL VVSAdd (mol[n].rv, 0.5 * deltaT, mol[n].ra);
  }
}


void ApplyBoundaryCond ()
{
  int n;

  DO_MOL VWrapAll (mol[n].r);
}


void ComputeForces()
{
  VecR dr;
  real fcVal, rr, rrCut, rri, rri3;
  int j1, j2, n;

  rrCut = Sqr(rCut);
  DO_MOL VZero (mol[n].ra);
  uSum = 0.;
  virSum = 0.;
  for (j1 = 0; j1 < nMol - 1; j1++){
    for (j2 = j1 + 1; j2 < nMol; j2++){
      VSub(dr, mol[j1].r, mol[j2].r);
      VWrapAll(dr);
      rr = VLenSq(dr);
      if (rr < rrCut){
        rri = 1. / rr;
        rri3 = Cube(rri);
        fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
        VVSAdd(mol[j1].ra, fcVal, dr);
        VVSAdd(mol[j2].ra, -fcVal, dr);
        uSum += 4. * rri3 * (rri3 - 1.) + 1.;
        virSum += fcVal * rr;
      }
    }
  }
}


void EvalProps()
{
  real vv;
  int n;

  VZero (vSum);
  vvSum = 0.;
  DO_MOL{
    VVAdd (vSum, mol[n].rv);
    vv = VLenSq (mol[n].rv);
    vvSum += vv;
  }
  kinEnergy.val = 0.5 * vvSum / nMol;
  totEnergy.val = kinEnergy.val + uSum / nMol;
  pressure.val = density * (vvSum + virSum) / (nMol * NDIM);
}


void AccumProps (int icode)
{
  if (icode == 0){
    PropZero (totEnergy);
    PropZero (kinEnergy);
    PropZero (pressure);
  }
  else if (icode == 1){
    PropAccum (totEnergy);
    PropAccum (kinEnergy);
    PropAccum (pressure);
  }
  else if (icode == 2){
    PropAvg (totEnergy, stepAvg);
    PropAvg (kinEnergy, stepAvg);
    PropAvg (pressure, stepAvg);
  }
}




void PrintSummary (FILE *fp)
{
  fprintf (fp,
      "%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
      stepCount, timeNow, VCSum (vSum) / nMol, PropEst (totEnergy),
      PropEst (kinEnergy), PropEst (pressure));
}
