#ifndef RAND_H
#define RAND_H

#include "..\global.h"

class SimpleRand
{
public:
// generate random number in [0,1)

// generate random number between [imin,imax)
static int rU_Int(int imin,int imax)
{
  return (rand()%(imax-imin))+imin;
}

static double rU_RightOpen()
{
  return ((double)rand())/(double(RAND_MAX)+1.0);
}

// generate random number in (0,1)
static double rU_Open()
{
  return ((double)rand()+0.5)/(double(RAND_MAX)+1.0);
}

// generate random number in [0,1]
static double rU_Closed()
{
  return ((double)rand())/(double(RAND_MAX));
}

// generate random number in [imin,imax]
static double rU_Intervall(double imin,double imax)
{
  return rU_Closed()*(imax-imin)+imin;
}


static double rExp(double lambda)
{
  return -log(rU_Open())/lambda;
}


// bernoulli-trial
static bool BoolEvent(double p)
{
  return (rU_RightOpen() < p);
}

static int Geometric(double p)
{
  return (ceil(log(1-SimpleRand::rU_Open())/log(1-p)));
}

};

#endif
