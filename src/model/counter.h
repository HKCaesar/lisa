#ifndef COUNTER_H
#define COUNTER_H

#include "model.h"

// stationary 16-bit counter
class StatCounter {
  public:
      StatCounter(int add=1,int limit=255):m_add(add),m_limit(limit){n0=n1=1;};
      int p1(){return ((int(n1)+1)*PSCALE)/(int(n0)+int(n1)+2);};
      void update(int bit)
      {
        bit?n1+=m_add:n0+=m_add;
        if (n1>m_limit || n0>m_limit) {n1=(n1+1)>>1;n0=(n0+1)>>1;};
      };
  protected:
      int m_add,m_limit;
      int n0,n1;
};
#endif // COUNTER_H
