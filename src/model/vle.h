#ifndef VLE_H
#define VLE_H

#include "../global.h"

// calculate optimal rice param k
class RiceParam {
  public:
   RiceParam(){
     cout << "Calculating optimal Rice-Param..." << endl;
     Calc();
   };
    int GetOptimalK(int val) {
      for (int k=0;k<kmax;k++) {
        if (val<=ub[k]) return k;
      }
      cout << " warning (Rice): optimal k for val: " << val << " not found!" << endl;
      return 0;
    }
    int CalcRiceOptimalK(int val)
    {
      for (int k=0;k<28;k++)
      {
         if (RiceCost(val,k)<RiceCost(val,k+1)) return k;
      }
      cout << "rp: could not find optimal k\n";
      return -1;
    }
    int CalcExpOptimalK(int val)
    {
      for (int k=0;k<28;k++)
      {
         if (ExpCost(val,k)<ExpCost(val,k+1)) return k;
      }
      cout << "rp: could not find optimal k\n";
      return -1;
    }
  private:
    void Calc() {
      kmax=0;
      vmax=1<<28;
      for (int val=0;val<vmax;val++) {
         if (RiceCost(val,kmax+1)<=RiceCost(val,kmax)) {
            cout << "upper bound for k=" << kmax << ": " << (val-1) << endl;
            ub[kmax]=val-1;
            kmax++;
         }
      }
    }
    int RiceCost(int val,int k)
    {
      int q=val>>k;
      return q+1+k;
    }
    int ExpCost(int val,int k) {
        int c=0;
        while (val >= (1<<k)) {
            c++;
            val-=(1<<k);
            k++;
        }
        c+=1;
        c+=k;
        return c;
    }
    int ub[32]; // upper bound for optimal k
    int vmax,kmax;
};


#endif // VLE_H
