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
    int ub[32]; // upper bound for optimal k
    int vmax,kmax;
};

extern RiceParam myRiceParam;



#endif // VLE_H
