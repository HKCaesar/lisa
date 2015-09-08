#ifndef MIXER_H
#define MIXER_H

#include "model.h"
#include "..\global.h"

class Mix2Linear
{
  public:
    Mix2Linear(){Init(WSCALEh);};
    void Init(int iw){w=iw;};
    //pm=(1-w)*p1+w*p2
    int Predict(int _p1,int _p2)
    {
      p1=_p1;p2=_p2;
      pm = p1+idiv_signed32((p2-p1)*w,WBITS);
      pm = clamp(pm,1,PSCALEm);
      return pm;
    }
    void Update(int bit,int rate)
    {
      int e=(bit<<PBITS)-pm;
      int d=idiv_signed32((p2-p1)*e,PBITS);
      upd_w(d,rate);
    }
    int w,p1,p2,pm;
  protected:
    inline int idiv_signed32(int val,int s){return val<0?-(((-val)+(1<<(s-1)))>>s):(val+(1<<(s-1)))>>s;};
    inline void upd_w(int d,int rate) {int wd=idiv_signed32(rate*d,PBITS);w=clamp(w+wd,0,int(WSCALE));};
};

#endif // MIXER_H
