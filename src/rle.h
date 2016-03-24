#ifndef RLE_H
#define RLE_H

#include "global.h"


// simple bitpack-compression for clusterfiles
class RLEPack2 {
  public:
      static void EncodeLabel(BitBuffer &bitout,int64_t diff,int nrun)
      {
        bitout.PutEliasGamma(std::abs(diff)+2);
        bitout.PutBit(diff<0);
        bitout.PutEliasGamma(nrun+1);
      }
      static int PackRow(int64_t *rowdata,uint32_t width,uint8_t *dstdata)
      {
        BitBuffer bitout(dstdata);
        int64_t lval=rowdata[0];
        int64_t llabel=0;
        int nrun=0;
        for (uint32_t i=1;i<width;i++) {
            int64_t val=rowdata[i];

            if (val==lval) nrun++;
            else {
              if (lval>0) {EncodeLabel(bitout,lval-llabel,nrun);llabel=lval;}
              else {bitout.PutEliasGamma(1);bitout.PutEliasGamma(nrun+1);}
              lval=val;
              nrun=0;
            }
        }
        if (lval>0) {EncodeLabel(bitout,lval-llabel,nrun);llabel=lval;}
        else {bitout.PutEliasGamma(1);bitout.PutEliasGamma(nrun+1);}
        bitout.Flush();
        return bitout.GetBytesProcessed();
      }
      static void UnpackRow(uint8_t *srcdata,uint32_t width,int64_t *dstdata)
      {
        memset(dstdata,0,width*sizeof(int64_t));
        BitBuffer bitin(srcdata);
        int nrun;
        uint32_t i=0;
        int64_t llabel=0;
        while (i<width) {
          int64_t diff=bitin.GetEliasGamma();
          if (diff>1) {
             diff-=2;
             bool sgn=bitin.GetBit();
             nrun=bitin.GetEliasGamma();
             if (sgn) diff=-diff;
             int64_t label=diff+llabel;
             for (int k=0;k<nrun;k++) dstdata[i+k]=label;
             llabel=label;
          } else nrun=bitin.GetEliasGamma();
          i+=nrun;
        }
      }
};

/*struct trun {
  int64_t label;
  uint32_t pos;
  uint32_t nrun;
};

class RLECount {
  public:
    static void GetRuns(int64_t *rowdata,uint32_t width,vector <trun>&runs,int &labeltype,int &runtype)
    {
        runs.clear();
        labeltype=runtype=0;

        int64_t maxlabel=0;
        uint32_t maxnrun=0;

        uint32_t i=0;
        while (i<width) { // count runs
         if (rowdata[i]) {
          uint32_t nrun=1;
          int64_t label=rowdata[i];
          while ((i+nrun<width && label==rowdata[i+nrun])) nrun++;
          trun myrun;
          myrun.label=label;
          myrun.pos=i;
          myrun.nrun=nrun;
          runs.push_back(myrun);
          if (nrun>maxnrun) maxnrun=nrun;
          if (label>maxlabel) maxlabel=label;
          i+=nrun;
         } else i++;
       }
       if (maxlabel<(1LL<<24)) labeltype=0;
       else if (maxlabel<(1LL<<32)) labeltype=1;
       else if (maxlabel<(1LL<<48)) labeltype=2;
       else labeltype=3;

       if (maxnrun>=(1<<16)) runtype=1;
    }
};


// simple bitpack-compression for clusterfiles
class RLEPack {
  public:
static int PackRow(int64_t *rowdata,uint32_t width,uint8_t *dstdata)
{
  std::vector <trun> runs;
  int labeltype,runtype;
  RLECount::GetRuns(rowdata,width,runs,labeltype,runtype);

  Utils::Put32LH(dstdata,runs.size());
  dstdata[4]=(uint8_t)labeltype;
  dstdata[5]=(uint8_t)runtype;
  int outidx=6;
  for (size_t i=0;i<runs.size();i++)
  {
     int64_t label=runs[i].label;
     uint32_t pos =runs[i].pos;
     uint32_t nrun=runs[i].nrun;
     switch (labeltype) {
       case 0:Utils::Put24LH(dstdata+outidx,label);outidx+=3;break;
       case 1:Utils::Put32LH(dstdata+outidx,label);outidx+=4;break;
       case 2:Utils::Put48LH(dstdata+outidx,label);outidx+=6;break;
       case 3:Utils::Put64LH(dstdata+outidx,label);outidx+=8;break;
       default:cout << "warning: unknown labeltype: " << labeltype << endl;break;
     }
     Utils::Put32LH(dstdata+outidx,pos);outidx+=4;
     switch (runtype) {
       case 0:Utils::Put16LH(dstdata+outidx,nrun);outidx+=2;break;
       case 1:Utils::Put32LH(dstdata+outidx,nrun);outidx+=4;break;
       default:cout << "warning: unknown runtype: " << runtype << endl;break;
     }
  }
  return outidx;
}
static void UnpackRow(uint8_t *srcdata,uint32_t width,int64_t *dstdata)
{
  memset(dstdata,0,width*sizeof(uint64_t));
  uint32_t totalruns=Utils::Get32LH(srcdata);
  int labeltype=srcdata[4];
  int runtype=srcdata[5];
  int inidx=6;
  for (uint32_t i=0;i<totalruns;i++)
  {
    uint32_t label,idx,nrun;
    switch (labeltype) {
      case 0:label=Utils::Get24LH(srcdata+inidx);inidx+=3;break;
      case 1:label=Utils::Get32LH(srcdata+inidx);inidx+=4;break;
      case 2:label=Utils::Get48LH(srcdata+inidx);inidx+=6;break;
      case 3:label=Utils::Get64LH(srcdata+inidx);inidx+=8;break;
      default:label=0;cout << "warning: unknown labeltype: " << labeltype << endl;break;
    }
    idx=Utils::Get32LH(srcdata+inidx);inidx+=4;
    switch (runtype) {
      case 0:nrun=Utils::Get16LH(srcdata+inidx);inidx+=2;break;
      case 1:nrun=Utils::Get32LH(srcdata+inidx);inidx+=4;break;
      default:nrun=0;cout << "warning: unknown runtype: " << runtype << endl;break;
    }
    for (uint32_t k=0;k<nrun;k++) dstdata[idx+k]=label;
  }
}
};*/
#endif // RLE_H
