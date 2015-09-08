#ifndef RLE_H
#define RLE_H

#include "global.h"

struct trun {
  int64_t label;
  uint32_t pos;
  uint32_t nrun;
};

class RLEPack {
  public:
static int PackRow(int64_t *rowdata,uint32_t width,uint8_t *dstdata)
{
  std::vector <trun> runs;
  uint32_t i=0;
  uint32_t maxnrun=0;
  int64_t maxlabel=0;
  while (i<width)
  {
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
  int labeltype=0;
  if (maxlabel<(1LL<<24)) labeltype=0;
  else if (maxlabel<(1LL<<32)) labeltype=1;
  else if (maxlabel<(1LL<<48)) labeltype=2;
  else labeltype=3;

  int runtype=0;
  if (maxnrun>=(1<<16)) runtype=1;
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
};
#endif // RLE_H
