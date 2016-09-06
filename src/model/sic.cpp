#include "sic.h"
#include "..\utils.h"

namespace RLEPack {
      void EncodeLabel(BitBufferSafe &bitout,int64_t diff,int nrun)
      {
        bitout.PutEliasGamma(std::abs(diff)+2);
        bitout.PutBit(diff<0);
        bitout.PutEliasGamma(nrun+1);
      }
      int PackRow(int64_t *rowdata,uint32_t width,vector <uint8_t> &dstdata)
      {
        BitBufferSafe bitout(dstdata);
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
      void UnpackRow(vector <uint8_t> &srcdata,uint32_t width,int64_t *dstdata)
      {
        memset(dstdata,0,width*sizeof(int64_t));
        BitBufferSafe bitin(srcdata);
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
             if (i+nrun>width) cerr << "rlepack: attempt to read over eol\n";
             else for (int k=0;k<nrun;k++) dstdata[i+k]=label;
             llabel=label;
          } else nrun=bitin.GetEliasGamma();
          i+=nrun;
        }
      }
};


SIC::SIC(int img_width,COMP_TYPE compression_type)
:width(img_width),comptype(compression_type)
{
  //cout << "sic: " << width << ", " << comptype << endl;
  tbuf=new uint8_t[10*width];
}

SIC::~SIC()
{
  delete []tbuf;
}

int SIC::CompressRowBinary(uint8_t *linebuf,uint8_t *outbuf)
{
  BitBuffer bitout(outbuf+4);
  bool sbit=linebuf[0];
  bitout.PutBits(sbit,1);
  int nrun=0;

  for (int i=1;i<width;i++)
  {
     const int val=linebuf[i];

     if (val!=0 && val!=1) {cout << "sic: invalid input value: " << (int)val << endl;return 0;};
     if (val==sbit) nrun++;
     else {
        bitout.PutEliasGamma(nrun+1);
        nrun=0;
        sbit=!sbit; // flip bit
     }
  }
  bitout.PutEliasGamma(nrun+1);
  bitout.Flush();
  Utils::Put32LH(outbuf,bitout.GetBytesProcessed());
  return bitout.GetBytesProcessed()+4;
}

int SIC::DecompressRowBinary(uint8_t *inbuf,uint8_t *linebuf)
{
  BitBuffer bitin(inbuf+4);
  bool sbit=bitin.GetBits(1);

  //int N=1,A=256;
  int i=0;
  while (i<width) {
    //int nrun=bitin.GetRice(BitBuffer::EstimateK(N,A));
    int nrun=bitin.GetEliasGamma();
    for (int r=0;r<nrun;r++) linebuf[i++]=sbit;
    //if (N>=256) {N>>=1;A>>=1;};
    //N++;A+=nrun;
    sbit=!sbit;
  }
  return bitin.GetBytesProcessed()+4;
}

void SIC::EncodeVal(BitBuffer &bitbuf,int val,int pred)
{
  int err=val-pred;
  if (err>=0) err=2*abs(err); // remap the error
  else err=2*abs(err)-1;

  /*int k=bitbuf.EstimateK(N,A);
  bitbuf.PutRice(err,k);*/

  /*A+=err;
  N+=1;

  if (N>256) {
    A>>=1;
    N>>=1;
  }*/
  bitbuf.PutEliasGamma(err+1);
}

int SIC::DecodeVal(BitBuffer &bitbuf,int pred)
{
  int err=bitbuf.GetEliasGamma()-1;
  if (err&1) err=-((err+1)/2);
  else err=err/2;

  return (err+pred);
}

int SIC::CompressRowGrey(uint8_t *linebuf,uint8_t *outbuf)
{
  BitBuffer bitout(outbuf+4);

  bool runmode=false;
  int runlen=0;
  vector <int>rhist(256);

  for (int i=0;i<width;i++) {
    int x=linebuf[i];
    int w=0,ww=0;
    if (i>0) w=linebuf[i-1];
    if (i>1) ww=linebuf[i-2];

    if (!runmode && w==ww) {runmode=1;runlen=0;};

    if (!runmode) {
      EncodeVal(bitout,x,w);
    } else {
      if (x==w) runlen++;
      else {
        bitout.PutEliasGamma(runlen+1);
        if (runlen<256) rhist[runlen]++;
        EncodeVal(bitout,x,w);
        runmode=false;
      }
    }
  }
  /*for (int i=0;i<32;i++) cout << rhist[i] << " ";
  cout << endl;*/
  if (runmode) bitout.PutEliasGamma(runlen+1);
  bitout.Flush();
  Utils::Put32LH(outbuf,bitout.GetBytesProcessed());
  return bitout.GetBytesProcessed()+4;
}

int SIC::DecompressRowGrey(uint8_t *inbuf,uint8_t *linebuf)
{
  BitBuffer bitin(inbuf+4);
  int i=0;
  bool runmode=false;
  while (i<width) {
    int w=0,ww=0;
    if (i>0) w=linebuf[i-1];
    if (i>1) ww=linebuf[i-2];

    if (!runmode && w==ww) runmode=true;

    if (!runmode) {
      linebuf[i]=DecodeVal(bitin,w);
      i++;
    } else {
      int runlen=bitin.GetEliasGamma()-1;
      for (int k=0;k<runlen;k++) linebuf[i+k]=w;
      i+=runlen;
      linebuf[i]=DecodeVal(bitin,w);
      i++;
      runmode=false;
    }
  }
  return bitin.GetBytesProcessed()+4;
}
