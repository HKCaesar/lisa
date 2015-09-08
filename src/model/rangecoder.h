#ifndef _RANGECODER_H
#define _RANGECODER_H 1

#include "..\global.h"
#include "model.h"


#define DO(n) for (uint32_t _=0;_<n;_++)

#define RANGE_ENC_NORMALIZE  while ((low ^ (low+range))<TOP || (range<BOT && ((range= -(int)low & (BOT-1)),1))) fputc(low>>24,file),range<<=8,low<<=8;
#define RANGE_DEC_NORMALIZE  while ((low ^ (low+range))<TOP || (range<BOT && ((range= -(int)low & (BOT-1)),1))) (code<<=8)+=fgetc(file),range<<=8,low<<=8;

#define SCALE_RANGE (((PSCALE-p1)*uint64_t(range)) >> PBITS) // 64 bit shift
//#define SCALE_RANGE ((uint64_t(range)*((PSCALE-p1)<<(32-PBITS)))>>32) // faster, because we only take the high word from qword

// Binary RangeCoder with Carry and 64-bit low
// derived from rc_v3 by Eugene Shelwien
class RangeCoderSH  {
  enum { NUM=4,TOP=0x01000000U,Thres=0xFF000000U};
  public:
    RangeCoderSH(FILE *iofile,int dec=0):file(iofile),decode(dec) {};
    void SetDecode(int dec){decode=dec;};
    void Init();
    void Stop();
    void EncodeBitOne(uint32_t p1,int bit);
    int  DecodeBitOne(uint32_t p1);
  protected:
    void ShiftLow();
    FILE *file;
    int decode;
    uint32_t range,code,FFNum,Cache;
    uint64_t lowc;
};

#endif
