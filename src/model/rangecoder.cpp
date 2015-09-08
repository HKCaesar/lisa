#include "RangeCoder.h"

// binary rangecoder
void RangeCoderSH::Init()
{
  range = 0xFFFFFFFF;
  lowc = FFNum = Cache = code = 0;
  if(decode==1) DO(NUM+1) (code <<=8) += fgetc(file);
}

void RangeCoderSH::Stop()
{
  if (decode==0) DO(NUM+1) ShiftLow();
}

void RangeCoderSH::EncodeBitOne(uint32_t p1,int bit)
{
  const uint32_t rnew=SCALE_RANGE;
  bit ? range-=rnew, lowc+=rnew : range=rnew;
  while(range<TOP) range<<=8,ShiftLow();
}

int RangeCoderSH::DecodeBitOne(uint32_t p1)
{
  const uint32_t rnew=SCALE_RANGE;
  int bit = (code>=rnew);
  bit ? range-=rnew, code-=rnew : range=rnew;
  while(range<TOP) range<<=8,(code<<=8)+=fgetc(file);
  return bit;
}

void RangeCoderSH::ShiftLow()
{
  uint32_t Carry = uint32_t(lowc>>32), low = uint32_t(lowc);
  if( low<Thres || Carry )
  {
     fputc( Cache+Carry, file);
     for (;FFNum != 0;FFNum--) fputc(Carry-1,file);
     Cache = low>>24;
   } else FFNum++;
  lowc = (low<<8);
}
