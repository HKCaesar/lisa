#ifndef CM_H
#define CM_H

#include "counter.h"
#include "mixer.h"
#include "rangecoder.h"

class CM {
  public:
      CM(int width,FILE *file,int decode=0);
      void Init();
      void EncodeRow(uint8_t *buf);
      void DecodeRow(uint8_t *buf);
      void Stop();
      ~CM();
  protected:
    uint8_t *linebuf;
    int m_width,rowcounter;
    StatCounter Counter[8];
    RangeCoderSH Coder;
};
#endif // CM_H
