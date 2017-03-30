#ifndef SIC_H
#define SIC_H

#include "../global.h"
#include "../file/bitio.h"

// sic (simple image coder) v0.1

// simple bitpack-compression for clusterfiles
namespace RLEPack {
      void EncodeLabel(BitBufferSafe &bitout,int64_t diff,int nrun);
      int PackRow(int64_t *rowdata,uint32_t width,vector <uint8_t> &dstdata);
      void UnpackRow(vector <uint8_t> &srcdata,uint32_t width,int64_t *dstdata);
};

class SIC {
  public:
    enum COMP_TYPE {COMP_NONE,COMP_BILEVEL,COMP_GRAY};
    SIC(int img_width,COMP_TYPE compression_type);
    ~SIC();
    int CompressRowBinary(uint8_t *linebuf,uint8_t *outbuf);
    static int CompressRowBinary(uint8_t *linebuf,uint8_t *outbuf,int width);
    static int DecompressRowBinary(uint8_t *inbuf,uint8_t *linebuf,int width);
    int DecompressRowBinary(uint8_t *inbuf,uint8_t *linebuf);
    int CompressRowGrey(uint8_t *linebuf,uint8_t *outbuf);
    int DecompressRowGrey(uint8_t *inbuf,uint8_t *linebuf);
  private:
    void EncodeVal(BitBuffer &bitbuf,int val,int pred);
    int DecodeVal(BitBuffer &bitbuf,int pred);
    int width;
    COMP_TYPE comptype;
    uint8_t *tbuf;
};
#endif // SIC_H
