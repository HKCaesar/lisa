/*
  Wrapper Class for libtiff
  able to read tiled tiffs in a scanline fashion
*/
#ifndef TIFF_H
#define TIFF_H

#include "img.h"
#include "tiffio.h"

class IMGTIFF : public IMG {
  public:
    IMGTIFF();
    int ReadHeader(const char *fname);
    void PrintInfo();
    void StartReader();
    void StopReader();
    IMGTYPE GetType(){return TYPETIFF;};
    int ReadRow(uint8_t*buf);
    int ReadRow(){return ReadRow(rowbuffer);};
    int WriteRow(){};
    void Close();
  private:
    TIFF *tif;
    uint8_t **tile_buf;
    uint32_t scanline_size,tiles_size,tilewidth,tileheight,nrow,nrow_mod;
    uint16_t comp,bps,spp;
    bool isTiled;
    tdata_t tiffbuf;
};

#endif // TIFF_H
