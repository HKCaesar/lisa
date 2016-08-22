#ifndef PGM_H
#define PGM_H

#include "..\global.h"
#include "img.h"

/*
  PGM ImageReader/Writer
  supports: 8-Bit Grayscale Reading/Writing, 24-Bit RGB Writing
*/

class IMGPGM : public IMG {
  public:
    enum PGMTYPE {PGMGRAY,PGMRGB};
    IMGPGM():max_gray(0),scanlinesize_(0),pgmtype_(PGMGRAY){};
    void PrintInfo();
    int GetMaxGray(){return max_gray;};
    off64_t GetSize(){return fsize;};
    void StartReader();
    void StopReader();
    void SetType(PGMTYPE pgmtype){pgmtype_=pgmtype;};
    int Create(const std::string &fname);
    int ReadRow(uint8_t *buf);
    IMGTYPE GetType(){return TYPEPGM;};
    int GetScanLineSize(){return scanlinesize_;};
    int ReadRow(){return ReadRow(rowbuffer);};
    int WriteRow();
    int ReadHeader();
  protected:
    void SetScanLineSize() {
      if (pgmtype_==PGMGRAY) scanlinesize_=width;
      else if (pgmtype_==PGMRGB) scanlinesize_=3*width;
      else scanlinesize_=0;
    }
    void ReadPGMLine(std::string &line);
    void SeekPGMLine(off64_t y);
    off64_t fsize,hdrpos;
    int max_gray,scanlinesize_;
    PGMTYPE pgmtype_;
};


#endif // PGM_H
