#ifndef PGM_H
#define PGM_H

#include "..\global.h"
#include "imgfile.h"

class PGM : public IMGFile {
  public:
    PGM():max_gray(0),valid(false){};
    void Close();
    void PrintInfo();
    int GetMaxGray(){return max_gray;};
    off64_t GetSize(){return fsize;};
    void StartReading();
    void StopReading();
    int Create(std::string &fname);
    bool ReadRow();
    bool WriteRow();
    int ReadHeader();
    uint8_t *linebuf;
  protected:
    void ReadPGMLine(std::string &line);
    void SeekPGMLine(off64_t y);
    off64_t fsize,hdrpos;
    int max_gray;
    bool valid;
};


#endif // PGM_H
