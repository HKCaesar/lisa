#ifndef PGM_H
#define PGM_H

#include "..\global.h"
#include "img.h"

class IMGPGM : public IMG {
  public:
    IMGPGM():max_gray(0),valid(false){};
    void Close();
    void PrintInfo();
    int GetMaxGray(){return max_gray;};
    off64_t GetSize(){return fsize;};
    void StartReader();
    void StopReader();
    int Create(std::string &fname);
    int ReadRow(uint8_t *buf);
    IMGTYPE GetType(){return TYPEPGM;};
    int ReadRow(){return ReadRow(rowbuffer);};
    int WriteRow();
    int ReadHeader();
  protected:
    void ReadPGMLine(std::string &line);
    void SeekPGMLine(off64_t y);
    off64_t fsize,hdrpos;
    int max_gray;
    bool valid;
};


#endif // PGM_H
