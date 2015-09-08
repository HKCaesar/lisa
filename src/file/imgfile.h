#ifndef IMGFILE_H
#define IMGFILE_H

#include "..\global.h"

class IMGFile {
  public:
    IMGFile():width(0),height(0),fsize(0){};
    bool OpenRead(const std::string &fname) {
      if ( (file=fopen(fname.c_str(),"rb"))!=NULL) return true;
      else return false;
    }
    void PrintInfo(){};
    void SetHandle(FILE *f){file=f;};
    void SetWidth(int w){width=w;};
    void SetHeight(int h){height=h;};
    int GetWidth(){return width;};
    int GetHeight(){return height;};
  protected:
    FILE *file;
    int width,height;
    off64_t fsize;
};
#endif // IMGFILE_H
