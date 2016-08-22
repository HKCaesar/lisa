#ifndef IMGFILE_H
#define IMGFILE_H

#include "..\global.h"

//base class for raster data
class IMG {
  public:
    enum IMGTYPE {UNKNOWN,TYPEPGM,TYPEBRI,TYPEASC,TYPETIFF};
    IMG():file(nullptr),width(0),height(0),fsize(0){};
    bool OpenRead(const std::string &fname) {
      if ( (file=fopen(fname.c_str(),"rb"))!=NULL) return true;
      else return false;
    }
    virtual void PrintInfo()=0;
    virtual IMGTYPE GetType()=0;
    void SetHandle(FILE *f){file=f;};
    void SetWidth(int w){width=w;};
    void SetHeight(int h){height=h;};
    virtual void StartReader()=0;
    virtual void StopReader()=0;
    virtual int ReadRow(uint8_t *buf)=0;
    virtual int ReadRow()=0;
    int GetWidth(){return width;};
    int GetHeight(){return height;};
    uint8_t *rowbuffer;
    virtual ~IMG(){Close();};
    void Close(){
      if (file!=nullptr) {fclose(file);file=nullptr;};
    }
  protected:

    FILE *file;
    int width,height;
    off64_t fsize;
};
#endif // IMGFILE_H
