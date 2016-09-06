#ifndef ASC_H
#define ASC_H

#include "img.h"
#include "..\utils.h"

// Functions for reading esri ASCII raster data

template <typename T>
 bool SplitTokens(const std::string &str,std::vector<T> &vec)
 {
  const char *pstr=(char*)str.c_str();
  const char delim=' ';
  const int  slen=(int)str.length();
  const int  vlen=(int)vec.size();
  char token[1025];
  int tidx=0;
  int sidx=0;
  int vidx=0;
  while (sidx<slen) {
    while (pstr[sidx]==delim && sidx<slen) sidx++;
    while (pstr[sidx]!=delim && sidx<slen) {
        if (tidx==1024) cout << "warning: token too long\n";
        else token[tidx++]=pstr[sidx];
        sidx++;
    };
    if (tidx) {
       token[tidx]='\0';
       if (vidx<vlen) {
            if (std::is_same<T,int>::value) vec[vidx]=atoi(token);
            else if (std::is_same<T,double>::value) vec[vidx]=atof(token);
       }
       else cout << "warning: array access out of bounds!\n";
       //if (std::string(token).compare("0.0")==0) cout << vidx << endl;
       vidx++;
       tidx=0;
    };
  }
  return (vidx==vlen);
 };

class IMGASC : public IMG
{
  public:
    IMGASC():nodata_avail(false){};
    double getTop(){return ttop;};
    double getLeft(){return tleft;};
    double getBottom(){return tbottom;};
    double getRight(){return tright;};
    int ReadHeader();

    IMGTYPE GetType(){return TYPEASC;};

    void PrintInfo();
    int WriteExtend(const std::string &fname,int prec);
    void StartReader();
    void StopReader();
    void SetExtend(const geoExtend &myExtend);
    int ReadRow(uint8_t *buf);
    int ReadRow(){return ReadRow(rowbuffer);};
    int getExtendWidth(){return pright-pleft;};
    int getExtendHeight(){return pbottom-ptop;};
  protected:
    bool getLine(std::string &line);
    double top,left,right,bottom,cellsize;
    std::string line;
    int linenum;
    bool nodata_avail;
    int nodata_value;
  private:
    std::vector <int>vecline;
    int getLongPos(double ref_pos,double left);
    int getLatPos(double ref_pos,double top);

    double ttop,tleft,tright,tbottom;
    int pleft,pright,ptop,pbottom;
    int64_t hist[2];
};

#endif // ASC_H
