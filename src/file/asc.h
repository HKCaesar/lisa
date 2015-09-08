#ifndef ASC_H
#define ASC_H

#include "imgfile.h"

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

class ASC : public IMGFile
{
  public:
    ASC(){};
    int ReadHeader();
    void PrintInfo();
    int WriteExtend(const std::string &fname,int prec);
    void StartReading();
    void StopReading();
    void SetExtend(const geoExtend &myExtend);
    uint8_t *linebuf;
    int ReadRow();
    int getExtendWidth(){return pright-pleft;};
    int getExtendHeight(){return pbottom-ptop;};
  protected:
    bool getLine(std::string &line);
    double top,left,right,bottom,cellsize;
    std::string line;
    int linenum;
  private:
    std::vector <int>vecline;
    int getLongPos(double ref_pos,double left);
    int getLatPos(double ref_pos,double top);

    double ttop,tleft,tright,tbottom;
    int nodataval;
    int pleft,pright,ptop,pbottom;
    int64_t hist[2];
};

#endif // ASC_H
