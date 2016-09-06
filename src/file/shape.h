#ifndef SHAPE_H
#define SHAPE_H

#include "asc.h"

// read classification files, can be merged with BM-class
class ShapeFile : public IMGASC
{
  public:
    ShapeFile(const std::vector <int>&vec_classes)
    :vec_classes(vec_classes),data(nullptr),referenced(false),ref_left(0.),ref_top(0.),ref_cellsize(0.) {
    };
    bool ReadShapeFile(const std::string &fname);
    ~ShapeFile() {
       if (data!=nullptr)  delete []data;
    }
    bool IsCLASS(int rx,int ry);
    void SetGeoRef(double left,double top,double cellsize)
    {
       ref_left=left;
       ref_top=top;
       ref_cellsize=cellsize;
       referenced=true;
    }
    bool IsReferenced(){return referenced;};
  private:
    int readToken(double &val);
    void GeoRef(int rx,int ry,int &xpos,int &ypos);
    void ReadFile();
    std::vector <int>vec_classes;
    int *data;
    bool referenced;
    double ref_left,ref_top,ref_cellsize;
};
#endif



