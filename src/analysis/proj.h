#ifndef PROJ_H
#define PROJ_H

#include "..\global.h"

// structure of an interpolation cell
struct inter_cell
{
  double pixel_width,pixel_height,pixel_area;
};

// defines functions for dealing with projections & coordinates (wgs-84)
class Projection {
  public:
    Projection(int w,int h):width(w),height(h),top(0),left(0),right(0),bottom(0),mean_pixelarea(0.0){};
    int ReadCoordinateFile(std::string &fname);
    double getLeft(){return left;};
    double getTop(){return top;};
    double getRight(){return right;};
    double getBottom(){return bottom;};
    double getCellsize(){return cellsize;};
    void PrintInfo() {
      GeoUtils::PrintInfo(top,left,right,bottom,cellsize);
    };
    void CalculateCellSize()
    {
      ycell=delta_lat/(double)height;
      xcell=delta_long/(double)width;
      if (fabs(xcell-ycell)>EPS)
      {
        cout << "warning: x/y cellsize different" << endl;
        cout << "xcell: " << xcell << endl;
        cout << "ycell: " << ycell << endl;
      }
      cellsize=xcell;
      cout << "cellsize: " << cellsize << " degree = " << GeoUtils::cellsize2arcsec(cellsize) << " arcsec = " << GeoUtils::cellsize2m(cellsize) << " m" << endl;
    }
    inter_cell &GetCellDim(int y)
    {
      return inter_matrix[y/dy];
    }
    void GenerateInterpolation(int cells);
    void GetLatLong(int x,int y,double &geo_width,double &geo_height);
    void SetDummyInterpolation(int len);
    double GetMeanPixelArea(){return mean_pixelarea;};
  protected:
    double CalcDist_Haversine(double lat1,double long1,double lat2,double long2);
    double CalcDist_Vincenty(double lat1,double long1,double lat2,double long2);
    void CalculateGeoDelta()
    {
      delta_lat=top-bottom;
      delta_long=right-left;
      cout << "dLat: " << delta_lat << " degree, dLong: " << delta_long << " degree" << endl;
    };
    int width,height,numcells;
    double top,left,right,bottom;
    double delta_lat; // geographische breite
    double delta_long;// geographische länge
    double xcell,ycell,cellsize,mean_pixelarea;
    int dy;
    std::vector < inter_cell> inter_matrix;
};

#endif // PROJ_H
