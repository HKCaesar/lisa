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
  const double WGS84_a=6378137.0;//length of semi-major axis
  const double WGS84_f=1.0/298.257223563; //flattening of the ellipsoid
  const double WGS84_b=(1.0-WGS84_f)*WGS84_a; //length of semi-minor axis of the ellipsoid
  public:
    Projection(int w,int h):width(w),height(h),top(0),left(0),right(0),bottom(0){};
    int ReadCoordinateFile(std::string &fname);
    void SetProjection(double ttop,double tleft,double tbottom,double tright);
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
      double delta_lat=top-bottom;
      double delta_long=right-left;

      double ycell=delta_lat/(double)height;
      double xcell=delta_long/(double)width;
      if (fabs(xcell-ycell)>EPS) {
        cout << "warning: x/y cellsize different" << endl;
        cout << "xcell: " << xcell << endl;
        cout << "ycell: " << ycell << endl;
      }
      cellsize=xcell;
      //cout << "cellsize: " << cellsize << " degree = " << GeoUtils::cellsize2arcsec(cellsize) << " arcsec = " << GeoUtils::cellsize2m(cellsize) << " m" << endl;
    }
    inter_cell &GetCellDim(int y)
    {
      return inter_matrix[y/dy];
    }
    void GenerateInterpolation(int cells);
    void SetDummyInterpolation(int len);
    //double GetMeanPixelArea(){return mean_pixelarea;};
  protected:
    void GetLatLong(int x,int y,double &geo_width,double &geo_height);
    double CalcDeltaLong(double lat);
    int GetPPosLong(double dLong);
    double GetLong(int w);
    double CalcDist_Haversine(double lat1,double long1,double lat2,double long2,double precision=0);
    double CalcDist_Vincenty(double lat1,double long1,double lat2,double long2,double precision=0);
    int width,height,numcells;
    double top,left,right,bottom;
    double cellsize;
    int dy;
    std::vector < inter_cell> inter_matrix;
};

#endif // PROJ_H
