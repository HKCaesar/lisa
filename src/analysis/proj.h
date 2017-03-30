#ifndef PROJ_H
#define PROJ_H

#include "..\global.h"
#include "..\utils.h"

// structure of an interpolation cell
struct inter_cell
{
  double pixel_width_bottom,pixel_width_top,pixel_height,pixel_area;
};

// defines functions for dealing with projections & coordinates (wgs-84)
class Projection {
  const double WGS84_a=6378137.0;//length of semi-major axis
  const double WGS84_f=1.0/298.257223563; //flattening of the ellipsoid
  const double WGS84_b=(1.0-WGS84_f)*WGS84_a; //length of semi-minor axis of the ellipsoid
//const double WGS84_e2=(WGS84_a*WGS84_a-WGS84_b*WGS84_b)/(WGS84_a*WGS84_a); //eccentricity squared
  const double WGS84_e2=(2.0-WGS84_f)*WGS84_f; //eccentricity squared
  const double WGS84_n=WGS84_f/(2.0-WGS84_f); //third flattening
  public:
    Projection(int w,int h):width(w),height(h),top(0),left(0),right(0),bottom(0){};
    int ReadCoordinateFile(const std::string &fname);
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
      return inter_matrix[y];
    }
    void GenerateInterpolation();
    void SetDummyInterpolation(int len);
    //double GetMeanPixelArea(){return mean_pixelarea;};
  protected:
    //inter_cell CalcArea_geodesic(double top,double left,double bottom,double right);
    void CalcArea_trapezoid_parallel(double top,double left,double bottom,double right,double &sheight,double &top_width,double &bottom_width,double &area);
    inter_cell CalcPixelArea_trapezoid(double top);
    inter_cell CalcPixelArea_exact(double top);
    double CalcArea_rectangle(double top,double left,double right);

    void GetLatLong(int x,int y,double &geo_width,double &geo_height);
    double CalcDeltaLong(double lat,double dLong);
    double CalcDeltaLat(double lat);
    double CalcDeltaLat(double lat1,double lat2);
    int GetPPosLong(double dLong);
    double GetLong(int w);
    double CalcDist_Haversine(double lat1,double long1,double lat2,double long2,double precision=0);
    double CalcDist_Vincenty(double lat1,double long1,double lat2,double long2,double precision=0);
    int width,height,numcells;
    double top,left,right,bottom;
    double cellsize;
    std::vector < inter_cell> inter_matrix;
};

#endif // PROJ_H
