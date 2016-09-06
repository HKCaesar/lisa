#ifndef BM_H
#define BM_H

#include "asc.h"

// Read Biomass files, extension of ASC-Class
class BM : public IMGASC
{
  public:
    BM(double mean_biomass_ha)
    :mean_biomass_m2(mean_biomass_ha/10000.),data(nullptr),referenced(false),ref_left(0.),ref_top(0.),ref_cellsize(0.) {
      cout << "mean biomass="<<mean_biomass_m2 << " [t/m^2]\n";
    };
    bool ReadAGBFile(const std::string &fname,double bthres);
    ~BM() {
       if (data!=nullptr)  delete []data;
    }
    double getBiomassRef(int rx,int ry);
    void SetGeoRef(double left,double top,double cellsize)
    {
       ref_left=left;
       ref_top=top;
       ref_cellsize=cellsize;
       referenced=true;
    }
  private:
    void GeoRef(int rx,int ry,int &xpos,int &ypos);
    void ReadFile(double bthres);
    double mean_biomass_m2;
    double *data;
    bool referenced;
    double ref_left,ref_top,ref_cellsize;
};
#endif // BM_H
