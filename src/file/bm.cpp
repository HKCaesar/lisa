#include "bm.h"

void BM::GeoRef(int rx,int ry,int &xpos,int &ypos)
{
  xpos=ypos=0;
  double ref_long=ref_left+(double)rx*ref_cellsize;
  double ref_lat=ref_top-(double)ry*ref_cellsize;

  xpos=round((ref_long-left)/cellsize);
  ypos=round((top-ref_lat)/cellsize);
}

double BM::getBiomassRef(int rx,int ry)
{
  if (referenced) {
    int x,y;
    GeoRef(rx,ry,x,y);
    if (x < 0 || x>=width) return mean_biomass_m2;
    if (y < 0 || y>=height) return mean_biomass_m2;
    double biomass_ha=data[y*width+x];
    if (biomass_ha<0) return mean_biomass_m2;
    return biomass_ha/10000.;
  } else return mean_biomass_m2;
}

void BM::ReadFile(double bthres)
{
    cout << "biomass threshold: " << bthres << " t/ha" << endl;

    // stats
    int64_t num_pixel_valid,num_pixel_total;
    num_pixel_valid=num_pixel_total=0;
    double val_total=0.;
    double min_value,max_value;
    max_value=std::numeric_limits<double>::min();
    min_value=std::numeric_limits<double>::max();

    std::vector <double>vecline;
    vecline.resize(width);

    linenum=0;
    while (linenum<height)
    {
        if (!getLine(line)) cout << "error reading file at line: "<<linenum << endl;
        else if (!SplitTokens<double>(line,vecline)) cout << "error parsing file at line: " << linenum << endl;
        else {
          for (size_t i=0;i<vecline.size();i++) {
            num_pixel_total++;
            double val=vecline[i];
            if (val>bthres) {
               data[linenum*width+i]=val;

               if (val<min_value) min_value=val;
               else if (val>max_value) max_value=val;

               num_pixel_valid++;val_total+=val;
            }
            else data[linenum*width+i]=mean_biomass_m2;
          }
        }
        linenum++;
        if (linenum%100==0) cout << "Reading: " << Utils::ConvertFixed(linenum*100/(double)height,1) << "%\r";
    }
    cout << "Reading: " << Utils::ConvertFixed(linenum*100/(double)height,1) << "%\n";
    cout << "total pixels: " << num_pixel_total << " ("<<width<<"x"<<height<<"="<<width*height<<")\n";
    cout << "valid pixels: " << num_pixel_valid << " ("<<Utils::ConvertFixed(num_pixel_valid*100./(double)num_pixel_total,1)<<"%)\n";
    cout << "mean: " << (val_total/(double)num_pixel_valid) << " t/ha\n";
    cout << "min : " << min_value << "\n";
    cout << "max : " << max_value << "\n";
}

bool BM::ReadAGBFile(const std::string &fname,double bthres)
{
  if (OpenRead(fname)) {
    cout << "reading AGB-file: '" << fname << "'" << endl;
    if (ReadHeader()==0) {
      PrintInfo();
      data=new double[width*height];
      ReadFile(bthres);
      fclose(file);
      return true;
    } else return false;
  } else return false;
}
