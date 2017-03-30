#ifndef CMDLINE_H
#define CMDLINE_H

#include "global.h"
#include "analysis\cluster.h"
#include "analysis/fractal.h"
#include "file/pgm.h"
#include "file/bri.h"
#include "file/asc.h"
#include "file/tiff.h"

class ComandLine {
  public:
    enum METHOD {NONE,ANALYZE,CONVERT,INFO,MAP,CLASSIFY,REDUCE,TEST,VERSION,FRACTAL};
    void Analyze(const std::string &str_ifile,const std::string &str_bfile,const std::string &shapefile,const std::vector<int>shape_mask,AnalyzeOptions &AnalyzeOptions,const geoExtend &myextend);
    void FractalAnalysis(const std::string &str_ifile);
    void Convert(const std::string &str_ifile,std::string &str_ofile,int cmode,bool globcover,bool overwrite,const geoExtend &myExtend);
    void TestConsistency();
    int OpenInputRaster(const std::string &str_ifile);
    void Create(const std::string &str_ifile,const std::string &str_ofile,int reduction_factor,int map_type,int edge_effect_dept,const geoExtend &myExtend);
    void Reduce(const std::string &str_ifile,const std::string &str_ofile,int reduction_factor);
    void Classify(const std::string &str_ifile,const std::string &str_maskfile,const std::string &str_ofile,int mapscale,int mapclass);
  private:
    int TestConsistence(int dimx,int dimy,double p,int verbose=0);
    IMGPGM myPGM;
    IMGBRI myBRI;
    IMGASC myASC;
    IMGTIFF myTIFF;
    IMG *myIMG;
    FILE *file;
};

#endif // CMDLINE_H
