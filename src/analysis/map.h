#ifndef MAP_H
#define MAP_H

#include "../global.h"
#include "../utils.h"
#include "../file/pgm.h"
#include "cluster.h"

class Map {
  struct MapLabel {
    int64_t label;
    double area;
    double edgelen;
    double closs;
  };
  public:
    Map(int reduction_factor,int maptype,int mapscale,int edge_effect_dept)
    :clusterfile_(NULL),reduction_factor_(reduction_factor),maptype_(maptype),mapscale_(mapscale),edge_effect_dept_(edge_effect_dept)
    {

    }
    void CalculateMap(const std::string &str_ifile,const std::string &str_ofile,const geoExtend &myExtend);
    static int TransformVal(double val,int nclasses);
  private:
    void ProcessRows(std::string stroutfile,int ptop,int pbottom,int pleft,int pright);
    double BlockAverage(double **rows,int width,int startx,int xblock,int yblock) const;
    void RowReduce(std::ofstream &stream,double **rows,int width,int xblock,int yblock);
    void RowReduce(IMGPGM &PGMFile,double **rows,int width,int xblock,int yblock);
    int ReadLabels(const std::string &strfile); // read all labels into "labels"
    void SkipRows(int nrows);
    void TransferLabels();
    void AllocMem();
    void FreeMem();
    int64_t maxlabel_;
    std::vector <MapLabel>labels_;
    std::vector<double>reflabels_;
    FILE *clusterfile_;
    int reduction_factor_,maptype_,mapscale_,edge_effect_dept_;
    uint32_t width_,height_;
    vector <uint8_t> rowdata_;
    int64_t *labelrow_;
    double **datarows_;
};

#endif // MAP_H
