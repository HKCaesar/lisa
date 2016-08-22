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
  struct RowLabel {
    double data;
    bool labeled;
  };
  public:
    Map(int reduction_factor,int maptype,int mapscale,int mapclass,int edge_effect_dept)
    :clusterfile_(NULL),reduction_factor_(reduction_factor),maptype_(maptype),mapscale_(mapscale),mapclass_(mapclass),edge_effect_dept_(edge_effect_dept)
    {

    }
    void Create(const std::string &str_ifile,const std::string &str_ofile,const geoExtend &myExtend);
    void Reduce(const std::string &str_ifile,const std::string &str_ofile);
    void Classify(const std::string &str_ifile,const std::string &str_maskfile,const std::string &str_ofile);
  private:
    static int TransformVal(double val,int offset,int nclasses);
    void PrintColorTable(int offset);
    double BlockAverageClassified(int width,int startx,int xblock,int yblock) const;
    int BlockAverageClass(const vector <vector<uint8_t>> &rows,int width,int startx,int xblock,int yblock) const;
    void RowReduceClass(IMGPGM &PGMFile,const vector <vector<uint8_t>> &rows,int width,int xblock,int yblock);
    int GetMapWidth(int pwidth);
    int GetMapHeight(int pwidth);
    void ProcessRows(const std::string &stroutfile,int ptop,int pbottom,int pleft,int pright);
    //double BlockAverage(double **rows,int width,int startx,int xblock,int yblock) const;
    void RowReduce(std::ofstream &stream,int width,int xblock,int yblock);
    void RowReduce(IMGPGM &PGMFile,int width,int xblock,int yblock);
    int ReadLabels(const std::string &strfile); // read all labels into "labels"
    void SkipRows(int nrows);
    void TransferLabels();
    void AllocMem();
    void FreeMem();
    int64_t maxlabel_;
    std::vector <MapLabel>labels_;
    std::vector<double>reflabels_;
    FILE *clusterfile_;
    int reduction_factor_,maptype_,mapscale_,mapclass_,edge_effect_dept_;
    uint32_t width_,height_;
    vector <uint8_t> rowdata_;
    int64_t *labelrow_;
    vector <vector<RowLabel>>datarows_;
    //double **datarows_;
};

#endif // MAP_H
