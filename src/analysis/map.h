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
    double edgeLen;
    double closs;
  };
  public:
    Map(int reductionFactor,int mapType,int mapScale,int edgeEffectDept)
    :m_clusterFile(NULL),m_reductionFactor(reductionFactor),m_mapType(mapType),m_mapScale(mapScale),m_edgeEffectDept(edgeEffectDept)
    {

    }
    void calculateMap(const std::string &str_ifile,const std::string &str_ofile,const geoExtend &myExtend);
    static int transformVal(double val,int nclasses);
  private:
    void processRows(std::string stroutfile,int ptop,int pbottom,int pleft,int pright);
    double blockAverage(double **rows,int width,int startx,int xblock,int yblock) const;
    void rowReduce(std::ofstream &stream,double **rows,int width,int xblock,int yblock);
    void rowReduce(IMGPGM &PGMFile,double **rows,int width,int xblock,int yblock);
    int readLabels(const std::string &strfile); // read all labels into "labels"
    void skipRows(int nrows);
    void transferLabels();
    void allocMem();
    void freeMem();
    int64_t m_maxLabel;
    std::vector <MapLabel>m_labels;
    std::vector<double>m_refLabels;
    FILE *m_clusterFile;
    int m_reductionFactor,m_mapType,m_mapScale,m_edgeEffectDept;
    uint32_t m_width,m_height,m_maxRowDataSize;
    uint8_t *m_rowData;
    int64_t *m_labelRow;
    double **m_dataRows;
};

#endif // MAP_H
