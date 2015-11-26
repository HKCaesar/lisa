#ifndef CLUSTER_H
#define CLUSTER_H

#include "proj.h"
#include "..\file\bm.h"
#include "..\file\bri.h"
#include "..\rle.h"

// implements connected component analyzation on .bri-files

struct tcelldata
{
  double area,border,biomass;
};

class BRIOptions
{
 public:
  BRIOptions(BRI &myBRI,Projection &myProj,BM &myBiomass)
  :BRIFile(myBRI),Proj(myProj),BMass(myBiomass)
  {
    pixel_len=1;
    write_clusterlabel=0;
    edge_dept=100.0;
    min_fragment_size=0.0;
    relative_carbon_loss=0.5;
    verbose=false;
  };
  BRI &BRIFile;
  Projection &Proj;
  BM &BMass;
  int pixel_len,write_clusterlabel;
  double edge_dept,relative_carbon_loss,min_fragment_size;
  bool verbose;
  std::string str_labelfile,str_clusterfile1,str_clusterfile2;
};

class ClusterBRI
{
  public:
    ClusterBRI(BRIOptions &opt);
    ~ClusterBRI();
    void ClusterAnalyzation();
    void SaveFullClusterData(std::string &fname);
    void SaveSmallClusterData(std::string &fname);
    void CheckClusters();
    cluster_stats &GetClusterStats(){return myStats;};
  protected:
    cluster_stats myStats;
    double CalculateCLossPerHA(int64_t label);
    double CalculateCLoss(double biomass,double area_m2,double edge_area);
    void WriteLabelFile();
    void WriteClusterfile();
    void WriteMarkedRow(int64_t *clusterow,uint32_t width,FILE *file);
    void DetectBorders(int row,int cur_row,int i,bool &bleft,bool &bright,bool &btop,bool &bbottom);
    double CalculateBorder(inter_cell &icell,bool left,bool right,bool top,bool bottom,int &border_pixel);
    void PrintProgress(int y,int height);
    void ProcessRow(int row,int cur_row);
    void CalculateStats();
    void PrintHist(std::vector <int64_t> &hist,std::string header);
    void PrintHist(std::vector <double> &hist,std::string header,std::string unit);
    void WriteHist(ofstream &file,std::vector <int64_t> &hist,std::string header);
    void WriteHist(ofstream &file,std::vector <double> &hist,std::string header);
    double CalculateEdgeAreaDE(double area,double edge_len);
    double CalculateEdgeAreaCircle(double area);
    int64_t FindRoot(int64_t clabel);
    int64_t FindCollapse(int64_t label);
    int64_t GetNumRoots();
    BRIOptions opt;

    int bri_width;
    uint8_t *rowbuffer;
    int64_t **wrows; // prev, cur, next
    std::vector <int64_t> cdata;
    std::vector <tcelldata> clusterdata;
    int64_t max_cluster_label,num_1pixel;
    int row;
    int64_t max_border_pixel;
    int lookahead_rows,bufrows,endrow;
    FILE *clusterfile1;
    int64_t *labelrow,*rowtmp,minlabel;
    uint8_t *clusterrowdata;
};

#endif // CLUSTER_H
