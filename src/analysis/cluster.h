#ifndef CLUSTER_H
#define CLUSTER_H

#include "proj.h"
#include "..\file\bm.h"
#include "..\file\img.h"

// implements connected component analyzation on .bri-files

class AnalyzeOptions {
  public:
  AnalyzeOptions()
  {
    bthres=0.; // biomass threshold: 0 t/ha
    mean_biomass=0.;
    relative_carbon_loss=0.5; // relative carbon loss in edge areas

    //ncells=1; // number of interpolation cells
    edge_dept=100; // edge effect dept 100m
    min_fragment_size=0; // minimum fragment size 0 ha
    pixel_len=1; // pixel number for edge detection

    write_mode=0;
    save_mode=0;

    forest_cover_threshold=0;
    nrows=0;

    check_consistency=false;
    flush_clusters=false;
    verbose=true;
  }
  double mean_biomass,bthres,relative_carbon_loss;
  int edge_dept,min_fragment_size,pixel_len,write_mode,save_mode;
  int forest_cover_threshold,nrows;
  bool check_consistency,flush_clusters,verbose;
};

struct tcelldata
{
  double area,border,biomass;
};


class BRIOptions
{
 public:
  BRIOptions(IMG &IMGFile,Projection &myProj,BM &myBiomass,AnalyzeOptions &analyze_options)
  :myIMG(IMGFile),Proj(myProj),BMass(myBiomass),analyze_opt(analyze_options)
  {
  };
  IMG &myIMG;
  Projection &Proj;
  BM &BMass;
  AnalyzeOptions analyze_opt;

  std::string str_labelfile,str_clusterfile1,str_clusterfile2,str_clusterflushfile;
};

class Cluster
{
  public:
    Cluster(BRIOptions &opt);
    ~Cluster();
    void ClusterAnalyzation();
    void SaveFullClusterData(std::string &fname);
    void SaveSmallClusterData(std::string &fname);
    void CheckClusters();
    void DeleteTempFiles();
    cluster_stats &GetClusterStats(){return myStats;};
    static double CalculateEdgeAreaDE(double area,double edge_len,double edge_effect_dept);
  protected:
    void FlushClusters(int cur_row);
    void AddClusterStats(int64_t parea,const tcelldata &cell);
    void AddClusterSmallStats(const tcelldata &cell,vector <int64_t>&hist_area,vector <double>&hist_totalarea,vector <double>&hist_totaledge,vector <double>&hist_biomass,vector <double>&hist_totalloss);
    int UnpackRow(int64_t *dstrow,uint8_t *srcrow,int len);
    cluster_stats myStats;
    void CompressTree(int cur_row);
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
    //uint8_t *rowbuffer;
    int64_t **wrows; // prev, cur, next
    std::vector <int64_t> cdata;
    std::vector <tcelldata> clusterdata;
    int64_t max_cluster_label,total_roots_written,num_1pixel;
    int row;
    int64_t max_border_pixel;
    int lookahead_rows,bufrows,endrow;
    FILE *clusterfile1;
    fstream ofs_clusterfile;
    int64_t *labelrow,*rowtmp,minlabel;
    uint8_t *clusterrowdata;
};

#endif // CLUSTER_H
