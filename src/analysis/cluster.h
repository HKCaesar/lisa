#ifndef CLUSTER_H
#define CLUSTER_H

#include "proj.h"
#include "..\file\bm.h"
#include "..\file\shape.h"
#include "..\file\img.h"

// implements connected component analyzation on .bri-files

class cluster_stats {
  public:
      cluster_stats()
      :fragment_state_total(4),fragment_state_area(4)
      {
        Reset();
      }
  void Reset(){
    cell_area=num_clusters=num_clusters10ha=num_clusters50ha=0;
    surface_area=total_area=mean_area=total_border_len=total_edge_area_de=total_edge_area_circle=0.0;
    total_biomass=total_closs=0.;

    std::fill(begin(fragment_state_total),end(fragment_state_total),0);
    std::fill(begin(fragment_state_area),end(fragment_state_area),0.0);

    max_area=std::numeric_limits<double>::min();
    min_area=std::numeric_limits<double>::max();
  };
  int64_t cell_area,num_clusters,num_clusters10ha,num_clusters50ha;
  double total_area,min_area,mean_area,total_border_len,total_edge_area_de,total_edge_area_circle,max_area;
  double total_biomass,total_closs,surface_area;
  vector <int64_t>fragment_state_total;
  vector <double>fragment_state_area;
};

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
    calc_surface_area=false;
  }
  double mean_biomass,bthres,relative_carbon_loss;
  int edge_dept,min_fragment_size,pixel_len,write_mode,save_mode;
  int forest_cover_threshold,nrows;
  bool check_consistency,flush_clusters,verbose,calc_surface_area;
};

struct tcelldata
{
  double area,border,biomass;
};


class BRIOptions
{
 public:
  BRIOptions(IMG &IMGFile,Projection &myProj,BM &myBiomass,ShapeFile &myShapeFile,AnalyzeOptions &analyze_options)
  :myIMG(IMGFile),Proj(myProj),BMass(myBiomass),SF(myShapeFile),analyze_opt(analyze_options)
  {
  };
  IMG &myIMG;
  Projection &Proj;
  BM &BMass;
  ShapeFile &SF;
  AnalyzeOptions analyze_opt;

  std::string str_labelfile,str_clusterfile1,str_clusterfile2,str_clusterflushfile;
};

class Cluster
{
  public:
    Cluster(BRIOptions &opt);
    ~Cluster();
    void ClusterAnalyzation(const geoExtend &myextent);
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
    void DetectBorders(int row,int cur_row,int i);
    void CalculateBorder(inter_cell &icell,double &border_len);
    void PrintProgress(int y,int height);
    void ProcessRow(int row,int row_offset,int cur_row,int mask_ptr);
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
    std::vector<std::vector<char>> mask_rows;
    std::vector <int64_t> cdata;
    std::vector <tcelldata> clusterdata;
    int64_t max_cluster_label,total_roots_written,num_1pixel;
    //int row;
    int64_t max_border_pixel;
    int lookahead_rows,bufrows,num_rows;
    FILE *clusterfile1;
    fstream ofs_clusterfile;
    int64_t *labelrow,*rowtmp,minlabel;
    vector <uint8_t> clusterrowdata_;
    uint32_t maxrowdatasize_;
    std::vector <bool>vborder;
};

#endif // CLUSTER_H
