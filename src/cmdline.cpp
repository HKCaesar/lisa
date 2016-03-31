#include "cmdline.h"
#include "common/rand.h"
#include "common/bmatrix.h"
#include "common/timer.h"
#include "analysis/map.h"

int ComandLine::OpenInputRaster(const std::string str_ifile)
{
  myIMG=nullptr;

  file=fopen(str_ifile.c_str(),"rb");
  cout << "open file: '" << str_ifile << "': ";
  if (file) {
    cout << "ok." << endl;

    if (myIMG==nullptr) {// input is pgm
      myPGM.SetHandle(file);
      if (myPGM.ReadHeader()==0) myIMG=&myPGM;
    };

    if (myIMG==nullptr) {// input is bri
      myBRI.SetHandle(file);
      if (myBRI.ReadHeader()==0) myIMG=&myBRI;
    };

    if (myIMG==nullptr) {// input is asc
      myASC.SetHandle(file);
      if (myASC.ReadHeader()==0) myIMG=&myASC;
    };

    if (myIMG==nullptr) { // use libtiff own file handling
        fclose(file);
        if (myTIFF.ReadHeader(str_ifile.c_str())==0) myIMG=&myTIFF;
    }
    return 0;
  } else return 1;
}

// connected component analysis of a supported raster file
void ComandLine::Analyze(const std::string &str_ifile,const std::string &str_bfile,AnalyzeOptions &AnalyzeOptions)
{
  if (OpenInputRaster(str_ifile)==0) {
    if (myIMG!=nullptr) {
      myIMG->PrintInfo();

      Projection myProj(myIMG->GetWidth(),myIMG->GetHeight());

      // if we have an asc file, use the coordiantes from meta-data
      if (myIMG->GetType()==IMG::TYPEASC) {
        myProj.SetProjection(myASC.getTop(),myASC.getLeft(),myASC.getBottom(),myASC.getRight());
      } else { // open coordinate file
        std::string str_pfile;
        Utils::ReplaceExt(str_ifile,str_pfile,".txt",true);
        if (myProj.ReadCoordinateFile(str_pfile)!=0) {
          cout << "missing coordinate file: '" << str_pfile << "'!" << endl;
          myIMG->Close();
          return;
        }
      }

      Timer myTimer;
      myTimer.Start();

      myProj.CalculateCellSize();
      myProj.PrintInfo();
      cout << endl;
      myProj.GenerateInterpolation();

      // read biomass file
      BM myBiomass(AnalyzeOptions.mean_biomass);

      if (str_bfile.length()) {
        Timer myTimer;
        myTimer.Start();
        if (!myBiomass.ReadAGBFile(str_bfile,AnalyzeOptions.bthres)) {
          cerr << "  warning: could not open file: '" << str_bfile << "'\n";
          return;
        }
        myTimer.Stop();
        cout << myTimer.ElapsedS() << endl;
        myBiomass.SetGeoRef(myProj.getLeft(),myProj.getTop(),myProj.getCellsize()); // setup geo-reference for Biomass-Card
      };

      BRIOptions options(*myIMG,myProj,myBiomass,AnalyzeOptions);

      if (AnalyzeOptions.write_mode>0) {
        Utils::ReplaceExt(str_ifile,options.str_clusterfile1,".tmp",true);
        Utils::ReplaceExt(str_ifile,options.str_clusterfile2,".bin",true);
        Utils::ReplaceExt(str_ifile,options.str_labelfile,".lab",true);
      }
      if (AnalyzeOptions.flush_clusters) {
        Utils::ReplaceExt(str_ifile,options.str_clusterflushfile,".clusters",true);
      }
      Cluster myCluster(options);
      myCluster.ClusterAnalyzation();

      if (AnalyzeOptions.check_consistency) myCluster.CheckClusters();

      if (AnalyzeOptions.save_mode>0) {
        std::string str_ofile;
        Utils::ReplaceExt(str_ifile,str_ofile,".csv",true);
        cout << "Saving clusters to '" << str_ofile << "'" << endl;
        if (AnalyzeOptions.save_mode==1) myCluster.SaveSmallClusterData(str_ofile);
        else if (AnalyzeOptions.save_mode==2) myCluster.SaveFullClusterData(str_ofile);
        else cout << "unknown save_mode: " << AnalyzeOptions.save_mode << endl;
      }
      myCluster.DeleteTempFiles();
      myTimer.Stop();
      cout << endl << "time: " << Utils::SecToTime(myTimer.ElapsedS()) << endl;

      myIMG->Close();
    } else cerr << "  warning: unsupported raster input format" << endl;
  } else cout << "not found!" << endl;
}


// convert given raster format pgm/asc/bri to bri
// output-pgm is currently broken
void ComandLine::Convert(const std::string &str_ifile,std::string &str_ofile,int cmode,bool globcover,bool force_overwrite,const geoExtend &myExtend)
{
  if (OpenInputRaster(str_ifile)==0) {
    if (myIMG!=nullptr) {
      myIMG->PrintInfo();

      if (cmode==ComandLine::CONVERT) {
        if (myIMG->GetType()==IMG::TYPEASC) { // save geo-referencing
          myASC.SetExtend(myExtend);
          std::string str_efile;
          Utils::ReplaceExt(str_ofile,str_efile,".txt");
          cout << "saving geo-extend to '" << str_efile << "'\n";
          if (Utils::OpenWriteCheck(str_efile,force_overwrite)) myASC.WriteExtend(str_efile,8);
        }
        IMGBRI outBRI;
        outBRI.ConvertToBRI(*myIMG,str_ofile,SIC::COMP_BILEVEL);
      }
      myIMG->Close();
    } else cerr << "  warning: unsupported raster input format" << endl;
  } else cout << "not found!" << endl;
}

int ComandLine::TestConsistence(int dimx,int dimy,double p,int verbose)
{
  int total_read_errors=0;

  cout << "test " << std::setw(4) << dimx << "x"<<std::setw(4) << dimy << ", p=" << std::fixed << std::setprecision(3) << p << ": ";
  std::string fname="test.bri";
  uint8_t linebuf[dimx];
  if (verbose) cout << "initializing random matrix of dim " << dimx << "x" << dimy << "...";
  BinMatrix m(dimy,dimx),t(dimy,dimx);
  for (int j=0;j<m.GetYDim();j++) // random matrix
    for (int i=0;i<m.GetXDim();i++)
       if (SimpleRand::BoolEvent(p)) m.Set(j,i,true);
  if (verbose) cout << "done."<< endl;
  if (verbose) m.Print();
  if (verbose) cout << "writing matrix to bri-file..." << endl;

  IMGBRI myBRI;
  myBRI.SetWidth(dimx);
  myBRI.SetHeight(dimy);
  if (myBRI.Create(fname)==0)
  {
    myBRI.Start();
    for (int j=0;j<m.GetYDim();j++)
    {
      for (int i=0;i<m.GetXDim();i++) linebuf[i]=(uint8_t)m.Get(j,i);
      myBRI.WriteRow(linebuf);
    }
    myBRI.Stop();
    myBRI.Close();
  } else cout << "error: could not create tempfile: '" << fname << "'\n";

  if (verbose) cout << "reading matrix from bri-file..." << endl;
  FILE *file=fopen(fname.c_str(),"rb");
  if (file) {
    myBRI.SetHandle(file);
    if (myBRI.ReadHeader()==0)
    {
      if (verbose) myBRI.PrintInfo();
      myBRI.Start();
      // check consitence of bri-file
      for (int j=0;j<t.GetYDim();j++)
      {
         myBRI.ReadRow(linebuf);
         for (int i=0;i<t.GetXDim();i++) t.Set(j,i,linebuf[i]);
      }
      for (int j=0;j<m.GetYDim();j++)
        for (int i=0;i<m.GetXDim();i++)
          if (m.Get(j,i)!=t.Get(j,i)) total_read_errors++;
      if (total_read_errors) cout << "error: " << total_read_errors << " elements mismatch.\n";
      myBRI.Stop();
      myBRI.SeekStart();

      BM myBiomass(0);

      Projection myProj(myBRI.GetWidth(),myBRI.GetHeight());
      myProj.SetDummyInterpolation(30); // set a dummy interpolation of 30m

      AnalyzeOptions analyze_opt;
      BRIOptions options(myBRI,myProj,myBiomass,analyze_opt);

      Cluster myCluster(options); // edge effect dept 100m
      myCluster.ClusterAnalyzation();
      cluster_stats &myStats=myCluster.GetClusterStats();

      ClusterLabel myLabel(t,30);
      cluster_stats &myStatsMatrix=myLabel.GetClusterStats();
      myLabel.LabelClusters();
      int error=0;
      int diff_num_clusters=std::abs(myStats.num_clusters-myStatsMatrix.num_clusters);
      int diff_cell_area=std::abs(myStats.cell_area-myStatsMatrix.cell_area);
      double diff_total_area=  fabs(myStats.total_area-myStatsMatrix.total_area);
      double diff_total_border=fabs(myStats.total_border_len-myStatsMatrix.total_border_len);
      double diff_max_area=fabs(myStats.max_area-myStatsMatrix.max_area);
      if (diff_num_clusters) {cout << "diff numclusters: "<<diff_num_clusters << endl;error++;};
      if (diff_cell_area) {cout << "diff cell area: "<<diff_cell_area << endl;error++;};
      if (diff_total_area>EPS)  {cout << "diff total area: "<<diff_total_area << endl;error++;};
      if (diff_total_border>EPS)  {cout << "diff total border: "<<diff_total_border << endl;error++;};
      if (diff_max_area>EPS)  {cout << "diff max area: "<<diff_max_area << endl;error++;};
      if (error) cout << "failed\n";
      else cout << "passed\n";
    } else cout << "error: tempfile no valid bri-file\n";
    fclose(file);
    remove(fname.c_str());
  } else cout << "error: could not open tempfile: '" << fname << "'\n";
  return 0;
}

void ComandLine::TestConsistency()
{
  cout << "Consistency test [numclusters,totalarea,maxarea,edgelen]...\n";
  srand(time(0));
  int dim=1000;
  TestConsistence(SimpleRand::rU_Int(1,dim),SimpleRand::rU_Int(1,dim),0,0);

  for (int i=0;i<25;i++)
    TestConsistence(SimpleRand::rU_Int(1,dim),SimpleRand::rU_Int(1,dim),SimpleRand::rU_Closed(),0);

  TestConsistence(SimpleRand::rU_Int(1,dim),SimpleRand::rU_Int(1,dim),1,0);
}

void ComandLine::CreateMap(const std::string &str_ifile,const std::string &str_ofile,int reduction_factor,int map_scale,int map_type,int edge_effect_dept,const geoExtend &myExtend)
{
  Map myDensityMap(reduction_factor,map_type,map_scale,edge_effect_dept);
  myDensityMap.CalculateMap(str_ifile,str_ofile,myExtend);
}
