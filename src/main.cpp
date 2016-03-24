#include "global.h"
#include "rand.h"
#include "file/pgm.h"
#include "file/bri.h"
#include "file/asc.h"
#include "file/tiff.h"
#include "file/bm.h"
#include "bmatrix.h"
#include "analysis/cluster.h"
#include "analysis/proj.h"
#include "model/rangecoder.h"
#include "model/counter.h"
#include "model/vle.h"
#include "cmdoptions.h"

void TestCM(std::string fname)
{
  FILE *file=fopen(fname.c_str(),"wb");
  RangeCoderSH myCoder(file);
  if (file) {
    myCoder.Init();
    myCoder.EncodeBitOne(10,1);
    myCoder.EncodeBitOne(9,0);
    myCoder.EncodeBitOne(8,1);
    myCoder.Stop();
    fclose(file);
    file=fopen(fname.c_str(),"rb");
    myCoder.SetDecode(1);
    myCoder.Init();
    int p;
    p=myCoder.DecodeBitOne(10);
    cout << p;
    p=myCoder.DecodeBitOne(9);
    cout << p;
    p=myCoder.DecodeBitOne(8);
    cout << p;
    myCoder.Stop();
    fclose(file);
  } else cout <<"could not open\n";
  //remove(fname.c_str());
}

int TestConsistence(int dimx,int dimy,double p,int verbose=0)
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

      ClusterBRI myCluster(options); // edge effect dept 100m
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


void SplitPath(const std::string &str,std::string &path,std::string &fname)
{
  unsigned found = str.find_last_of("/\\");
  path=str.substr(0,found);
  fname=str.substr(found+1);
}

double BlockReduce(double **rows,int width,int startx,int xblock,int yblock)
{
  double sum=0.0;
  int64_t k=0;
  for (int j=0;j<yblock;j++) {
    for (int i=startx;i<std::min(width,startx+xblock);i++) {
      sum+=(double)rows[j][i];
      k++;
    }
  }
  return sum/(double)k;
}

void RowReduce(ofstream &stream,double **rows,int width,int xblock,int yblock)
{
  int nblocks=ceil(width/(double)xblock);
  for (int block=0;block<nblocks;block++)
  {
    double avg=BlockReduce(rows,width,block*xblock,xblock,yblock);
    stream << std::to_string(avg);
    if (block<nblocks-1) stream<<",";
  }
  stream << std::endl;
}

void PrintVersion(int mode=0)
{
  cout << LISA_VERSION << endl;
  if (mode) {
  cout << "compiled ";
  #ifdef __GNUC__
    cout << "with GNU " << __VERSION__;
  #endif
  #if __x86_64__
    cout << " (64-bit)";
  #else
    cout << " (32-bit)";
  #endif
  cout << " on " << __DATE__ << endl;
  }
  #ifdef TIFF_SUPPORT
    cout << endl << TIFFGetVersion() << endl;
  #endif
}



class ComandLine {
  struct tlabel {
    int64_t label;
    double area;
    double edge_len;
    double closs;
  };
  public:
    enum METHOD {NONE,ANALYZE,CONVERT,INFO,MAP,TEST,VERSION};
    void Analyze(const std::string &str_ifile,const std::string &str_bfile,AnalyzeOptions &AnalyzeOptions);
    void Convert(const std::string &str_ifile,std::string &str_ofile,int cmode,bool globcover,bool overwrite,const geoExtend &myExtend);
    void TestConsistency();
    void Map(const std::string &str_ifile,int reduction_factor,const geoExtend &myExtend);
    int OpenInputRaster(const std::string str_ifile);
  private:
    int ReadLabels(const string &strfile,int64_t &maxlabel,vector <tlabel>&labels);
    void TransferLabels(const vector <tlabel>&labels,int64_t maxlabel,int ltype,vector <double>&reflabels);
    IMGPGM myPGM;
    IMGBRI myBRI;
    IMGASC myASC;
    IMGTIFF myTIFF;
    IMG *myIMG;
    FILE *file;
};

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
        cout << myTimer.SecElapsed() << endl;
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
      ClusterBRI myCluster(options);
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
      cout << endl << "time: " << Utils::SecToTime(myTimer.SecElapsed()) << endl;

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

int ComandLine::ReadLabels(const string &strfile,int64_t &maxlabel,vector <tlabel>&labels)
{
  maxlabel=0;
  ifstream ifile(strfile);
  if (ifile.is_open()) {
    int64_t flength=Utils::GetStreamSize(ifile);
    double total_area=0.;
    double total_edgelen=0.;
    std::string line;
    int cnt=0;
    while ( getline (ifile,line) )
    {
      if ( (++cnt)%10000==0) {
        float fratio=ifile.tellg()*100.0/(double)flength;
        cout << "reading labels: '" << strfile << "': ";
        cout << Utils::ConvertFixed(fratio,1) << "%\r";
        cnt=0;
      }
      std::vector <std::string> tokens;
      StringUtils::Tokenize(line,tokens,", \t\n");
      if (tokens.size()!=5) cout << "warning: undefined number of tokens: " << tokens.size() << "\n";
      else {
        tlabel labelentry;
        labelentry.label=strtoll(tokens[0].c_str(),NULL,10);
        labelentry.area=atof(tokens[2].c_str());
        labelentry.edge_len=atof(tokens[3].c_str());
        labelentry.closs=atof(tokens[4].c_str());

        total_area+=labelentry.area;
        total_edgelen+=labelentry.edge_len;
        if (labelentry.label>maxlabel) maxlabel=labelentry.label;
        labels.push_back(labelentry);
      }
    }
    ifile.close();
    cout << endl;
    cout << "number of labels: " << labels.size() << "\n";
    cout << "total area:       " << Utils::ConvertFixed(Utils::SqMetre_To_MillHa(total_area),2) << " 10^6 ha" << endl;
    cout << "edge len:         " << Utils::ConvertFixed(Utils::Metre_To_MillKm(total_edgelen),2) << " 10^6 km" << endl;
    return 0;
  } else return 1;
}

void ComandLine::TransferLabels(const vector <tlabel>&labels,int64_t maxlabel,int ltype,vector <double>&reflabels)
{
  cout << "transfering labels: ";
  switch (ltype) {
    case 0: cout << "closs";break;
    case 1: cout << "core/area";break;
    default:cout << "unknown";break;
  }
  cout << endl;
  reflabels.resize(maxlabel+1);
  for (size_t i=0;i<labels.size();i++) {
    double val=0.;
    if (ltype==0) val=labels[i].closs;
    else if (ltype==1) {
      double edge_area=ClusterBRI::CalculateEdgeAreaDE(labels[i].area,labels[i].edge_len,100);
      double core_area=labels[i].area-edge_area;
      val=core_area/labels[i].area;
    }
    reflabels[labels[i].label]=val;
  }
}

void ComandLine::Map(const std::string &str_ifile,int reduction_factor,const geoExtend &myExtend)
{
  std::string str_lfile,str_rfile,str_pfile;
  Utils::ReplaceExt(str_ifile,str_lfile,".lab",true);
  Utils::ReplaceExt(str_ifile,str_rfile,".dat",true);
  Utils::ReplaceExt(str_ifile,str_pfile,".txt",true);

  std::vector <tlabel> labels;
  std::vector <double> reflabels;
  int64_t maxlabel;
  ReadLabels(str_lfile,maxlabel,labels);
  TransferLabels(labels,maxlabel,1,reflabels);
  labels.clear();

  cout << "reading clusters from '" << str_ifile << "': ";
  FILE *clusterfile=fopen(str_ifile.c_str(),"rb");
  if (clusterfile!=NULL) {
    uint8_t tbuf[8];
    fread(tbuf,1,8,clusterfile);
    uint32_t width=Utils::Get32LH(tbuf);
    uint32_t height=Utils::Get32LH(tbuf+4);
    cout << width << "x" << height << endl;

    Projection Proj(width,height);
    Proj.ReadCoordinateFile(str_pfile);
    Proj.CalculateCellSize();

    int pleft,pright,ptop,pbottom;
    Frame::SetExtend(Proj.getLeft(),Proj.getTop(),Proj.getCellsize(),myExtend,width,height,pleft,ptop,pright,pbottom);

    uint8_t *rowdata=new uint8_t[width*8];
    int64_t *labelrow=new int64_t[width];
    size_t tread;
    //vector <vector<double>>datarows(reduction_factor,vector<double>(width));
    double **datarows=new double*[reduction_factor];
    for (int i=0;i<reduction_factor;i++) datarows[i]=new double[width];

    cout << "reduction factor: 1:" << reduction_factor << ", array: " << (((int64_t)reduction_factor*(int64_t)width*sizeof(double))>>20) << " mb\n";

    ofstream rfile(str_rfile);

    if (ptop) {
      for (int i=0;i<ptop-1;i++) {
        tread=fread(tbuf,1,4,clusterfile);
        uint32_t rowsize=Utils::Get32LH(tbuf);
        fseek(clusterfile,rowsize,SEEK_CUR);
        if (i%10==0) cout << "skipping " << (ptop-1) << " lines: " << Utils::ConvertFixed(i*100/(double)(ptop-1),1) << "%\r";
      }
      cout << endl;
    }

    int64_t total_cells=0;
    int rowcnt=0;
    for (int row=ptop;row<pbottom;row++) {
        tread=fread(tbuf,1,4,clusterfile);
        uint32_t rowsize=Utils::Get32LH(tbuf);
        tread=fread(rowdata,1,rowsize,clusterfile); // possible buffer overflow if rowsize>8*width
        if (tread!=rowsize) cout << "warning: could not read\n";
        else {
          RLEPack2::UnpackRow(rowdata,width,labelrow);
          double *drow=datarows[rowcnt];
          int j=0;
          for (int i=pleft;i<pright;i++) {
            drow[j]=0.;
            if (labelrow[i]) {
              total_cells++;
              if (labelrow[i]>maxlabel) cout << "warning: label outside range\n";
              else drow[j]=reflabels[labelrow[i]];
            };
            j++;
          }
        }
        rowcnt++;
        cout << (row+1) << "/" << height <<"\r";
        if (rowcnt>=reduction_factor) {
          RowReduce(rfile,datarows,pright-pleft,reduction_factor,rowcnt);
          rowcnt=0;
        }
    }
    if (rowcnt) RowReduce(rfile,datarows,pright-pleft,reduction_factor,rowcnt);
    cout << "total cells:   " << total_cells << endl;
    rfile.close();
    for (int i=0;i<reduction_factor;i++) delete []datarows[i];
    delete []datarows;
    delete []rowdata;
    delete []labelrow;
    fclose(clusterfile);
  } else cout << "error: could not open: '" << str_ifile << "'\n";
}

void TestBitIO()
{
    srand(time(0));
    uint8_t buf[256];
    BitBuffer bitout(buf),bitin(buf);
    for (int i=0;i<17;i++) {
      int r=rand() & 1;
      cout << r;
      bitout.PutBit(r);
    }
    int val;
    val=113;bitout.PutBits(val,7);cout << "," << val;
    val=0;bitout.PutRice(val,2);cout << "," << val;
    val=(1<<23);bitout.PutRice(val,31);cout << "," << val;
    bitout.Flush();
    cout << ":" << bitout.GetBytesProcessed() << endl;
    for (int i=0;i<17;i++) {
      int bit;
      bit = bitin.GetBit();
    cout << bit;
    }
    val = bitin.GetBits(7);cout << "," << val;
    val = bitin.GetRice(2);cout << "," << val;
    val = bitin.GetRice(31);cout << "," << val;
    cout << endl;
}

const std::string LISA_USAGE={
"lisa [-/--options]\n\n"
"-a,--analyze  analyze connected components of a raster file (#-area cells)\n"
"-c,--convert  convert raster file into [bri] file (g=globcover)\n"
"-d,--dept     edge effect dept of d [m], default: 100\n"
"-e,--extend   top,left,right,bottom\n"
"-f,--fragment statistics using minimum fragment size in [ha] (default: 0)\n"
"-m,--map      produce a density map, out of bin/lab file\n"
"-p,--pixel    minimum pixel length for edge detection, default: 1\n"
"-s,--save     save results to .csv file, 1=small, 2=large\n"
"-t,--test     test consistence of lisa\n"
"-v,--verbose  verbosity level [0-2] (default: 1)\n"
"-w,--write    write clusterlabel data, 1=clusters+labels, 2=labels\n"
"--info        info about raster file\n"
"--input       inputfile\n"
"--output      outputfile\n"
"--version     print version info\n"
"--nrows       number of rows to process\n"
"--threshold   forest cover threshold for forest/nonforest map\n"
"--flush       flush clusters to use fixed amount of memory\n"
"--agb-file    saatchi agb biomass file [t/ha]\n"
"--check       check consistency of component analysis\n"
"--bthres      biomass threshold [t/ha] (default: 0 t/ha)\n"
"--rloss       relative carbon loss in edge areas, default: 0.5\n"
"--force       force overwrite of files\n"
"-r[#]   reduction factor for use with --map option, default: 500\n"
"-b[#]   mean biomass for use with --analyze\n"
"supported raster file formats: asc, pgm, tiff, bri\n"
};

int main(int argc,char *argv[])
{
    #if 0
    srand(time(0));
    uint8_t buf[1024];
    BitBuffer64LH bout(buf),bin(buf);
    int num=32;
    for (int i=0;i<num;i++) {
       int b=(rand()%15)+1;
       bout.PutEliasGamma(b);
       cout << b << " ";
    }
    bout.Flush();
    cout << endl << bout.GetBytesProcessed() << endl;
    for (int i=0;i<num;i++) {
       int b=bin.GetEliasGamma();
       cout << b << " ";
    }
    cout << endl << bin.GetBytesProcessed() << endl;
    return 0;
    #endif
    AnalyzeOptions AnalyzeOptions;

    int verbosity_level=1;
    ComandLine::METHOD cmode=ComandLine::ANALYZE;

    bool globcover=false;
    int reduction_factor=500;
    bool force_overwrite=false;
    std::string str_ifile,str_ofile;
    std::string str_bfile;
    geoExtend myGeoExtend;

    if (argc < 2) {
      cout << LISA_USAGE << endl;
      return 1;
    }

    CmdOptions myCmdOpt(argc,argv);

    if (myCmdOpt.SearchOption("-e","--extend")) {
      std::vector <double>vextend;
      myCmdOpt.getopt(vextend);
      if (vextend.size()!=4) cout << "warning: unexpected size of extend=" << vextend.size() << endl;
      else {
        myGeoExtend.top=vextend[0];
        myGeoExtend.left=vextend[1];
        myGeoExtend.right=vextend[2];
        myGeoExtend.bottom=vextend[3];
      }
    }
    if (myCmdOpt.SearchOption("-a","--analyze")) cmode=ComandLine::ANALYZE;
    if (myCmdOpt.SearchOption("","--info")) cmode=ComandLine::INFO;
    if (myCmdOpt.SearchOption("-c","--convert")) cmode=ComandLine::CONVERT;
    if (myCmdOpt.SearchOption("-t","--test")) cmode=ComandLine::TEST;
    if (myCmdOpt.SearchOption("-m","--map")) cmode=ComandLine::MAP;
    if (myCmdOpt.SearchOption("","--version")) cmode=ComandLine::VERSION;
    if (myCmdOpt.SearchOption("","--force")) force_overwrite=true;
    if (myCmdOpt.SearchOption("-v","--verbose")) myCmdOpt.getopt(verbosity_level);
    if (myCmdOpt.SearchOption("-d","--dept")) {myCmdOpt.getopt(AnalyzeOptions.edge_dept);};
    if (myCmdOpt.SearchOption("","--threshold")) {myCmdOpt.getopt(AnalyzeOptions.forest_cover_threshold);};
    if (myCmdOpt.SearchOption("","--nrows")) {myCmdOpt.getopt(AnalyzeOptions.nrows);};
    if (myCmdOpt.SearchOption("","--flush")) {AnalyzeOptions.flush_clusters=true;};
    if (myCmdOpt.SearchOption("-f","--fragment")) {myCmdOpt.getopt(AnalyzeOptions.min_fragment_size);};
    if (myCmdOpt.SearchOption("-w","--write")) myCmdOpt.getopt(AnalyzeOptions.write_mode);
    if (myCmdOpt.SearchOption("-s","--save")) myCmdOpt.getopt(AnalyzeOptions.save_mode);
    if (myCmdOpt.SearchOption("-b","")) myCmdOpt.getopt(AnalyzeOptions.mean_biomass);
    if (myCmdOpt.SearchOption("-r","")) myCmdOpt.getopt(reduction_factor);
    if (myCmdOpt.SearchOption("-p","--pixel")) myCmdOpt.getopt(AnalyzeOptions.pixel_len);
    if (myCmdOpt.SearchOption("-i","--input")) myCmdOpt.getopt(str_ifile);
    if (myCmdOpt.SearchOption("-o","--output")) myCmdOpt.getopt(str_ofile);
    if (myCmdOpt.SearchOption("","--agb-file")) myCmdOpt.getopt(str_bfile);
    if (myCmdOpt.SearchOption("","--bthres")) myCmdOpt.getopt(AnalyzeOptions.bthres);
    if (myCmdOpt.SearchOption("","--rloss")) myCmdOpt.getopt(AnalyzeOptions.relative_carbon_loss);
    if (myCmdOpt.SearchOption("","--check")) AnalyzeOptions.check_consistency=true;

    if (verbosity_level>1) {
       cout << "mode:          " << cmode << endl;
       cout << "edge dept:     " << AnalyzeOptions.edge_dept << endl;
       cout << "mean biosmass: " << AnalyzeOptions.mean_biomass << endl;
       cout << "min fragment:  " << AnalyzeOptions.min_fragment_size << endl;
       cout << "pixel len:     " << AnalyzeOptions.pixel_len << endl;
       cout << "rel. c-loss:   " << AnalyzeOptions.relative_carbon_loss << endl;
       cout << "savemode:      " << AnalyzeOptions.save_mode << endl;
       cout << "writemode:     " << AnalyzeOptions.write_mode << endl;
       cout << "bthres:        " << AnalyzeOptions.bthres << endl;
       cout << "extend:        " << myGeoExtend.top<<","<<myGeoExtend.left<<","<<myGeoExtend.right<<","<<myGeoExtend.bottom<<endl;
       cout << "infile:        '" << str_ifile << "'" << endl;
       cout << "outfile:       '" << str_ofile << "'" << endl;
       cout << "agb-file:      '" << str_bfile << "'" << endl;
    }

    ComandLine myCmdLine;
    switch (cmode) {
      case ComandLine::ANALYZE:myCmdLine.Analyze(str_ifile,str_bfile,AnalyzeOptions);break;
      case ComandLine::CONVERT:
      case ComandLine::INFO:myCmdLine.Convert(str_ifile,str_ofile,cmode,globcover,force_overwrite,myGeoExtend);break;
      case ComandLine::TEST:myCmdLine.TestConsistency();break;
      case ComandLine::MAP:myCmdLine.Map(str_ifile,reduction_factor,myGeoExtend);break;
      case ComandLine::VERSION:PrintVersion(1);break;
      default: cout << "unknown mode: " << cmode << endl;break;
    }
    return 0;
}
