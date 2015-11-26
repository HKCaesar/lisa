#include "global.h"
#include "rand.h"
#include "file/pgm.h"
#include "file/bri.h"
#include "file/asc.h"
#include "file/bm.h"
#include "bmatrix.h"
#include "analysis/cluster.h"
#include "analysis/proj.h"
#include "model/rangecoder.h"
#include "model/counter.h"

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

  BRI myBRI;
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

      BRIOptions options(myBRI,myProj,myBiomass);

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

int StringToInt(std::string &s)
{
  int val=0;
  try
  {
    val=std::stoi(s);
  }
  catch(std::invalid_argument&) //or catch(...) to catch all exceptions
  {
    std::cout << "invalid argument: '" << s << "'" << std::endl;
  }
  return val;
}

double StringToDouble(std::string &s)
{
  double val=0;
  try
  {
    val=std::stod(s);
  }
  catch(std::invalid_argument&) //or catch(...) to catch all exceptions
  {
    std::cout << "invalid argument: '" << s << "'" << std::endl;
  }
  return val;
}

bool isDouble( string myString ) {
    std::istringstream iss(myString);
    double f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail();
}

void SplitPath(const std::string &str,std::string &path,std::string &fname)
{
  unsigned found = str.find_last_of("/\\");
  path=str.substr(0,found);
  fname=str.substr(found+1);
}

double BlockReduce(float **rows,int width,int startx,int xblock,int yblock)
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

void RowReduce(ofstream &stream,float **rows,int width,int xblock,int yblock)
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
}

class CmdOptions {
  public:
    CmdOptions(int argc,char *argv[]):argc(argc),argv(argv) {};
    bool SearchOption(const std::string &sshort,const std::string &&sslong); // check if option exists
    void getopt(std::string &val);
    void getopt(int &val);
    void getopt(double &val);
    void getopt(std::vector <double>&dtokens);
private:
    int argc;
    char **argv;
    bool optfound;
    std::string optvalue;
};

void CmdOptions::getopt(std::string &val)
{
  if (optfound && optvalue.length()) {
    val=optvalue;
  };
}

void CmdOptions::getopt(int &val)
{
  if (optfound && optvalue.length()) {
     try {
       val=std::stoi(optvalue);
     } catch (std::invalid_argument&) {
     }
  };
}

void CmdOptions::getopt(double &val)
{
  if (optfound && optvalue.length()) {
     try {
       val=std::stod(optvalue);
     } catch (std::invalid_argument&) {
     }
  };
}

void CmdOptions::getopt(std::vector <double>&dtokens)
{
  std::vector <std::string> stokens;
  StringUtils::Tokenize(optvalue,stokens,",");
  for (size_t i=0;i<stokens.size();i++) {
      double val=std::stod(stokens[i]);
      dtokens.push_back(val);
  }
}

bool CmdOptions::SearchOption(const std::string &sshort,const std::string &&sslong)
{
  std::string usshort=StringUtils::toupper(sshort);
  std::string usslong=StringUtils::toupper(sslong);
  optvalue="";
  optfound=false;
  int i=1;
  while (i<argc) {
    std::string sarg=StringUtils::toupper(std::string(argv[i]));
    std::string arglong(sarg);
    std::string argshort(sarg.substr(0,2));

    if (usshort.length() && argshort.compare(usshort)==0) {
      if (sarg.length()>2) optvalue=sarg.substr(2);
      optfound=true;
      return optfound;
    } else if (usslong.length() && arglong.compare(usslong)==0) {
      if (i<argc-1) optvalue=std::string(argv[i+1]);
      optfound=true;
      return optfound;
    }
    i++;
  }
  return optfound;
}

struct tAnalyzeOptions {
  double mean_biomass,bthres,rloss;
  int ncells,edge_dept,min_fragment,pixel_len,write_mode,save_mode;
  bool check_consistency;
};

class ComandLine {
  struct tlabel {
    int64_t label;
    float closs;
  };
  public:
      static void Analyze(const std::string &str_ifile,const std::string &str_bfile,const tAnalyzeOptions &AnalyzeOptions);
      static void Convert(const std::string &str_ifile,std::string &str_ofile,int cmode,bool globcover,bool overwrite,const geoExtend &myExtend);
      static void TestConsistency();
      static void Reduce(const std::string &str_ifile,int reduction_factor,const geoExtend &myExtend);
};

void ComandLine::Analyze(const std::string &str_ifile,const std::string &str_bfile,const tAnalyzeOptions &AnalyzeOptions)
{
  FILE *file=fopen(str_ifile.c_str(),"rb");
  cout << "open file: '" << str_ifile << "': ";
  if (file)  {
        cout << "ok." << endl;
        BRI myBRI;
        myBRI.SetHandle(file);
        if (myBRI.ReadHeader()==0)
        {
          std::string str_pfile;
          Utils::ReplaceExt(str_ifile,str_pfile,".txt",true);

          Projection myProj(myBRI.GetWidth(),myBRI.GetHeight());
          myBRI.PrintInfo();
          if (myProj.ReadCoordinateFile(str_pfile)==0) {
            cout << endl;
            Timer myTimer;
            myTimer.Start();
            myProj.CalculateCellSize();
            myProj.GenerateInterpolation(AnalyzeOptions.ncells);

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

            BRIOptions options(myBRI,myProj,myBiomass);

            options.edge_dept=AnalyzeOptions.edge_dept;
            options.relative_carbon_loss=AnalyzeOptions.rloss;
            options.min_fragment_size=AnalyzeOptions.min_fragment;
            options.pixel_len=AnalyzeOptions.pixel_len;
            options.write_clusterlabel=AnalyzeOptions.write_mode;
            options.verbose=true;
            if (options.write_clusterlabel>0) {
                Utils::ReplaceExt(str_ifile,options.str_clusterfile1,".tmp",true);
                Utils::ReplaceExt(str_ifile,options.str_clusterfile2,".bin",true);
                Utils::ReplaceExt(str_ifile,options.str_labelfile,".lab",true);
            }
            ClusterBRI myCluster(options);
            myCluster.ClusterAnalyzation();

            if (AnalyzeOptions.check_consistency) myCluster.CheckClusters();

            std::string str_ofile;
            Utils::ReplaceExt(str_ifile,str_ofile,".csv",true);
            if (AnalyzeOptions.save_mode>0) {
              cout << "Saving clusters to '" << str_ofile << "'" << endl;
              if (AnalyzeOptions.save_mode==1) myCluster.SaveSmallClusterData(str_ofile);
              else if (AnalyzeOptions.save_mode==2) myCluster.SaveFullClusterData(str_ofile);
              else cout << "unknown save_mode: " << AnalyzeOptions.save_mode << endl;
            }
            myTimer.Stop();
            cout << endl << "time: " << Utils::SecToTime(myTimer.SecElapsed()) << endl;
          } else cout << "missing coordinate file: '" << str_pfile << "'!" << endl;

        } else cout << "no valid file type found!" << endl;
        fclose(file);
      } else cout << "not found!";
}

void ComandLine::Convert(const std::string &str_ifile,std::string &str_ofile,int cmode,bool globcover,bool force_overwrite,const geoExtend &myExtend)
{
      FILE *file=fopen(str_ifile.c_str(),"rb");
      cout << "open file: '" << str_ifile << "': ";
      if (file)
      {
         cout << "ok." << endl;
         PGM myPGM;
         BRI myBRI;
         ASC myASC;
         myPGM.SetHandle(file);
         if (myPGM.ReadHeader()==0)
         {
           myPGM.PrintInfo();
           if (cmode==2) { // convert pgm to bri
            Utils::ReplaceExt(str_ifile,str_ofile,".bri");
            cout << "converting to '" << str_ofile << "'" << endl;
            myBRI.ConvertFromPGM(myPGM,str_ofile,globcover);
           }
         } else {
           myBRI.SetHandle(file);
           if (myBRI.ReadHeader()==0)
           {
             myBRI.PrintInfo();
             if (cmode==2) { // convert bri to pgm
                Utils::ReplaceExt(str_ifile,str_ofile,".bri");
                cout << "converting to '" << str_ofile << "'" << endl;
                myBRI.ConvertToPGM(myPGM,str_ofile);
             }
           } else {
               myASC.SetHandle(file);
               if (myASC.ReadHeader()==0) {
                  myASC.PrintInfo();
                  if (cmode==2) {
                    myASC.SetExtend(myExtend);
                    std::string str_efile;
                    Utils::ReplaceExt(str_ofile,str_efile,".txt");
                    cout << "saving geo-extend to '" << str_efile << "'\n";
                    if (Utils::OpenWriteCheck(str_efile,force_overwrite)) myASC.WriteExtend(str_efile,8);
                    else return;

                    Utils::ReplaceExt(str_ifile,str_ofile,".bri");
                    cout << "converting to '" << str_ofile << "'\n";
                    myBRI.ConvertFrom(myASC,str_ofile);
                  }
               } else cout << "error: unknown file format\n";
           }
         }
         fclose(file);
      } else cout << "not found!";
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

void ComandLine::Reduce(const std::string &str_ifile,int reduction_factor,const geoExtend &myExtend)
{
      cout << str_ifile << endl;
      std::vector <tlabel> vlabel;
      std::vector <float> reflabel;
      std::string str_lfile,str_rfile,str_pfile;
      Utils::ReplaceExt(str_ifile,str_lfile,".lab",true);
      Utils::ReplaceExt(str_ifile,str_rfile,".dat",true);
      Utils::ReplaceExt(str_ifile,str_pfile,".txt",true);


      ifstream lfile(str_lfile);
      int64_t maxlabel=0;
      int cnt=0;
      if (lfile.is_open()) {
        int64_t flength=Utils::GetStreamSize(lfile);
        std::string line;
        while ( getline (lfile,line) )
        {
          cnt++;
          if (cnt%10000==0) {
            float fratio=lfile.tellg()*100.0/(double)flength;
            cout << "reading labels: '" << str_lfile << "': ";
            cout << Utils::ConvertFixed(fratio,1) << "%\r";
            cnt=0;
          }
          std::vector <std::string> tokens;
          StringUtils::Tokenize(line,tokens,", \t\n");
          if (tokens.size()!=2) cout << "warning: undefined number of tokens: " << tokens.size() << "\n";
          else {
            tlabel labelentry;
            labelentry.label=strtoll(tokens[0].c_str(),NULL,10);
            labelentry.closs=atof(tokens[1].c_str());
            if (labelentry.label>maxlabel) maxlabel=labelentry.label;
            vlabel.push_back(labelentry);
          }
        }
        lfile.close();
        cout << "\nok. number of labels: " << vlabel.size() << "\n";

        cout << "transfering labels...";
        reflabel.resize(maxlabel+1);
        for (size_t i=0;i<vlabel.size();i++) {
            reflabel[vlabel[i].label]=vlabel[i].closs;
        }
        cout << "ok\n";
        vlabel.resize(0);

        cout << "reading clusters from '" << str_ifile << "'...";
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
            cout << Proj.getCellsize() << endl;

            int pleft,pright,ptop,pbottom;
            Frame::SetExtend(Proj.getLeft(),Proj.getTop(),Proj.getCellsize(),myExtend,width,height,pleft,ptop,pright,pbottom);

            uint8_t *rowdata=new uint8_t[width*8];
            int64_t *labelrow=new int64_t[width];
            size_t tread;
            float **datarows=new float*[reduction_factor];
            for (int i=0;i<reduction_factor;i++) datarows[i]=new float[width];
            cout << "reduction factor: 1:" << reduction_factor << ", array: " << (((int64_t)reduction_factor*(int64_t)width*sizeof(float))>>20) << " mb\n";

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
                //cout << rowsize << endl;
                tread=fread(rowdata,1,rowsize,clusterfile); // possible buffer overflow if rowsize>8*width
                if (tread!=rowsize) cout << "warning: could not read\n";
                else {
                  RLEPack::UnpackRow(rowdata,width,labelrow);
                  float *drow=datarows[rowcnt];
                  int j=0;
                  for (int i=pleft;i<pright;i++) {
                    drow[j]=0.;
                    if (labelrow[i]) {
                        total_cells++;
                        if (labelrow[i]>maxlabel) cout << "warning: label outside range\n";
                        else drow[j]=reflabel[labelrow[i]];
                    };
                    j++;
                  }
                  rowcnt++;
                  cout << (row+1) << "/" << height <<"\r";
                  if (rowcnt>=reduction_factor) {
                    /*int w=pright-pleft;
                    for (int i=0;i<w;i++) {
                        rfile << std::to_string(drow[i]);
                        if (i<w-1) rfile<<",";
                    }
                    rfile << endl;*/
                    RowReduce(rfile,datarows,pright-pleft,reduction_factor,rowcnt);
                    rowcnt=0;
                  }
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
      } else cout << "error: could not open: '" << str_lfile << "'\n";

}


int main(int argc,char *argv[])
{
    tAnalyzeOptions AnalyzeOptions;
    AnalyzeOptions.bthres=0.; // biomass threshold: 0 t/ha
    AnalyzeOptions.edge_dept=100; // edge effect dept 100m
    AnalyzeOptions.mean_biomass=0.;
    AnalyzeOptions.min_fragment=0; // minimum fragment size 0 ha
    AnalyzeOptions.ncells=1; // number of interpolation cells
    AnalyzeOptions.pixel_len=1; // pixel number for edge detection
    AnalyzeOptions.rloss=0.5; // relative carbon loss in edge areas
    AnalyzeOptions.save_mode=0;
    AnalyzeOptions.write_mode=0;
    AnalyzeOptions.check_consistency=false;

    int verbosity_level=1;
    int cmode=0;

    bool globcover=false;
    int reduction_factor=500;
    bool force_overwrite=false;
    std::string str_ifile,str_ofile;
    std::string str_bfile;
    geoExtend myGeoExtend;

    if (argc < 2)
    {
       cout << "lisa [-/--options]" << endl << endl;
       cout << "-a,--analyze  analyze connected components of a bri file (#-area cells)" << endl;
       cout << "-c,--convert  convert [pgm/bri/asc] file into [bri] file (g=globcover)" << endl;
       cout << "-d,--dept     edge effect dept of d [m], default: 100" << endl;
       cout << "-e,--extend   top,left,right,bottom" << endl;
       cout << "-f,--fragment statistics using minimum fragment size in [ha] (default: 0)" << endl;
       cout << "-i,--info     info about [pgm/bri/asc] file" << endl;
       cout << "-m,--map      produce a density map, out of bin/lab file" << endl;
       cout << "-p,--pixel    minimum pixel length for edge detection, default: 1" << endl;
       cout << "-s,--save     save results to .csv file, 1=small, 2=large" << endl;
       cout << "-t,--test     test consistence of lisa" << endl;
       cout << "-v,--verbose  verbosity level [0-2] (default: 1)" << endl;
       cout << "-w,--write    write clusterlabel data, 1=labels+closs, 2=closs" << endl;
       cout << "--input       inputfile" << endl;
       cout << "--output      outputfile" << endl;
       cout << "--version     print version info" << endl;
       cout << "--agb-file    saatchi agb biomass file [t/ha]" << endl;
       cout << "--check       check consistency of component analysis" << endl;
       cout << "--bthres      biomass threshold [t/ha] (default: 0 t/ha)"<<endl;
       cout << "--rloss       relative carbon loss in edge areas, default: 0.5" << endl;
       cout << "--force       force overwrite of files"<<endl;
       cout << "-r[#]   reduction factor for use with --map option, default: 500"<<endl;
       cout << "-b[#]   mean biomass for use with --analyze"<<endl;
       return 1;
    } else {
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
    if (myCmdOpt.SearchOption("-a","--analyze")) {cmode=0;myCmdOpt.getopt(AnalyzeOptions.ncells);};
    if (myCmdOpt.SearchOption("","--info")) cmode=1;
    if (myCmdOpt.SearchOption("-c","--convert")) cmode=2;
    if (myCmdOpt.SearchOption("-t","--test")) cmode=3;
    if (myCmdOpt.SearchOption("-m","--map")) cmode=4;
    if (myCmdOpt.SearchOption("","--version")) cmode=5;
    if (myCmdOpt.SearchOption("","--force")) force_overwrite=true;
    if (myCmdOpt.SearchOption("-v","--verbose")) myCmdOpt.getopt(verbosity_level);
    if (myCmdOpt.SearchOption("-d","--dept")) {myCmdOpt.getopt(AnalyzeOptions.edge_dept);};
    if (myCmdOpt.SearchOption("-f","--fragment")) {myCmdOpt.getopt(AnalyzeOptions.min_fragment);};
    if (myCmdOpt.SearchOption("-w","--write")) myCmdOpt.getopt(AnalyzeOptions.write_mode);
    if (myCmdOpt.SearchOption("-s","--save")) myCmdOpt.getopt(AnalyzeOptions.save_mode);
    if (myCmdOpt.SearchOption("-b","")) myCmdOpt.getopt(AnalyzeOptions.mean_biomass);
    if (myCmdOpt.SearchOption("-r","")) myCmdOpt.getopt(reduction_factor);
    if (myCmdOpt.SearchOption("-p","--pixel")) myCmdOpt.getopt(AnalyzeOptions.pixel_len);
    if (myCmdOpt.SearchOption("-i","--input")) myCmdOpt.getopt(str_ifile);
    if (myCmdOpt.SearchOption("-o","--output")) myCmdOpt.getopt(str_ofile);
    if (myCmdOpt.SearchOption("","--agb-file")) myCmdOpt.getopt(str_bfile);
    if (myCmdOpt.SearchOption("","--bthres")) myCmdOpt.getopt(AnalyzeOptions.bthres);
    if (myCmdOpt.SearchOption("","--rloss")) myCmdOpt.getopt(AnalyzeOptions.rloss);
    if (myCmdOpt.SearchOption("","--check")) AnalyzeOptions.check_consistency=true;

    if (verbosity_level>1) {
       cout << "mode:          " << cmode << endl;
       cout << "ncells:        " << AnalyzeOptions.ncells << endl;
       cout << "edge dept:     " << AnalyzeOptions.edge_dept << endl;
       cout << "mean biosmass: " << AnalyzeOptions.mean_biomass << endl;
       cout << "min fragment:  " << AnalyzeOptions.min_fragment << endl;
       cout << "pixel len:     " << AnalyzeOptions.pixel_len << endl;
       cout << "rel. c-loss:   " << AnalyzeOptions.rloss << endl;
       cout << "savemode:      " << AnalyzeOptions.save_mode << endl;
       cout << "writemode:     " << AnalyzeOptions.write_mode << endl;
       cout << "bthres:        " << AnalyzeOptions.bthres << endl;
       cout << "extend:        " << myGeoExtend.top<<","<<myGeoExtend.left<<","<<myGeoExtend.right<<","<<myGeoExtend.bottom<<endl;
       cout << "infile:        '" << str_ifile << "'" << endl;
       cout << "outfile:       '" << str_ofile << "'" << endl;
       cout << "agb-file:      '" << str_bfile << "'" << endl;
    }

    }

    if (cmode==0)
    {
      ComandLine::Analyze(str_ifile,str_bfile,AnalyzeOptions);
    } else if (cmode==1 || cmode==2)
    {
      ComandLine::Convert(str_ifile,str_ofile,cmode,globcover,force_overwrite,myGeoExtend);
    } else if (cmode==3)
    {
      ComandLine::TestConsistency();
    } else if (cmode==4)
    {
      ComandLine::Reduce(str_ifile,reduction_factor,myGeoExtend);
    } else if (cmode==5) {
      PrintVersion(1);
    }

    //m.Print();
    return 0;
}
