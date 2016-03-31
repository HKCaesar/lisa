#include "map.h"
#include "../model/sic.h"
#include "../common/timer.h"

// wrapper for strtok, str gets destroyed
class StrTok {
  public:
    StrTok()
    {
    }
    StrTok(std::string &str,const char *delim)
    {
      FirstToken(str,delim);
    }
    // requires the iternal string buffer to end on '\0'
    void FirstToken(std::string &str,const char *delim)
    {
      delim_=delim;
      token_=strtok(&str[0],delim_);
    }
    void NextToken()
    {
      token_=strtok(NULL,delim_);
    }
    bool valid(){return token_!=NULL;};
    double GetDouble(){return atof(token_);};
    long long int GetLongLong(){return strtoll(token_,NULL,10);};
  private:
    const char *delim_;
    char *token_;
};

int Map::ReadLabels(const string &strfile)
{
  maxlabel_=0;
  labels_.clear();

  ifstream ifile(strfile);
  if (ifile.is_open()) {
    int64_t filelength=Utils::GetStreamSize(ifile);
    double total_area=0.;
    double total_edgelen=0.;

    Timer myTimer;
    myTimer.Start();

    MapLabel labelentry;
    std::string line;
    int cnt=0;
    const char delim[]=", \t\n";
    //std::ios::sync_with_stdio(false);
    while ( getline (ifile,line) )
    {
      if ( (++cnt)==1<<16) {
        cout << "reading labels: '" << strfile << "': " << Utils::ConvertFixed(ifile.tellg()*100.0/(double)filelength,1) << "%\r";
        cnt=0;
      }

      #if 0
        std::vector <std::string> tokens;
        StringUtils::Tokenize(line,tokens,delim);

        if (tokens.size()!=5) cout << "warning: undefined number of tokens: " << tokens.size() << "\n";
        else {
          labelentry.label=strtoll(tokens[0].c_str(),NULL,10);
          labelentry.area=atof(tokens[2].c_str());
          labelentry.edge_len=atof(tokens[3].c_str());
          labelentry.closs=atof(tokens[4].c_str());
        }
      #else // twice as fast
        StrTok tok(line,delim);
        labelentry.label=tok.GetLongLong();
        tok.NextToken(); // skip area in pixels
        tok.NextToken();
        labelentry.area=tok.GetDouble();
        tok.NextToken();
        labelentry.edgelen=tok.GetDouble();
        tok.NextToken();
        labelentry.closs=tok.GetDouble();
        if (!tok.valid()) cout << "map: unexpected number of tokens in file: '" << strfile << "'\n";
      #endif

      total_area+=labelentry.area;
      total_edgelen+=labelentry.edgelen;
      if (labelentry.label>maxlabel_) maxlabel_=labelentry.label;
      labels_.push_back(labelentry);
    }
    ifile.clear();
    cout << "reading labels: '" << strfile << "': " << Utils::ConvertFixed(ifile.tellg()*100.0/(double)filelength,1) << "%\n";
    ifile.close();
    myTimer.Stop();
    cout << "number of labels: " << labels_.size() << "\n";
    cout << "total area:       " << Utils::ConvertFixed(Utils::SqMetre_To_MillHa(total_area),2) << " 10^6 ha" << endl;
    cout << "edge len:         " << Utils::ConvertFixed(Utils::Metre_To_MillKm(total_edgelen),2) << " 10^6 km" << endl;
    cout << "time elapsed:     " << myTimer.ElapsedS() << " seconds" << endl;
    return 0;
  } else return 1;
}

void Map::TransferLabels()
{
  reflabels_.resize(maxlabel_+1);
  for (size_t i=0,ilen=labels_.size();i<ilen;i++) {
    double val=0.;
    if (maptype_==0) val=labels_[i].closs;
    else if (maptype_==1) {
      double edge_area=Cluster::CalculateEdgeAreaDE(labels_[i].area,labels_[i].edgelen,edge_effect_dept_);
      double core_area=labels_[i].area-edge_area;
      val=core_area/labels_[i].area;
    }
    reflabels_[labels_[i].label]=val;
  }
  labels_.clear();
}

void Map::AllocMem()
{
  labelrow_=new int64_t[width_];
  datarows_=new double*[reduction_factor_];
  for (int i=0;i<reduction_factor_;i++) datarows_[i]=new double[width_];
}

void Map::FreeMem()
{
  for (int i=0;i<reduction_factor_;i++) delete []datarows_[i];
  delete []datarows_;
  delete []labelrow_;
}

void Map::SkipRows(int nrows)
{
  uint8_t tbuf[4];
  for (int i=0;i<nrows;i++) {
    size_t tread=fread(tbuf,1,4,clusterfile_);
    uint32_t rowsize=Utils::Get32LH(tbuf);
    fseeko64(clusterfile_,rowsize,SEEK_CUR);
    if (i%1000==0) cout << "skipping " << nrows << " lines: " << Utils::ConvertFixed(i*100/(double)(nrows),1) << "%\r";
  }
  cout << endl;
}

double Map::BlockAverage(double **rows,int width,int startx,int xblock,int yblock) const
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

void Map::RowReduce(ofstream &stream,double **rows,int width,int xblock,int yblock)
{
  int nblocks=ceil(width/(double)xblock);
  for (int block=0;block<nblocks;block++)
  {
    double avg=BlockAverage(rows,width,block*xblock,xblock,yblock);
    stream << std::to_string(avg);
    if (block<nblocks-1) stream<<",";
  }
  stream << std::endl;
}

// transforms a value from the range 0..100 to nclasses in 0..255
int Map::TransformVal(double val,int nClasses)
{
  const double classWidth=256/(double)nClasses;
  const double outVal=(val*255.)/100.;
  int iVal=int(int(outVal/classWidth)*classWidth+classWidth/2);
  return iVal;
}

void Map::RowReduce(IMGPGM &PGMFile,double **rows,int width,int xblock,int yblock)
{
  int nblocks=ceil(width/(double)xblock);
  for (int block=0;block<nblocks;block++)
  {
    double avg=BlockAverage(rows,width,block*xblock,xblock,yblock);
    int ival=TransformVal(avg*100.,mapscale_);
    if (ival<0 || ival>255) cout << "map: warning: ival is out of range\n";

    PGMFile.rowbuffer[block]=uint8_t(ival);
  }
  PGMFile.WriteRow();
}

void Map::ProcessRows(std::string stroutfile,int ptop,int pbottom,int pleft,int pright)
{
  bool write_pgm=false;

  // check file extension for pgm
  std::size_t found=stroutfile.find_last_of('.');
  if (found!=std::string::npos) write_pgm=StringUtils::toupper(stroutfile.substr(found+1))=="PGM";

  int pwidth=pright-pleft;
  int pheight=pbottom-ptop;
  int rwidth=ceil(pwidth/(double)reduction_factor_);
  int rheight=ceil(pheight/(double)reduction_factor_);

  IMGPGM myPGM;
  ofstream rfile;
  if (write_pgm) {
    cout << "PGM: " << rwidth << "x" << rheight << endl;
    myPGM.SetWidth(rwidth);
    myPGM.SetHeight(rheight);
    myPGM.Create(stroutfile);
    myPGM.StartReader(); //my api sucks
  } else {
    cout << "DAT: " << rwidth << "x" << rheight << endl;
    rfile.open(stroutfile);
  }

  uint8_t tbuf[4];
  int64_t total_cells=0;
  int rowcnt=0;
  int outputrows=0;

  for (int row=ptop;row<pbottom;row++) {
    size_t tread=fread(tbuf,1,4,clusterfile_);
    uint32_t rowsize=Utils::Get32LH(tbuf);
    if (rowdata_.size()<rowsize) rowdata_.resize(rowsize);
    tread=fread(&rowdata_[0],1,rowsize,clusterfile_);
    if (tread!=rowsize) cout << "warning: could not read\n";
    else {
      RLEPack::UnpackRow(rowdata_,width_,labelrow_);
      double *drow=datarows_[rowcnt];
      int j=0;
      for (int i=pleft;i<pright;i++) {
        drow[j]=0.;
        if (labelrow_[i]) {
          total_cells++;

          if (labelrow_[i]>maxlabel_) cout << "warning: label outside range\n";
          else drow[j]=reflabels_[labelrow_[i]];
        };
        j++;
      }
    }
    rowcnt++;
    if (rowcnt>=reduction_factor_) {
      cout << (row+1) << "/" << height_ <<"\r";
      if (write_pgm) RowReduce(myPGM,datarows_,pwidth,reduction_factor_,rowcnt);
      else RowReduce(rfile,datarows_,pwidth,reduction_factor_,rowcnt);
      rowcnt=0;outputrows++;
    }
  }
  if (rowcnt) {
    if (write_pgm) RowReduce(myPGM,datarows_,pwidth,reduction_factor_,rowcnt);
    else RowReduce(rfile,datarows_,pwidth,reduction_factor_,rowcnt);
    outputrows++;
  };
  cout << "total cells:    " << total_cells << endl;
  cout << "rows processed: " << outputrows << endl;

  if (write_pgm) {
    myPGM.StopReader();
    myPGM.Close();
  } else {
    rfile.close();
  }
}

void Map::CalculateMap(const std::string &str_ifile,const std::string &str_ofile,const geoExtend &myExtend)
{
  std::string str_lfile,str_pfile,str_outfile;
  Utils::ReplaceExt(str_ifile,str_lfile,".lab",true);
  Utils::ReplaceExt(str_ifile,str_pfile,".txt",true);

  if (!str_ofile.length()) { // if no outputfile given, assume .dat
    Utils::ReplaceExt(str_ifile,str_outfile,".dat",true);
  } else str_outfile=str_ofile;

  cout << "reading clusters from '" << str_ifile << "': ";
  clusterfile_=fopen(str_ifile.c_str(),"rb");
  if (clusterfile_!=NULL) {
    uint8_t tbuf[8];
    fread(tbuf,1,8,clusterfile_);

    width_=Utils::Get32LH(tbuf);
    height_=Utils::Get32LH(tbuf+4);
    cout << width_ << "x" << height_ << endl;

    Projection Proj(width_,height_);
    Proj.ReadCoordinateFile(str_pfile);
    Proj.CalculateCellSize();

    int pleft,pright,ptop,pbottom;
    Frame::SetExtend(Proj.getLeft(),Proj.getTop(),Proj.getCellsize(),myExtend,width_,height_,pleft,ptop,pright,pbottom);

    cout << "map: writing to '" << str_outfile << "'\n";
    Timer myTimer;
    myTimer.Start();
    AllocMem();
    cout << "map: reduction factor: 1:" << reduction_factor_ << ", array: " << (((int64_t)reduction_factor_*(int64_t)width_*sizeof(double))>>20) << " mb\n";
    cout << "map: processing window (" << ptop << "," << pleft << ")x(" << pbottom << "," << pright << ")\n";
    cout << "map: scale " << mapscale_ << ", type ";
    switch (maptype_) {
      case 0: cout << "closs";break;
      case 1: cout << "core/area";break;
      default:cout << "unknown";break;
    }
    if (maptype_==1) cout << ", edge effect dept " << edge_effect_dept_ << "m";
    cout << endl;

    ReadLabels(str_lfile);
    TransferLabels();

    if (ptop) SkipRows(ptop-1);

    ProcessRows(str_outfile,ptop,pbottom,pleft,pright);

    fclose(clusterfile_);
    FreeMem();
    myTimer.Stop();
    cout << "time elapsed:   " << myTimer.ElapsedS() << " seconds" << endl;
  } else cout << "error: could not open: '" << str_ifile << "'\n";
}

