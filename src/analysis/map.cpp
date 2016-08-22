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
    cout << "time elapsed:     " << Utils::Sec2Time(round(myTimer.ElapsedS())) << endl;
    return 0;
  } else return 1;
}

void Map::TransferLabels()
{
  reflabels_.resize(maxlabel_+1);
  for (auto &label : labels_) {
    double val=0.;
    if (maptype_==0) val=label.closs;
    else if (maptype_==1) {
      double edge_area=Cluster::CalculateEdgeAreaDE(label.area,label.edgelen,edge_effect_dept_);
      double core_area=label.area-edge_area;
      val=core_area/label.area;
    }
    reflabels_[label.label]=val;
  }
  labels_.clear();
}

void Map::AllocMem()
{
  labelrow_=new int64_t[width_];
  datarows_.resize(reduction_factor_,vector<RowLabel>(width_));
  /*datarows_=new double*[reduction_factor_];
  for (int i=0;i<reduction_factor_;i++) datarows_[i]=new double[width_];*/
}

void Map::FreeMem()
{
  //for (int i=0;i<reduction_factor_;i++) delete []datarows_[i];
  //delete []datarows_;
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

void Map::RowReduce(ofstream &stream,int width,int xblock,int yblock)
{
  int nblocks=ceil(width/(double)xblock);
  for (int block=0;block<nblocks;block++)
  {
    double avg=BlockAverageClassified(width,block*xblock,xblock,yblock);
    stream << std::to_string(avg);
    if (block<nblocks-1) stream<<",";
  }
  stream << std::endl;
}

// transforms a value from the range 0..100 to nclasses in offset..255
int Map::TransformVal(double val,int offset,int nClasses)
{
  if (val>0.) {
    const double classWidth=(256-offset)/(double)nClasses;
    const double outVal=(val*(255-offset))/100.;
    int iVal=offset+int(int(outVal/classWidth)*classWidth+classWidth/2);
    return iVal;
  } else return 0;
}

double Map::BlockAverageClassified(int width,int startx,int xblock,int yblock) const
{
  double sum=0.0;
  int num_below=0;
  int num_above=0;
  for (int j=0;j<yblock;j++) {
    const vector <RowLabel> &drow=datarows_[j];
    for (int i=startx;i<std::min(width,startx+xblock);i++) {
      if (drow[i].labeled) { // is it classified as forest?
        sum+=drow[i].data;
        num_above++;
      } else num_below++;
    }
  }
  if (num_above>=num_below) {
    return sum/(double)num_above;
  } else return 0.;
}

void Map::RowReduce(IMGPGM &PGMFile,int width,int xblock,int yblock)
{
  int nblocks=ceil(width/(double)xblock);
  for (int block=0;block<nblocks;block++)
  {
    double avg=BlockAverageClassified(width,block*xblock,xblock,yblock);
    int ival=round(avg*100.);
    if (ival<0 || ival>100) cout << "map: warning: ival is out of range\n";

    PGMFile.rowbuffer[block]=uint8_t(ival);
  }
  PGMFile.WriteRow();
}

int Map::GetMapWidth(int pwidth)
{
  return ceil(pwidth/(double)reduction_factor_);
}

int Map::GetMapHeight(int pheight)
{
  return ceil(pheight/(double)reduction_factor_);
}

void Map::ProcessRows(const std::string &stroutfile,int ptop,int pbottom,int pleft,int pright)
{
  bool write_pgm=false;

  // check file extension for pgm
  std::size_t found=stroutfile.find_last_of('.');
  if (found!=std::string::npos) write_pgm=StringUtils::toupper(stroutfile.substr(found+1))=="PGM";

  int pwidth=pright-pleft;
  int pheight=pbottom-ptop;

  int rwidth=GetMapWidth(pwidth);
  int rheight=GetMapHeight(pheight);

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

  int64_t zero_cells=0;

  for (int row=ptop;row<pbottom;row++) {
    size_t tread=fread(tbuf,1,4,clusterfile_);
    uint32_t rowsize=Utils::Get32LH(tbuf);
    if (rowdata_.size()<rowsize) rowdata_.resize(rowsize);
    tread=fread(&rowdata_[0],1,rowsize,clusterfile_);
    if (tread!=rowsize) cout << "warning: could not read\n";
    else {
      RLEPack::UnpackRow(rowdata_,width_,labelrow_);
      vector <RowLabel> &drow=datarows_[rowcnt];
      int j=0;
      for (int i=pleft;i<pright;i++) {
        drow[j].labeled=false;
        drow[j].data=0.;
        if (labelrow_[i]) {
          drow[j].labeled=true;
          total_cells++;

          if (labelrow_[i]>maxlabel_) cout << "warning: label outside range\n";
          else drow[j].data=reflabels_[labelrow_[i]];
          if (drow[j].data==0.) zero_cells++;
        };
        j++;
      }
    }
    rowcnt++;
    if (rowcnt>=reduction_factor_) {
      cout << (row+1) << "/" << height_ <<"\r";
      if (write_pgm) RowReduce(myPGM,pwidth,reduction_factor_,rowcnt);
      else RowReduce(rfile,pwidth,reduction_factor_,rowcnt);
      rowcnt=0;outputrows++;
    }
  }
  if (rowcnt) {
    if (write_pgm) RowReduce(myPGM,pwidth,reduction_factor_,rowcnt);
    else RowReduce(rfile,pwidth,reduction_factor_,rowcnt);
    outputrows++;
  };
  cout << "total cells:    " << total_cells << " (null cells: " << zero_cells << ")" << endl;
  cout << "rows processed: " << outputrows << endl;

  if (write_pgm) {
    myPGM.StopReader();
  } else {
    rfile.close();
  }
}

void Map::PrintColorTable(int offset)
{
    cout << "colortable" << endl;
    int val=0;
    int color_code=TransformVal(val,offset,mapscale_);
    cout << "[" << val;
    for (val=1;val<=100;val++) {
       int color_code_next=TransformVal(val,offset,mapscale_);
       if (color_code_next!=color_code) {
          cout << "," << (val-1) << "]:" << color_code << endl;
          color_code=color_code_next;
          cout << "[" << val;
       }
    }
    cout << "," << (val-1) << "]:" << color_code << endl;
}

void Map::Create(const std::string &str_ifile,const std::string &str_ofile,const geoExtend &myExtend)
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
    cout << "map: reduction factor: 1:" << reduction_factor_ << ", array: " << (((int64_t)reduction_factor_*(int64_t)width_*sizeof(RowLabel))>>20) << " mb"<<endl;
    cout << "map: processing window (" << ptop << "," << pleft << ")x(" << pbottom << "," << pright << ")" << endl;
    cout << "map: type ";
    switch (maptype_) {
      case 0: cout << "closs";break;
      case 1: cout << "core/area";break;
      default:cout << "unknown";break;
    }
    if (maptype_==1) cout << " (edge effect dept " << edge_effect_dept_ << "m)";
    cout << ",size " << GetMapWidth(pright-pleft) << "x" << GetMapHeight(pbottom-ptop)<<endl;
    cout << endl;

    Timer myTimer;
    myTimer.Start();
    AllocMem();

    if (ptop) SkipRows(ptop-1);

    ReadLabels(str_lfile);
    TransferLabels();

    ProcessRows(str_outfile,ptop,pbottom,pleft,pright);

    FreeMem();
    myTimer.Stop();
    fclose(clusterfile_);
    cout << "time elapsed:   " << Utils::Sec2Time(round(myTimer.ElapsedS())) << endl;
  } else cout << "error: could not open: '" << str_ifile << "'\n";
}


#include "../file/tiff.h"

int Map::BlockAverageClass(const vector <vector<uint8_t>> &rows,int width,int startx,int xblock,int yblock) const
{
  int nclasses=3;
  vector <int>nclass(nclasses);
  for (int j=0;j<yblock;j++) {
    for (int i=startx;i<std::min(width,startx+xblock);i++) {
      int c=rows[j][i];
      if (c>=0 && c<nclasses) nclass[c]++;
      else cerr << "map: warning class index out of bounds.\n";
    }
  }
  int maxclass=0;
  int maxval=nclass[0];
  for (int i=1;i<nclasses;i++) {
    if (nclass[i]>=maxval) {
        maxval=nclass[i];
        maxclass=i;
    }
  }
  return maxclass;
}

void Map::RowReduceClass(IMGPGM &PGMFile,const vector <vector<uint8_t>> &rows,int width,int xblock,int yblock)
{
  int nblocks=ceil(width/(double)xblock);
  for (int block=0;block<nblocks;block++)
    PGMFile.rowbuffer[block]=BlockAverageClass(rows,width,block*xblock,xblock,yblock);

    /*if (avg_class==1) PGMFile.rowbuffer[block]=0;
    else PGMFile.rowbuffer[block]=255;*/
  PGMFile.WriteRow();
}

void Map::Reduce(const std::string &str_ifile,const std::string &str_ofile)
{
  IMGTIFF myTIFF;
  IMG *myIMG=nullptr;
  cout << "reading '" << str_ifile << "'" << endl;
  if (myTIFF.ReadHeader(str_ifile.c_str())==0)
  {
    myIMG=&myTIFF;
    myIMG->StartReader();
    int width=myIMG->GetWidth();
    int height=myIMG->GetHeight();
    myIMG->PrintInfo();

    int rwidth=GetMapWidth(width);
    int rheight=GetMapHeight(height);

    IMGPGM myPGM;
    cout << "PGM: " << rwidth << "x" << rheight << endl;
    myPGM.SetWidth(rwidth);
    myPGM.SetHeight(rheight);
    myPGM.Create(str_ofile);
    myPGM.StartReader(); //my api sucks

    std::vector <std::vector<uint8_t>> datarows(reduction_factor_,vector <uint8_t>(width));
    int rowcnt=0;
    for (int row=0;row<height;row++) {
        if (row%1000==0) cout << row << '\r';
        uint8_t *rowbuf=&datarows[rowcnt][0];
        myIMG->ReadRow(rowbuf);

        rowcnt++;
        if (rowcnt>=reduction_factor_) {
          RowReduceClass(myPGM,datarows,width,reduction_factor_,rowcnt);
          rowcnt=0;
        }
        /*for (int i=0;i<width;i++) {
            cout << (int)rowbuf[i];
        }*/
        //cout << endl;
    }
    if (rowcnt) RowReduceClass(myPGM,datarows,width,reduction_factor_,rowcnt);
    myPGM.StopReader();

    myIMG->StopReader();
  }
}

void Map::Classify(const std::string &str_ifile,const std::string &str_maskfile,const std::string &str_ofile)
{
  IMGPGM dataPGM,maskPGM,outPGM;
  if (!dataPGM.OpenRead(str_ifile)) {cerr << "could not open: '" << str_ifile << "'" << endl;return;};
  if (!maskPGM.OpenRead(str_maskfile)) {cerr << "could not open: '" << str_maskfile << "'" << endl;return;};
  if (dataPGM.ReadHeader()!=0) {cerr << "no valid image file: '" << str_ifile << "'" << endl;return;};
  if (maskPGM.ReadHeader()!=0) {cerr << "no valid image file: '" << str_maskfile << "'" << endl;return;};
  if ( (dataPGM.GetWidth()!=maskPGM.GetWidth()) || (dataPGM.GetHeight()!=maskPGM.GetHeight())) {cerr << "image dimensions do not match" << endl;return;};

  cout << "reading datafile: '" << str_ifile << "'";
  dataPGM.PrintInfo();
  cout << "reading maskfile: '" << str_maskfile << "'";
  maskPGM.PrintInfo();
  cout << "writing:          '" << str_ofile << "'" << endl;

  int width=dataPGM.GetWidth();
  int height=dataPGM.GetHeight();

  outPGM.SetWidth(width);
  outPGM.SetHeight(height);
  outPGM.SetType(IMGPGM::PGMRGB);
  outPGM.Create(str_ofile);

  dataPGM.StartReader();
  maskPGM.StartReader();
  outPGM.StartReader();

  int map_offset=64;
  if (mapclass_==0) cout << "map: 4 classes" << endl;
  else if (mapclass_==1) {
    cout << "map: scale " << mapscale_ << ", offset " << map_offset << endl;
    PrintColorTable(map_offset);
  }

  for (int row=0;row<height;row++) {
    if (row%1000==0) cout << row << '\r';
    dataPGM.ReadRow();
    maskPGM.ReadRow();
    for (int i=0;i<width;i++) {
        int dataval=dataPGM.rowbuffer[i];
        int maskval=maskPGM.rowbuffer[i];

        int rval=0;
        int gval=0;
        int bval=0;

        if (dataval>0) {
          if (dataval>100) cerr << "map: input val out of range: " << dataval << endl;
          if (mapclass_==0) {
            if (dataval<=25) {rval=255;gval=0;bval=0;}
            else if (dataval<=50) {rval=255;gval=128;bval=0;}
            else if (dataval<=75) {rval=255;gval=255;bval=0;}
            else {rval=0;gval=255;bval=0;}
          } else if (mapclass_==1) {
            gval=TransformVal(dataval,map_offset,mapscale_);
            if (gval<0||gval>255) cerr << "map: output val is out of range\n";
          }
          /*if (dataval==88) {rval=255;gval=0;bval=0;}
          else if (dataval==136) {rval=255;gval=128;bval=0;}
          else if (dataval==184) {rval=255;gval=255;bval=0;}
          else if (dataval==232) {rval=0;gval=255;bval=0;};*/
          //gval=dataval;
        } else {
          if (maskval==0 || maskval==2) { // blue for 0-nodata or 2-permanent water bodies
            rval=50;gval=150;bval=255;
          }
        }

        outPGM.rowbuffer[3*i]=rval;
        outPGM.rowbuffer[3*i+1]=gval;
        outPGM.rowbuffer[3*i+2]=bval;
    }
    outPGM.WriteRow();
  }
  outPGM.StopReader();
  maskPGM.StopReader();
  dataPGM.StopReader();
  cout << height << '\n';
}

