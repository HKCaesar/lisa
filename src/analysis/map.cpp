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
      firstToken(str,delim);
    }
    // requires the iternal string buffer to end on '\0'
    void firstToken(std::string &str,const char *delim)
    {
      m_delim=delim;
      m_token=strtok(&str[0],delim);
    }
    void nextToken()
    {
      m_token=strtok(NULL,m_delim);
    }
    bool valid(){return m_token!=NULL;};
    double getDouble(){return atof(m_token);};
    long long int getLongLong(){return strtoll(m_token,NULL,10);};
  private:
    const char *m_delim;
    char *m_token;
};

int Map::readLabels(const string &strfile)
{
  m_maxLabel=0;
  m_labels.clear();

  ifstream ifile(strfile);
  if (ifile.is_open()) {
    int64_t fileLength=Utils::GetStreamSize(ifile);
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
        cout << "reading labels: '" << strfile << "': " << Utils::ConvertFixed(ifile.tellg()*100.0/(double)fileLength,1) << "%\r";
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
        labelentry.label=tok.getLongLong();
        tok.nextToken(); // skip area in pixels
        tok.nextToken();
        labelentry.area=tok.getDouble();
        tok.nextToken();
        labelentry.edgeLen=tok.getDouble();
        tok.nextToken();
        labelentry.closs=tok.getDouble();
        if (!tok.valid()) cout << "map: unexpected number of tokens in file: '" << strfile << "'\n";
      #endif

      total_area+=labelentry.area;
      total_edgelen+=labelentry.edgeLen;
      if (labelentry.label>m_maxLabel) m_maxLabel=labelentry.label;
      m_labels.push_back(labelentry);
    }
    ifile.clear();
    cout << "reading labels: '" << strfile << "': " << Utils::ConvertFixed(ifile.tellg()*100.0/(double)fileLength,1) << "%\n";
    ifile.close();
    myTimer.Stop();
    cout << "number of labels: " << m_labels.size() << "\n";
    cout << "total area:       " << Utils::ConvertFixed(Utils::SqMetre_To_MillHa(total_area),2) << " 10^6 ha" << endl;
    cout << "edge len:         " << Utils::ConvertFixed(Utils::Metre_To_MillKm(total_edgelen),2) << " 10^6 km" << endl;
    cout << "time elapsed:     " << myTimer.ElapsedS() << " seconds" << endl;
    return 0;
  } else return 1;
}

void Map::transferLabels()
{
  m_refLabels.resize(m_maxLabel+1);
  for (size_t i=0,ilen=m_labels.size();i<ilen;i++) {
    double val=0.;
    if (m_mapType==0) val=m_labels[i].closs;
    else if (m_mapType==1) {
      double edge_area=Cluster::CalculateEdgeAreaDE(m_labels[i].area,m_labels[i].edgeLen,m_edgeEffectDept);
      double core_area=m_labels[i].area-edge_area;
      val=core_area/m_labels[i].area;
    }
    m_refLabels[m_labels[i].label]=val;
  }
  m_labels.clear();
}

void Map::allocMem()
{
  m_rowData=new uint8_t[m_maxRowDataSize];
  m_labelRow=new int64_t[m_width];
  m_dataRows=new double*[m_reductionFactor];
  for (int i=0;i<m_reductionFactor;i++) m_dataRows[i]=new double[m_width];
}

void Map::freeMem()
{
  for (int i=0;i<m_reductionFactor;i++) delete []m_dataRows[i];
  delete []m_dataRows;
  delete []m_rowData;
  delete []m_labelRow;
}

void Map::skipRows(int nrows)
{
  uint8_t tbuf[4];
  for (int i=0;i<nrows;i++) {
    size_t tread=fread(tbuf,1,4,m_clusterFile);
    uint32_t rowsize=Utils::Get32LH(tbuf);
    fseek(m_clusterFile,rowsize,SEEK_CUR);
    if (i%1000==0) cout << "skipping " << nrows << " lines: " << Utils::ConvertFixed(i*100/(double)(nrows),1) << "%\r";
  }
  cout << endl;
}

double Map::blockAverage(double **rows,int width,int startx,int xblock,int yblock) const
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

void Map::rowReduce(ofstream &stream,double **rows,int width,int xblock,int yblock)
{
  int nblocks=ceil(width/(double)xblock);
  for (int block=0;block<nblocks;block++)
  {
    double avg=blockAverage(rows,width,block*xblock,xblock,yblock);
    stream << std::to_string(avg);
    if (block<nblocks-1) stream<<",";
  }
  stream << std::endl;
}

// transforms a value from the range 0..100 to nclasses in 0..255
int Map::transformVal(double val,int nClasses)
{
  const double classWidth=256/(double)nClasses;
  const double outVal=(val*255.)/100.;
  int iVal=int(int(outVal/classWidth)*classWidth+classWidth/2);
  return iVal;
}

void Map::rowReduce(IMGPGM &PGMFile,double **rows,int width,int xblock,int yblock)
{
  int nblocks=ceil(width/(double)xblock);
  for (int block=0;block<nblocks;block++)
  {
    double avg=blockAverage(rows,width,block*xblock,xblock,yblock);
    int ival=transformVal(avg*100.,m_mapScale);
    if (ival<0 || ival>255) cout << "map: warning: ival is out of range\n";

    PGMFile.rowbuffer[block]=uint8_t(ival);
  }
  PGMFile.WriteRow();
}

void Map::processRows(std::string stroutfile,int ptop,int pbottom,int pleft,int pright)
{
  bool write_pgm=false;

  // check file extension for pgm
  std::size_t found=stroutfile.find_last_of('.');
  if (found!=std::string::npos) write_pgm=StringUtils::toupper(stroutfile.substr(found+1))=="PGM";

  int pwidth=pright-pleft;
  int pheight=pbottom-ptop;
  int rwidth=ceil(pwidth/(double)m_reductionFactor);
  int rheight=ceil(pheight/(double)m_reductionFactor);

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
    size_t tread=fread(tbuf,1,4,m_clusterFile);
    uint32_t rowsize=Utils::Get32LH(tbuf);
    if (rowsize>=m_maxRowDataSize) {
      cout << "map: rowsize >= maxrowdatasize, skipping line\n";
      fseek(m_clusterFile,rowsize,SEEK_CUR);
    } else {
      tread=fread(m_rowData,1,rowsize,m_clusterFile);
      if (tread!=rowsize) cout << "warning: could not read\n";
      else {
        RLEPack::UnpackRow(m_rowData,m_width,m_labelRow);
        double *drow=m_dataRows[rowcnt];
        int j=0;
        for (int i=pleft;i<pright;i++) {
          drow[j]=0.;
          if (m_labelRow[i]) {
            total_cells++;

            if (m_labelRow[i]>m_maxLabel) cout << "warning: label outside range\n";
            else drow[j]=m_refLabels[m_labelRow[i]];
          };
          j++;
        }
      }
      rowcnt++;
      if (rowcnt>=m_reductionFactor) {
        cout << (row+1) << "/" << m_height <<"\r";
        if (write_pgm) rowReduce(myPGM,m_dataRows,pwidth,m_reductionFactor,rowcnt);
        else rowReduce(rfile,m_dataRows,pwidth,m_reductionFactor,rowcnt);
        rowcnt=0;outputrows++;
      }
    }
  }
  if (rowcnt) {
    if (write_pgm) rowReduce(myPGM,m_dataRows,pwidth,m_reductionFactor,rowcnt);
    else rowReduce(rfile,m_dataRows,pwidth,m_reductionFactor,rowcnt);
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

void Map::calculateMap(const std::string &str_ifile,const std::string &str_ofile,const geoExtend &myExtend)
{
  std::string str_lfile,str_pfile,str_outfile;
  Utils::ReplaceExt(str_ifile,str_lfile,".lab",true);
  Utils::ReplaceExt(str_ifile,str_pfile,".txt",true);

  if (!str_ofile.length()) { // if no outputfile given, assume .dat
    Utils::ReplaceExt(str_ifile,str_outfile,".dat",true);
  } else str_outfile=str_ofile;

  cout << "reading clusters from '" << str_ifile << "': ";
  m_clusterFile=fopen(str_ifile.c_str(),"rb");
  if (m_clusterFile!=NULL) {
    uint8_t tbuf[8];
    fread(tbuf,1,8,m_clusterFile);

    m_width=Utils::Get32LH(tbuf);
    m_height=Utils::Get32LH(tbuf+4);
    m_maxRowDataSize=m_width*10;
    cout << m_width << "x" << m_height << endl;

    Projection Proj(m_width,m_height);
    Proj.ReadCoordinateFile(str_pfile);
    Proj.CalculateCellSize();

    int pleft,pright,ptop,pbottom;
    Frame::SetExtend(Proj.getLeft(),Proj.getTop(),Proj.getCellsize(),myExtend,m_width,m_height,pleft,ptop,pright,pbottom);

    cout << "map: writing to '" << str_outfile << "'\n";
    Timer myTimer;
    myTimer.Start();
    allocMem();
    cout << "map: reduction factor: 1:" << m_reductionFactor << ", array: " << (((int64_t)m_reductionFactor*(int64_t)m_width*sizeof(double))>>20) << " mb\n";
    cout << "map: processing window (" << ptop << "," << pleft << ")x(" << pbottom << "," << pright << ")\n";
    cout << "map: scale " << m_mapScale << ", type ";
    switch (m_mapType) {
      case 0: cout << "closs";break;
      case 1: cout << "core/area";break;
      default:cout << "unknown";break;
    }
    if (m_mapType==1) cout << ", edge effect dept " << m_edgeEffectDept << "m";
    cout << endl;

    readLabels(str_lfile);
    transferLabels();

    if (ptop) skipRows(ptop-1);

    processRows(str_outfile,ptop,pbottom,pleft,pright);

    fclose(m_clusterFile);
    freeMem();
    myTimer.Stop();
    cout << "time elapsed:   " << myTimer.ElapsedS() << " seconds" << endl;
  } else cout << "error: could not open: '" << str_ifile << "'\n";
}

