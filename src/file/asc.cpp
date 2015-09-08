#include "asc.h"

bool ASC::getLine(std::string &line)
{
  int ch;
  line.clear();
  while (1)
  {
    ch=getc(file);
    switch (ch) {
      case 13:break;
      case 10:return true;
      case EOF:return false;
      default: line+=ch;break;
    }
  };
  return true;
}

void ASC::PrintInfo()
{
  std::cout << "ASC: " << width << "x" << height << endl;
  GeoUtils::PrintInfo(top,left,right,bottom,cellsize);
}

int ASC::ReadHeader()
{
  std::string ukey,key,val;
  width=height=pleft=pright=ptop=pbottom=0;

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("NCOLS")==0) width=std::stoi(val);
  else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("NROWS")==0) height=std::stoi(val);
  else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("XLLCORNER")==0) left=std::stod(val);
  else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("YLLCORNER")==0) bottom=std::stod(val);
  else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("CELLSIZE")==0) cellsize=std::stod(val);
  else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  // should be fixed, because NODATA_VALUE is optional
  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("NODATA_VALUE")==0) nodataval=std::stoi(val);
  else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  const double delta_lat=cellsize*(double)height;
  const double delta_long=cellsize*(double)width;
  right=left+delta_long;
  top=bottom+delta_lat;

  pleft=0;pright=width;
  ptop=0;pbottom=height;

  tleft=left;tright=right;ttop=top;tbottom=bottom;

  return 0;
}

int ASC::getLongPos(double ref_pos,double left)
{
  return round((ref_pos-left)/cellsize);
}

int ASC::getLatPos(double ref_pos,double top)
{
  return round((-ref_pos+top)/cellsize);
}

int ASC::WriteExtend(const std::string &fname,int prec)
{
  ofstream ofile(fname);
  if (ofile.is_open()) {
    ofile << "ref=wgs84\n";
    ofile << "top=" << Utils::ConvertFixed(ttop,prec)<<"\n";
    ofile << "left=" << Utils::ConvertFixed(tleft,prec)<<"\n";
    ofile << "right=" << Utils::ConvertFixed(tright,prec)<<"\n";
    ofile << "bottom=" << Utils::ConvertFixed(tbottom,prec)<<"\n";
    return 0;
  } else return 1;
}

// reframe the extend
void ASC::SetExtend(const geoExtend &myExtend)
{
if ( (myExtend.top-myExtend.bottom>0.) && (myExtend.right-myExtend.left>0.))
  {
    std::cout << "\nreframing extend..." << endl;
    pleft=getLongPos(myExtend.left,left);
    pright=getLongPos(myExtend.right,left);
    ptop=getLatPos(myExtend.top,top);
    pbottom=getLatPos(myExtend.bottom,top);
    cout << "x: [" << pleft << "," << pright << "], y: [" << ptop << "," << pbottom << "]\n";
    double tleft=(left+pleft*cellsize);
    double tright=(left+pright*cellsize);
    double ttop=(top-ptop*cellsize);
    double tbottom=(top-pbottom*cellsize);
    cout << std::fixed << std::setprecision(8) << "Top: " << ttop << ", Left: " << tleft << ", Right: " << tright << ", Bottom: " << tbottom << endl;
    if (pleft<0 || pleft>width || pright<0 || pright>width || ptop<0 || ptop>height || pbottom<0 || pbottom>height) {
      cout << "error: cannot match extends\n";
      pleft=0;pright=width;
      ptop=0;pbottom=height;
    }
  }
}

void ASC::StopReading()
{
  if (linebuf) delete []linebuf,linebuf=0;
}

void ASC::StartReading()
{
  hist[0]=hist[1]=0;
  linebuf=new uint8_t[getExtendWidth()];
  vecline.resize(width);
  linenum=0;

  for (int i=0;i<ptop-1;i++) {
    getLine(line);
    if (linenum%10==0) cout << "skipping " << (ptop-1) << " lines: " << Utils::ConvertFixed(i*100/(double)(ptop-1),1) << "%\r";
    linenum++;
  }
}

int ASC::ReadRow()
{
  getLine(line);
  if (!Utils::SplitTokenInt(line,vecline)) cout << "error reading file at line: " << linenum << endl;
  int j=0;
  for (int i=pleft;i<pright;i++)
  {
      linebuf[j]=0;
      const int val=vecline[i];
      // binary read
      if (val==0 || val==1) {linebuf[j]=val;hist[val]++;}
      else if (val!=nodataval) cout << "unknown pixel value: " << vecline[i] << " at line " << linenum << ", column " << i << endl;

      #if 0
        if (val==99) linebuf[j]=0;
        else if (val==11) linebuf[j]=1;
        else if (val!=nodataval) cout << "unknown pixel value: " << vecline[i] << "at line " << linenum << ", column " << i << endl;
      #endif // 0const int val=vecline[i];
      j++;
  }
  linenum++;
  return 0;
}
