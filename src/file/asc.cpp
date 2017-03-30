#include "asc.h"

bool IMGASC::getLine(std::string &line)
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

void IMGASC::PrintInfo()
{
  std::cout << "ASC: " << width << "x" << height << endl;
  //GeoUtils::PrintInfo(top,left,right,bottom,cellsize);
  if (nodata_avail) cout << "nodata val: " << nodata_value << endl;
  GeoUtils::PrintInfo(top,left,right,bottom,cellsize);
}

int IMGASC::ReadHeader()
{
  std::string ukey,key,val;
  width=height=pleft=pright=ptop=pbottom=0;

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("NCOLS")==0) width=std::stoi(val);
  else return 1;
  //else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("NROWS")==0) height=std::stoi(val);
  else return 1;
  //else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("XLLCORNER")==0) left=std::stod(val);
  else return 1;
  //else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("YLLCORNER")==0) bottom=std::stod(val);
  else return 1;
  //else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("CELLSIZE")==0) cellsize=std::stod(val);
  else return 1;
  //else {Utils::PrintWarning("unknown key: '" + key + "'");return 1;};

  // if we came this far, it's probably an .asc file
  off_t posnodata=ftello64(file);

  getLine (line);Utils::Split(line,key,val);ukey=StringUtils::toupper(key);
  if (ukey.compare("NODATA_VALUE")==0) {nodata_avail=true;nodata_value=std::stoi(val);}
  else fseeko64(file,posnodata,0); // no nodata value

  const double delta_lat=cellsize*(double)height;
  const double delta_long=cellsize*(double)width;
  right=left+delta_long;
  top=bottom+delta_lat;

  pleft=0;pright=width;
  ptop=0;pbottom=height;

  tleft=left;tright=right;ttop=top;tbottom=bottom;

  return 0;
}

int IMGASC::getLongPos(double ref_pos,double left)
{
  return round((ref_pos-left)/cellsize);
}

int IMGASC::getLatPos(double ref_pos,double top)
{
  return round((-ref_pos+top)/cellsize);
}

int IMGASC::WriteExtend(const std::string &fname,int prec)
{
  ofstream ofile(fname);
  if (ofile.is_open()) {
    ofile << "ref=wgs84\n";
    ofile << "top=" << Utils::ConvertFixed(ttop,prec)<<"\n";
    ofile << "left=" << Utils::ConvertFixed(tleft,prec)<<"\n";
    ofile << "right=" << Utils::ConvertFixed(tright,prec)<<"\n";
    ofile << "bottom=" << Utils::ConvertFixed(tbottom,prec)<<"\n";
    ofile.close();
    return 0;
  } else return 1;
}

// reframe the extend
void IMGASC::SetExtend(const geoExtend &myExtend)
{
  Frame::SetExtend(left,top,cellsize,myExtend,width,height,pleft,ptop,pright,pbottom);
}

void IMGASC::StopReader()
{
  if (rowbuffer) delete []rowbuffer,rowbuffer=0;
}

void IMGASC::StartReader()
{
  hist[0]=hist[1]=0;
  rowbuffer=new uint8_t[getExtendWidth()];
  vecline.resize(width);
  linenum=0;

  for (int i=0;i<ptop-1;i++) {
    getLine(line);
    if (linenum%10==0) cout << "skipping " << (ptop-1) << " lines: " << Utils::ConvertFixed(i*100/(double)(ptop-1),1) << "%\r";
    linenum++;
  }
}

int IMGASC::ReadRow(uint8_t *buf)
{
  getLine(line);
  if (!Utils::SplitTokenInt(line,vecline)) cout << "error reading file at line: " << linenum << endl;
  int j=0;
  for (int i=pleft;i<pright;i++)
  {
      buf[j]=0;
      const int val=vecline[i];
      // binary read
      if (val==0 || val==1) {buf[j]=val;hist[val]++;}
      else if (!nodata_avail || (nodata_avail && val!=nodata_value)) cout << "unknown pixel value: " << vecline[i] << " at line " << linenum << ", column " << i << endl;
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
