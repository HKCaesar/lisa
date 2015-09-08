#include "pgm.h"

void PGM::ReadPGMLine(std::string &line)
{
  char c;
  line.clear();
  do {
    c=getc(file);
  } while (Utils::iswhite(c));
  do
  {
    line+=c;
    c=getc(file);
    //cout << (int)c << " ";
  } while (!Utils::iswhite(c));

}

void PGM::SeekPGMLine(off64_t y)
{
  off64_t pos=hdrpos+((off64_t)width*y);
  fseeko64(file,pos,SEEK_SET);
}


int PGM::ReadHeader()
{
  std::string line;
  fsize=Utils::GetFileSize(file);
  ReadPGMLine(line);
  if (line.compare("P5")==0) { //PGM file found
    ReadPGMLine(line);
    width=std::stoi(line);
    ReadPGMLine(line);
    height=std::stoi(line);
    ReadPGMLine(line);
    max_gray=std::stoi(line);
    hdrpos=ftello64(file);
       // last line must end at filesize
    SeekPGMLine(height);
    valid=(ftello64(file)==fsize);
    fseeko64(file,hdrpos,SEEK_SET);
    return 0;
  } else {
    fseeko64(file,0,SEEK_SET);
    return 1;
  }
}

int PGM::Create(std::string &str_ofile)
{
  file=fopen(str_ofile.c_str(),"wb");
  if (file)
  {
    std::string s;
    fputc('P',file);
    fputc('5',file);
    fputc(10,file);
    s=to_string(width);
    for (size_t i=0;i<s.length();i++) fputc(s[i],file);
    fputc(32,file);
    s=to_string(height);
    for (size_t i=0;i<s.length();i++) fputc(s[i],file);
    fputc(10,file);
    s=string("255");
    for (size_t i=0;i<s.length();i++) fputc(s[i],file);
    fputc(10,file);
    hdrpos=ftello64(file);
    return 0;
  } else return 1;
}

void PGM::Close()
{
  fclose(file);
}

void PGM::PrintInfo()
{
  cout << endl << "PGM";
  if (valid) cout << "(valid)";
  else cout << "(invalid)";
  cout << ": " << width << "x" << height << ": " << max_gray << endl;
  cout << "filesize: " << (fsize>>20) << " MB, data start: " << hdrpos << " bytes" << endl;
}

void PGM::StartReading()
{
  linebuf=new uint8_t[width];
  if (ftello64(file)!=hdrpos) fseeko64(file,hdrpos,SEEK_SET);
}

void PGM::StopReading()
{
  if (linebuf) delete []linebuf,linebuf=0;
}

bool PGM::ReadRow()
{
  size_t nread=fread(linebuf,1,width,file);
  if (nread!=(size_t)width) {cout << "error: could not read (" << nread << " from " << width << ")" << endl;return false;};
  return true;
}

bool PGM::WriteRow()
{
  size_t nwrite=fwrite(linebuf,1,width,file);
  if (nwrite!=(size_t)width) {cout << "error: could not write (" << nwrite << " from " << width << ")" << endl;return false;};
  return true;
}
