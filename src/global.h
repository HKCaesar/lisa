#ifndef GLOBAL_H
#define GLOBAL_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <list>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <sys/time.h>

using namespace std;

#define LISA_VERSION "Large Image Spatial Analysis v0.91 (c) 2014/2015 - Sebastian Lehmann"

#define BOOL(x) (!(!(x)))

#define BitSet(arg,posn) ((arg) | (1L << (posn)))
#define BitClear(arg,posn) ((arg) & ~(1L << (posn)))
#define BitTest(arg,posn) BOOL((arg) & (1L << (posn)))
#define BitFlip(arg,posn) ((arg) ^ (1L << (posn)))

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define EPS (1E-12)

template <class T> T clamp(T val,T min,T max) {return val<min?min:val>max?max:val;};
template <class T> T div_signed(T val,T s){return val<0?-(((-val)+(1<<(s-1)))>>s):(val+(1<<(s-1)))>>s;};

class geoExtend {
  public:
    geoExtend():top(0.),left(0.),right(0.),bottom(0.){};
    double top,left,right,bottom;
};

class cluster_stats {
  public:
      cluster_stats(){Reset();}
  void Reset(){
    cell_area=num_clusters=num_clusters10ha=num_clusters50ha=0;
    total_area=mean_area=total_border_len=total_edge_area_de=total_edge_area_circle=max_area=0.0;
    total_biomass=0.;
  };
  int64_t cell_area,num_clusters,num_clusters10ha,num_clusters50ha;
  double total_area,mean_area,total_border_len,total_edge_area_de,total_edge_area_circle,max_area;
  double total_biomass;
};

class Timer {
  public:
      void Start()
      {
        gettimeofday(&tim, NULL);
        t1=tim.tv_sec+(tim.tv_usec/1000000.0);
      }
      void Stop()
      {
        gettimeofday(&tim, NULL);
        t2=tim.tv_sec+(tim.tv_usec/1000000.0);
      }
      double SecElapsed(){return t2-t1;};
   protected:
      double t1,t2;
      struct timeval tim;
};

class GeoUtils {
  public:
      static void PrintInfo(double top,double left,double right,double bottom,double cellsize)
      {
        cout << "Top: " << top << ", Left: " << left << ", Right: " << right << ", Bottom: " << bottom << endl;
        cout << "dLat: " << deltalat(top,bottom) << " degree, dLong: " << deltalong(left,right) << " degree" << endl;
        cout << "Cellsize: " << cellsize << " degree = " << GeoUtils::cellsize2arcsec(cellsize) << " arcsec = " << GeoUtils::cellsize2m(cellsize) << " m" << endl;
      }
      static double deltalat(double top,double bottom)
      {
         return top-bottom;
      }
      static double deltalong(double left,double right)
      {
         return right-left;
      }
      static double cellsize2arcsec(double cellsize) {
        return cellsize*3600.;
      }
      // at the equator one arcsec approx. 1/60th of a nautical mile (or 101.27 feet or 30.87 meters)
      static double cellsize2m(double cellsize) {
        return cellsize*3600.*30.87;
      }
};

class Utils {
  public:
    static void PrintWarning(const std::string &s) {
      std::cerr << "warning: " << s << "\n";
    }
    static int PromptYesNo(const std::string &prompt)
    {
      char type;
      do {
        cout << prompt;
        cin >> type;
        type=std::tolower(type);
      } while (type!='y' && type!='n');
      if (type=='y') return 1;
      else return 0;
    }
    static FILE *OpenWriteFILE(const std::string &fname,bool force=false)
    {
      FILE *fin=fopen(fname.c_str(),"rb");
      if (fin!=NULL) {
        fclose(fin);
        if (force || PromptYesNo(std::string("warning: file '" + fname + "' already exists! Overwrite (y/n)?"))) fin=fopen(fname.c_str(),"wb");
        else fin=NULL;
      } else fin=fopen(fname.c_str(),"wb");
      return fin;
    }
    static int OpenWriteCheck(const std::string &fname,bool force=false)
    {
      FILE *fi=fopen(fname.c_str(),"rb");
      if (fi!=NULL) {
        fclose(fi);
        if (force || PromptYesNo(std::string("warning: file '" + fname + "' already exists! Overwrite (y/n)?"))) return 1;
        else return 0;
      } else return 1;
    }
    static int SMod(int val,int m) // signed modulo
    {
       val=val<0?(val%m) +m : val%m;
       return val;
    }
    static double Metre_To_MillKm(double len)
    {
       return len/(1000.*1000000.);
    }
    static double SqMetre_To_MillHa(double area)
    {
       return area/(10000.*1000000.);
    }
    static double ToRadians(double deg)
    {
       return (deg*M_PI/180.0);
    }
    static void ToUpper(std::string &str)
    {
      std::transform(str.begin(), str.end(),str.begin(), ::toupper);
    }
    static void Split(std::string &str,std::string &key,std::string &val)
    {
      std::string ws="\r\t\n ";
      std::size_t found = str.find_first_of(ws);
      if (found!=std::string::npos) {
        std::string tkey=str.substr(0,found);
        std::string tval=str.substr(found+1);
        key=Utils::RemoveWhite(tkey);
        val=Utils::RemoveWhite(tval);
      }
    }
static bool SplitDoubleFast(const std::string &str,std::vector<double>&vec)
{
  const char *pstr=(char*)str.c_str();
  const int slen=str.length();
  const char delim=' ';
  int idx=0;
  int start=0;
  int state=0;
  char token[1025];
  for (int i=0;i<slen;i++)
  {
    if (pstr[i]==delim) {
      if (state==1) {
         int len=i-start;
         if (len>1024){cout << "token too long\n";len=1024;};
         strncpy(token,&pstr[start],len);token[len]='\0';
         if (idx<(int)vec.size()) vec[idx++]=atof(token);
      }
      state=0;
    } else {
       if (state==0) start=i;
       state=1;
    };
  }
  return (idx==(int)vec.size());
}
static bool SplitTokenInt(const std::string &str,std::vector<int> &vec)
{
  const char *pstr=(char*)str.c_str();
  const char delim=' ';
  const int  slen=(int)str.length();
  const int  vlen=(int)vec.size();
  char token[1025];
  int tidx=0;
  int sidx=0;
  int vidx=0;
  while (sidx<slen) {
    while (pstr[sidx]==delim && sidx<slen) sidx++;
    while (pstr[sidx]!=delim && sidx<slen) {
        if (tidx==1024) cout << "warning: token too long\n";
        else token[tidx++]=pstr[sidx];
        sidx++;
    };
    if (tidx) {
       token[tidx]='\0';
       if (vidx<vlen) vec[vidx]=atoi(token);
       else cout << "warning: array access out of bounds!\n";
       //if (std::string(token).compare("0.0")==0) cout << vidx << endl;
       vidx++;
       tidx=0;
    };
  }
  return (vidx==vlen);
}
    static std::string RemoveWhite(std::string &in)
    {
       std::string out;
       for (size_t i=0;i<in.length();i++) {
         if (!iswhite(in[i])) out+=in[i];
       }
      return out;
    }
    static bool iswhite(int c)
    {
      return (c==10 || c==13 || c==32 || c==9);
    };
    static void Put16LH(uint8_t *buf,uint32_t val)
    {
      buf[0]=val&0xff;
      buf[1]=(val>>8)&0xff;
    }
    static void Put24LH(uint8_t *buf,uint32_t val)
    {
      buf[0]=val&0xff;
      buf[1]=(val>>8)&0xff;
      buf[2]=(val>>16)&0xff;
    }
    static void Put32LH(uint8_t *buf,uint32_t val)
    {
      buf[0]=val&0xff;
      buf[1]=(val>>8)&0xff;
      buf[2]=(val>>16)&0xff;
      buf[3]=(val>>24)&0xff;
    }
    static void Put48LH(uint8_t *buf,uint64_t val)
    {
      buf[0]=val&0xff;
      buf[1]=(val>>8)&0xff;
      buf[2]=(val>>16)&0xff;
      buf[3]=(val>>24)&0xff;
      buf[4]=(val>>32)&0xff;
      buf[5]=(val>>40)&0xff;
    }
    static void Put64LH(uint8_t *buf,uint64_t val)
    {
      buf[0]=val&0xff;
      buf[1]=(val>>8)&0xff;
      buf[2]=(val>>16)&0xff;
      buf[3]=(val>>24)&0xff;
      buf[4]=(val>>32)&0xff;
      buf[5]=(val>>40)&0xff;
      buf[6]=(val>>48)&0xff;
      buf[7]=(val>>56)&0xff;
    }
    static int Get16LH(uint8_t *buf)
    {
      uint32_t val;
      val=buf[0];
      val+=(buf[1]<<8);
      return val;
    }
    static uint32_t Get24LH(uint8_t *buf)
    {
      uint32_t val;
      val=buf[0];
      val+=(buf[1]<<8);
      val+=(buf[2]<<16);
      return val;
    }
    static uint32_t Get32LH(uint8_t *buf)
    {
      uint32_t val;
      val=buf[0];
      val+=(buf[1]<<8);
      val+=(buf[2]<<16);
      val+=(buf[3]<<24);
      return val;
    }
    static uint64_t Get48LH(uint8_t *buf)
    {
      uint64_t val;
      val=(uint64_t)buf[0];
      val+=((uint64_t)buf[1]<<8);
      val+=((uint64_t)buf[2]<<16);
      val+=((uint64_t)buf[3]<<24);
      val+=((uint64_t)buf[4]<<32);
      val+=((uint64_t)buf[5]<<40);
      return val;
    }
    static uint64_t Get64LH(uint8_t *buf)
    {
      uint64_t val;
      val=(uint64_t)buf[0];
      val+=((uint64_t)buf[1]<<8);
      val+=((uint64_t)buf[2]<<16);
      val+=((uint64_t)buf[3]<<24);
      val+=((uint64_t)buf[4]<<32);
      val+=((uint64_t)buf[5]<<40);
      val+=((uint64_t)buf[6]<<48);
      val+=((uint64_t)buf[7]<<56);
      return val;
    }
    static off64_t GetFileSize(FILE *ifile)
    {
      fseeko64(ifile, 0L, SEEK_END);
      off64_t fsize = ftello64(ifile);
      fseeko64(ifile, 0L, SEEK_SET);
      return fsize;
    }
    static uint64_t GetStreamSize(ifstream &stream)
    {
      stream.seekg (0, stream.end);
      int64_t ssize = stream.tellg();
      stream.seekg (0, stream.beg);
      return ssize;
    }
static void ReplaceExt(const std::string &ifile,std::string &ofile,std::string ext,bool forced=false)
{
  if (!ofile.length() || forced) // replace ext of infile
  {
    size_t found=ifile.find_last_of(".");
    if (found==string::npos) {
      ofile=ifile;
      ofile+=ext;
    } else {
      ofile=ifile.substr(0,found);
      ofile+=ext;
    }
  }
}
static std::string SecToTime(double time)
{
  std::stringstream str_time;
  int itime=round(time);
  int hours=itime/3600;
  itime=itime%3600;
  int minutes=itime/60;
  itime=itime%60;
  int seconds=itime;
  str_time << setw(2) << setfill('0') << hours << ":";
  str_time << setw(2) << setfill('0') << minutes << ":";
  str_time << setw(2) << setfill('0') << seconds;
  return str_time.str();
}
      static std::string ConvertFixed(double number,int decimals)
      {
        std::ostringstream buff;
        buff<<setprecision(decimals)<<fixed<<number;
        return buff.str();
      };
};

class StringUtils {
  public:
    // works for ASCII
    static std::string toupper(const std::string &s)
    {
      std::string ret(s.size(), char());
        for(size_t i = 0; i < s.size(); ++i)
          ret[i] = (s[i] <= 'z' && s[i] >= 'a') ? s[i]-('a'-'A') : s[i];
      return ret;
    }

static void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
};

#endif // GLOBAL_H
