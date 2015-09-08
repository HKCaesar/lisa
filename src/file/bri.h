#ifndef BRI_H
#define BRI_H

#include "pgm.h"
#include "asc.h"
//#include "model/cm.h"
#include "imgfile.h"

/*#include "model/rangecoder.h"
#include "model/counter.h"
#include "model/mixer.h"*/

/*
  BRI Binary Run Image format
  represents binary images for arbitrary sizes

  header format (size: 12 byte)
  4b: magic number: 99991
  4b: width
  4b: height
  4b: flags

  image is composed into rows 1..height
  every row starts with
  4b: number of runs, bit 31 is set if len of runs<(1<<16), bit 30 is set if pos of runs<(1<<16)
  (number of runs)*(size of pos+size of len) row-size in bytes (pairs of pos,len)
*/

struct run
{
  int pos,len;
};

class BRI : public IMGFile
{
  public:
    BRI():m_globcover(false),num_invalidpixel(0)
    {
      for (int i=0;i<256;i++) pixel_hist[i]=0;
    };
    void ConvertFromPGM(PGM &myPGM,std::string &str_ofile,bool globcover=false);
    void ConvertToPGM(PGM &myPGM,std::string &str_ofile);
    void ConvertFrom(ASC &myASC,std::string &str_ofile);
    int ReadHeader();
    void PrintInfo();
    void Start();
    void Stop();
    void SeekStart();
    int Create(std::string &str_ofile);
    void Close();
    int WriteRow(uint8_t *linebuf);
    int ReadRow(uint8_t *buf);
  protected:
    void PrintProgress(int line,int maxline,FILE *file);
    void PushRun(int idx,int pos,int len);
    void ConvertToBytes(uint8_t *buf,int nruns);
    int GetRuns(uint8_t *linebuf,int &flags);
    int CompressRuns(uint8_t *buf,int nruns,int flags);
    int DecompressRuns(uint8_t *buf,int nruns,bool small_pos,bool small_runs);
    std::vector <run> runs;
    off64_t hdrpos;
    uint8_t* tbuf;
    int flags;
    bool m_globcover;
    int nlines_smallpos,nlines_smalllen;
    uint64_t num_invalidpixel;
    uint64_t pixel_hist[256];
};

#endif // BRI_H
