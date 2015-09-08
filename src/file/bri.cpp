#include "bri.h"

// adds a run to the list of runs
void BRI::PushRun(int idx,int pos,int len)
{
  run trun;
  trun.pos=pos;
  trun.len=len;
  runs[idx]=trun;
}

// determine all 1-runs in a row
int BRI::GetRuns(uint8_t *linebufs,int &flags)
{
  flags=0;
  int idx=0;
  int nrun=0;
  int max_run=0;
  int max_pos=0;
  int last_epos=0;
  for (int i=0;i<width;i++)
  {
     int val=linebufs[i];
     pixel_hist[val]++;
     if (m_globcover)
     {
       if (val==40 || val==160) val=1;
       else val=0;
     } else if (val!=0 && val != 1) {
         //cout << num_warnings << '\n';
         num_invalidpixel++;
         //cout << endl << "invalid pixel value (" << i << ":" << (int)val << ")" << endl;
         val=0; // set to no-forest
     };

     if (val==1) nrun++;
     else if (nrun) {
       int pos=i-nrun;
       int rpos=pos-last_epos;
       PushRun(idx,rpos,nrun);
       if (nrun>max_run)max_run=nrun;
       if (rpos>max_pos)max_pos=rpos;
       idx++;last_epos=pos+nrun;nrun=0;
     };
  }
  if (nrun) {
    int pos=width-nrun;
    int rpos=pos-last_epos;
    PushRun(idx,rpos,nrun);
    if (nrun>max_run)max_run=nrun;
    if (rpos>max_pos)max_pos=rpos;
    idx++;last_epos=pos+nrun;nrun=0;
  };
  if (max_run<(1<<16)) flags=BitSet(flags,1);
  if (max_pos<(1<<16)) flags=BitSet(flags,2);
  return idx;
}

int BRI::CompressRuns(uint8_t *buf,int nruns,int flags)
{
  int idx=0;
  int row_data=nruns;
  bool small_len=BitTest(flags,1);
  bool small_pos=BitTest(flags,2);

  if (small_len) {row_data=BitSet(row_data,31);nlines_smalllen++;};
  if (small_pos) {row_data=BitSet(row_data,30);nlines_smallpos++;};

  Utils::Put32LH(buf+idx,row_data);idx+=4;
  for (int i=0;i<nruns;i++)
  {
      if (small_pos) {Utils::Put16LH(buf+idx,runs[i].pos);idx+=2;}
      else {Utils::Put32LH(buf+idx,runs[i].pos);idx+=4;}

      if (small_len) {Utils::Put16LH(buf+idx,runs[i].len);idx+=2;}
      else {Utils::Put32LH(buf+idx,runs[i].len);idx+=4;}
  }
  return idx;
}

int BRI::DecompressRuns(uint8_t *buf,int nruns,bool small_pos,bool small_runs)
{
  int idx=0;
  int last_pos=0;
  for (int i=0;i<nruns;i++)
  {
     run trun;
     if (small_pos) {trun.pos=last_pos+Utils::Get16LH(buf+idx);idx+=2;}
     else {trun.pos=last_pos+Utils::Get32LH(buf+idx);idx+=4;}

     if (small_runs) {trun.len=Utils::Get16LH(buf+idx);idx+=2;}
     else {trun.len=Utils::Get32LH(buf+idx);idx+=4;};
     runs[i]=trun;
     last_pos=trun.pos+trun.len;
  }
  return idx;
}

void BRI::PrintInfo()
{
  cout << endl << "BRI";
  cout << ": " << width << "x" << height << endl;
  cout << "filesize: " << (fsize>>20) << " MB, data start: " << hdrpos << " bytes" << endl;
}

int BRI::Create(std::string &str_ofile)
{
  file=fopen(str_ofile.c_str(),"wb");
  if (file)
  {
    uint8_t buf[16];
    Utils::Put32LH(buf,99991); // magic number
    Utils::Put32LH(buf+4,width); // width
    Utils::Put32LH(buf+8,height); // height
    Utils::Put32LH(buf+12,0); // flags
    fwrite(buf,1,16,file);
    return 0;
  } else return 1;
}

int BRI::ReadHeader()
{
  fsize=Utils::GetFileSize(file);
  uint8_t buf[16];
  fread(buf,1,16,file);
  int magic=Utils::Get32LH(buf);
  if (magic==99991)
  {
    width=Utils::Get32LH(buf+4);
    height=Utils::Get32LH(buf+8);
    flags=Utils::Get32LH(buf+12);
    hdrpos=ftello64(file);
    return 0;
  } else {
    fseeko64(file,0,SEEK_SET);
    return 1;
  }
}

void BRI::SeekStart()
{
  fseeko64(file,hdrpos,SEEK_SET);
}

void BRI::Close()
{
  fclose(file);
}

int BRI::WriteRow(uint8_t *linebuf)
{
  int flags;
  int nruns=GetRuns(linebuf,flags);
  int maxcompsize=(nruns*8)+4;
  if (maxcompsize>=10*width) { // we get no compression
    cout << "error: negative compression, buffer overflow\n";
    return 0;
  } else {
    //if (maxcompsize>=width) cout << "error: cannot compress: " << maxcompsize << ">" << width << endl;
    int ncompressed=CompressRuns(tbuf,nruns,flags);
    int nwritten=fwrite(tbuf,1,ncompressed,file);
    if (nwritten!=ncompressed) {cout << "error: can not write" << endl;return 0;};
    return nwritten;
  }

}

void BRI::Start()
{
  tbuf=new uint8_t[10*width];
  runs.resize(width);
}

void BRI::Stop()
{
  if (tbuf) delete []tbuf,tbuf=0;
}

void BRI::PrintProgress(int y,int maxline,FILE *file)
{
  cout << (y) << "/" << maxline << " (" << (ftello64(file)>>20) << " MB)\r";
}

void BRI::ConvertFromPGM(PGM &myPGM,std::string &str_ofile,bool globcover)
{
  m_globcover=globcover;
  nlines_smalllen=nlines_smallpos=0;
  SetWidth(myPGM.GetWidth());
  SetHeight(myPGM.GetHeight());
  if (Create(str_ofile)==0)
  {
     Start();
     myPGM.StartReading();

     uint64_t total_written=0;
     int lines_to_code=myPGM.GetHeight();

     //CM myCM(width,file);
     //myCM.Init();

     for (int y=0;y<lines_to_code;y++)
     {
       if(!myPGM.ReadRow()) return;
       #if 0
       myCM.EncodeRow(myPGM.linebuf);
       #else
       total_written+=WriteRow(myPGM.linebuf);
       #endif
       if ((y+1)%100==0) PrintProgress(y+1,height,file);
     }
     PrintProgress(height,height,file);
     if (myPGM.GetHeight()) {
       double sp=(double)nlines_smallpos/(double)myPGM.GetHeight()*100.0;
       double sl=(double)nlines_smalllen/(double)myPGM.GetHeight()*100.0;
       cout << endl << " 16-bit pos: " << sp << ", 16-bit len: " << sl;
     }
     myPGM.StopReading();
     //myCM.Stop();
     Stop();
     cout << "\ntotal written: " << ftello64(file) << " bytes\n";
     if (num_invalidpixel) cerr << "  number of invalid pixels: " << Utils::ConvertFixed((num_invalidpixel*100.)/((double)myPGM.GetWidth()*lines_to_code),2)<<"%\n";
     for (int i=0;i<256;i++) {
        if (pixel_hist[i]) cout << std::setw(4) << i << ": " << Utils::ConvertFixed((pixel_hist[i]*100.)/((double)myPGM.GetWidth()*lines_to_code),2)<<"%\n";
     };
     Close();
  }
}

void BRI::ConvertFrom(ASC &myASC,std::string &str_ofile)
{
  nlines_smalllen=nlines_smallpos=0;

  int lines_to_code=myASC.getExtendHeight();
  SetWidth(myASC.getExtendWidth());
  SetHeight(lines_to_code);
  if (Create(str_ofile)==0)
  {
     uint64_t total_written=0;
     Start();
     myASC.StartReading();
     for (int y=0;y<lines_to_code;y++)
     {
       myASC.ReadRow();
       if ((y+1)%100==0) PrintProgress(y+1,lines_to_code,file);
       total_written+=WriteRow(myASC.linebuf);
     }
     PrintProgress(lines_to_code,lines_to_code,file);
     myASC.StopReading();
     Stop();
     cout << "\ntotal written: " << ftello64(file) << " bytes\n";
     Close();
  }

}


void BRI::ConvertToBytes(uint8_t *buf,int nruns)
{
  memset(buf,0,width*sizeof(uint8_t));
  for (int i=0;i<nruns;i++)
  {
    for (int k=0;k<runs[i].len;k++) buf[runs[i].pos+k]=1;
  }
}

int BRI::ReadRow(uint8_t *buf)
{
  int nread=fread(tbuf,1,4,file);
  uint32_t row_data=Utils::Get32LH(tbuf);
  bool small_len=BitTest(row_data,31);
  bool small_pos=BitTest(row_data,30);

  row_data=BitClear(row_data,31);
  row_data=BitClear(row_data,30);

  int nruns=row_data;

  int bytes_to_read=0;
  if (small_pos) bytes_to_read+=2*nruns;
  else bytes_to_read+=4*nruns;

  if (small_len) bytes_to_read+=2*nruns;
  else bytes_to_read+=4*nruns;

  nread+=fread(tbuf,1,bytes_to_read,file);
  if (nread!=bytes_to_read+4) {cout << endl << "could not read!" << endl;return 1;};
  DecompressRuns(tbuf,nruns,small_pos,small_len);
  ConvertToBytes(buf,nruns);
  return nread;
}

void BRI::ConvertToPGM(PGM &myPGM,std::string &str_ofile)
{
  myPGM.SetWidth(width);
  myPGM.SetHeight(height);
  if (myPGM.Create(str_ofile)==0)
  {
    Start();
    //CM myCM(width,file,1);
    //myCM.Init();

    myPGM.StartReading();
    uint64_t total_read=0;
    //height=10000;
    for (int y=0;y<height;y++)
    {
      //myCM.DecodeRow(myPGM.linebuf);
      total_read+=ReadRow(myPGM.linebuf);
      myPGM.WriteRow();
      if ((y+1)%100==0) PrintProgress(y+1,height,file);
      //if (y%100==0) cout << (y+1) << "/" << height << " (" << (total_read>>20) << " MB)\r";
    }
    PrintProgress(height,height,file);
    myPGM.StopReading();
    //myCM.Stop();
    Stop();
    myPGM.Close();
  }
}
