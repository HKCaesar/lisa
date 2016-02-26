#include "bri.h"

void IMGBRI::Start()
{
  tbuf=new uint8_t[10*width];
  tbuf2=new uint8_t[10*width];
  tbuf3=new uint8_t[10*width];
  rowbuffer=new uint8_t[width];
  runs.resize(width);
}

void IMGBRI::Stop()
{
  if (tbuf) delete []tbuf,tbuf=0;
  if (tbuf2) delete []tbuf2,tbuf2=0;
  if (tbuf3) delete []tbuf3,tbuf3=0;
  if (rowbuffer) delete []rowbuffer,rowbuffer=0;
}

void IMGBRI::PrintInfo()
{
  cout << endl << "BRI";
  cout << ": " << width << "x" << height << ", flags: " << flags << endl;
  cout << "filesize: " << (fsize>>20) << " MiB, data start: " << hdrpos << " bytes" << endl;
}

int IMGBRI::Create(std::string &str_ofile,int flags)
{
  file=fopen(str_ofile.c_str(),"wb");
  if (file)
  {
    uint8_t buf[16];
    Utils::Put32LH(buf,99991); // magic number
    Utils::Put32LH(buf+4,width); // width
    Utils::Put32LH(buf+8,height); // height
    Utils::Put32LH(buf+12,flags); // flags
    comp_type=flags&15; // set compression type
    fwrite(buf,1,16,file);
    return 0;
  } else return 1;
}

int IMGBRI::ReadHeader()
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
    comp_type=flags&15;
    hdrpos=ftello64(file);
    return 0;
  } else {
    fseeko64(file,0,SEEK_SET);
    return 1;
  }
}

// Compression type - 0
// runs of 0/1 are bitcompressed into fixed-size (16 or 32-bit) tokens
// obsolete

/*
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
}*/

int IMGBRI::DecompressRuns(uint8_t *buf,int nruns,bool small_pos,bool small_runs)
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


// Compression type - 1
// runs of 0/1 are coded via VLE RiceCodes

int IMGBRI::CompressRunsVLE(uint8_t *linebuf,uint8_t *outbuf)
{
  BitBuffer bitout(outbuf+4);
  bool sbit=linebuf[0];
  bitout.PutBits(sbit,1);
  int nrun=0;
  //int N=1,A=256;
  for (int i=1;i<width;i++)
  {
     if (linebuf[i]==sbit) nrun++;
     else {
        //bitout.PutRice(nrun,BitBuffer::EstimateK(N,A));
        bitout.PutEliasGamma(nrun+1);
        //if (N>=256) {N>>=1;A>>=1;};
        //N++;A+=nrun;
        nrun=0;
        sbit=!sbit; // flip bit
     }
  }
  //bitout.PutRice(nrun,BitBuffer::EstimateK(N,A));
  bitout.PutEliasGamma(nrun+1);
  bitout.Flush();
  Utils::Put32LH(outbuf,bitout.GetBytesProcessed());
  return bitout.GetBytesProcessed()+4;
}

int IMGBRI::DecompressRunsVLE(uint8_t *inbuf,uint8_t *linebuf)
{
  BitBuffer bitin(inbuf+4);
  bool sbit=bitin.GetBits(1);

  //int N=1,A=256;
  int i=0;
  while (i<width) {
    //int nrun=bitin.GetRice(BitBuffer::EstimateK(N,A));
    int nrun=bitin.GetEliasGamma()-1;
    for (int r=0;r<(nrun+1);r++) linebuf[i++]=sbit;
    //if (N>=256) {N>>=1;A>>=1;};
    //N++;A+=nrun;
    sbit=!sbit;
  }
  return bitin.GetBytesProcessed()+4;
}


void IMGBRI::SeekStart()
{
  fseeko64(file,hdrpos,SEEK_SET);
}

void IMGBRI::Close()
{
  fclose(file);
}

void IMGBRI::ConvertToBytes(uint8_t *buf,int nruns)
{
  memset(buf,0,width*sizeof(uint8_t));
  for (int i=0;i<nruns;i++)
  {
    const int rlen=runs[i].len;
    const int rpos=runs[i].pos;
    for (int k=0;k<rlen;k++) buf[rpos+k]=1;
  }
}

int IMGBRI::ReadRow(uint8_t *buf)
{
  int nread=0;
  if (comp_type==0) { // obsolete
    nread=fread(tbuf,1,4,file);
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
  } else if (comp_type==1) {
    nread=fread(tbuf,1,4,file);
    uint32_t bytes_to_read=Utils::Get32LH(tbuf);
    nread+=fread(tbuf+4,1,bytes_to_read,file);
    int comp_bufpos=DecompressRunsVLE(tbuf,buf);
    if (nread!=comp_bufpos) cout << " error: invalid decompression: " << nread << " != " << comp_bufpos << endl;
  }
  return nread;
}

//#define LISA_DEBUG

int IMGBRI::WriteRow(uint8_t *linebuf)
{
   int ncompressed=CompressRunsVLE(linebuf,tbuf);
   #ifdef LISA_DEBUG
     int dcompressed=DecompressRunsVLE(tbuf,tbuf2);
     int nerror=0;
     for (int i=0;i<width;i++)
       if (tbuf2[i]!=linebuf[i]) nerror++;
     if (nerror) cout << nerror << " number of errors in line" << endl;
     if (dcompressed!=ncompressed) cout << " decompression error: " << ncompressed << ":" << dcompressed << endl;
    #endif
    int nwritten=fwrite(tbuf,1,ncompressed,file);
    if (nwritten!=ncompressed) {cout << "error: can not write" << endl;return 0;};
    return nwritten;
}


void IMGBRI::PrintProgress(int y,int maxline,FILE *file)
{
  cout << (y) << "/" << maxline << " (" << (ftello64(file)>>20) << " MB)\r";
}

void IMGBRI::PrintCompression(uint64_t insize,FILE *infile,FILE *outfile)
{
  uint64_t inpos=ftello64(infile);
  uint64_t outpos=ftello64(outfile);
  double ratio1=0.,ratio2=0.;
  if (insize) ratio1 =  inpos*100. / (double)insize;
  if (inpos) ratio2 = outpos*100. / (double)inpos;
  cout << (inpos>>20) << " MB (" << Utils::ConvertFixed(ratio1,1) << "%)->" << (outpos>>20) << " MB: " << Utils::ConvertFixed(ratio2,1) << "% \r";
}

void IMGBRI::ConvertToBRI(IMG &myIMG,std::string &str_outfile)
{
  SetWidth(myIMG.GetWidth());
  SetHeight(myIMG.GetHeight());
  if (Create(str_outfile,1)==0) {
    uint8_t *rowbuffer=new uint8_t[myIMG.GetWidth()];
    Start();
    myIMG.StartReader();
    uint64_t total_written=0;
    for (int y=0;y<myIMG.GetHeight();y++)
    {
      myIMG.ReadRow(rowbuffer);
      for (int i=0;i<width;i++) {
        if (rowbuffer[i]!=0 && rowbuffer[i]!=1) {cout << "  warning: invalid input value: " << rowbuffer[i];delete []rowbuffer;return;};
      }
      if ((y+1)%100==0) PrintProgress(y+1,myIMG.GetHeight(),file);
      total_written+=WriteRow(rowbuffer);
    }
    PrintProgress(myIMG.GetHeight(),myIMG.GetHeight(),file);
    myIMG.StopReader();
    Stop();
    cout << "\ntotal written: " << ftello64(file) << " bytes\n";
    Close();
    delete []rowbuffer;
  }
}

void IMGBRI::ConvertToPGM(IMGPGM &myPGM,std::string &str_ofile)
{
  myPGM.SetWidth(width);
  myPGM.SetHeight(height);
  if (myPGM.Create(str_ofile)==0)
  {
    Start();
    //CM myCM(width,file,1);
    //myCM.Init();

    myPGM.StartReader();
    uint64_t total_read=0;
    //height=10000;
    for (int y=0;y<height;y++)
    {
      //myCM.DecodeRow(myPGM.linebuf);
      total_read+=ReadRow(myPGM.rowbuffer);
      myPGM.WriteRow();
      if ((y+1)%100==0) PrintProgress(y+1,height,file);
      //if (y%100==0) cout << (y+1) << "/" << height << " (" << (total_read>>20) << " MB)\r";
    }
    PrintProgress(height,height,file);
    myPGM.StopReader();
    //myCM.Stop();
    Stop();
    myPGM.Close();
  }
}

