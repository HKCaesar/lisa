#include "bri.h"

void IMGBRI::Start()
{
  tbuf=new uint8_t[10*width];
  tbuf2=new uint8_t[10*width];
  tbuf3=new uint8_t[10*width];
  rowbuffer=new uint8_t[width];
  runs.resize(width);
  mySIC = new SIC(width,comptype);
}

void IMGBRI::Stop()
{
  if (tbuf) delete []tbuf,tbuf=0;
  if (tbuf2) delete []tbuf2,tbuf2=0;
  if (tbuf3) delete []tbuf3,tbuf3=0;
  if (rowbuffer) delete []rowbuffer,rowbuffer=0;
  if (mySIC) delete mySIC;
}

void IMGBRI::PrintInfo()
{
  cout << endl << "BRI";
  cout << ": " << width << "x" << height << ", flags: " << flags << ", comp: ";
  switch (comptype) {
    case SIC::COMP_NONE:cout << "none";break;
    case SIC::COMP_BILEVEL:cout << "bilevel";break;
    case SIC::COMP_GRAY:cout << "gray level";break;
    default: cout << "unknown";break;
  }
  cout << endl;
  cout << "filesize: " << (fsize>>20) << " MiB, data start: " << hdrpos << " bytes" << endl;
}

int IMGBRI::Create(std::string &str_ofile)
{
  file=fopen(str_ofile.c_str(),"wb");
  if (file)
  {
    uint8_t buf[16];
    Utils::Put32LH(buf,99991); // magic number
    Utils::Put32LH(buf+4,width); // width
    Utils::Put32LH(buf+8,height); // height
    int flags=0;
    int t;
    switch (comptype) {
      case SIC::COMP_NONE:t=0;break;
      case SIC::COMP_BILEVEL:t=1;break;
      case SIC::COMP_GRAY:t=2;break;
      default: t=0;break;
    }
    flags|=t;
    Utils::Put32LH(buf+12,flags); // flags
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
    switch (flags&15) {
      case 0:comptype=SIC::COMP_NONE;break;
      case 1:comptype=SIC::COMP_BILEVEL;break;
      case 2:comptype=SIC::COMP_GRAY;break;
    }
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
  if (comptype==SIC::COMP_NONE) { // obsolete
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
  } else if (comptype==SIC::COMP_BILEVEL) {
    nread=fread(tbuf,1,4,file);
    uint32_t bytes_to_read=Utils::Get32LH(tbuf);
    nread+=fread(tbuf+4,1,bytes_to_read,file);
    int comp_bufpos=mySIC->DecompressRowBinary(tbuf,buf);
    if (nread!=comp_bufpos) cout << " error: invalid decompression: " << nread << " != " << comp_bufpos << endl;
  }
  return nread;
}

int IMGBRI::WriteRow(uint8_t *linebuf)
{
  int ncompressed=0;
  if (comptype==SIC::COMP_BILEVEL) {
       ncompressed=mySIC->CompressRowBinary(linebuf,tbuf);
       mySIC->DecompressRowBinary(tbuf,tbuf2);
       int nerr=0;
       for (int i=0;i<width;i++) if (linebuf[i]!=tbuf2[i]) nerr++;
       if (nerr) cout << nerr << " errors in decompression\n";
  } else if (comptype==SIC::COMP_GRAY) {
    ncompressed=mySIC->CompressRowGrey(linebuf,tbuf);
    /*mySIC->DecompressRowGrey(tbuf,tbuf2);
    for (int i=0;i<width;i++) {
        if (linebuf[i]!=tbuf2[i]) cout << "error in decompression at pos: " << i << endl;
    }*/
  }
  int nwritten=fwrite(tbuf,1,ncompressed,file);
  if  (nwritten!=ncompressed) {cout << "error: can not write" << endl;return 0;};
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

void IMGBRI::ConvertToBRI(IMG &myIMG,std::string &str_outfile,SIC::COMP_TYPE compression_type)
{
  comptype=compression_type;

  SetWidth(myIMG.GetWidth());
  SetHeight(myIMG.GetHeight());
  if (Create(str_outfile)==0) {
    uint8_t *rowbuffer=new uint8_t[myIMG.GetWidth()];
    Start();
    myIMG.StartReader();
    uint64_t total_written=0;
    for (int y=0;y<myIMG.GetHeight();y++)
    {
      myIMG.ReadRow(rowbuffer);
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

