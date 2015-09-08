#include "cm.h"


CM::CM(int width,FILE *file,int decode)
:m_width(width),Coder(file,decode)
{
  linebuf=new uint8_t[m_width];
}

CM::~CM()
{
  if (linebuf) delete []linebuf,linebuf=0;
}

void CM::Init()
{
   for (int i=0;i<8;i++) Counter[i]=StatCounter(32,1<<13);
   Coder.Init();
   rowcounter=0;
}

void CM::Stop()
{
  Coder.Stop();
}

void CM::EncodeRow(uint8_t *buf)
{
  int n,w,nw;
  for (int i=0;i<m_width;i++)
  {
    n=w=nw=0;
    if (rowcounter>0) {
      n=linebuf[i];
      if (i>0) nw=linebuf[i-1];
    }
    if (i>0) w=buf[i-1];

    int ctx4=w+(n<<1)+(nw<<2);
    int bit=buf[i];
    Coder.EncodeBitOne(Counter[ctx4].p1(),bit);
    Counter[ctx4].update(bit);
  }
  memcpy(linebuf,buf,m_width*sizeof(uint8_t));
  rowcounter++;
}

void CM::DecodeRow(uint8_t *buf)
{
  int n,w,nw;
  for (int i=0;i<m_width;i++)
  {
    n=w=nw=0;
    if (rowcounter>0) {
      n=linebuf[i];
      if (i>0) nw=linebuf[i-1];
    }
    if (i>0) w=buf[i-1];

    int ctx4=w+(n<<1)+(nw<<2);
    int bit=Coder.DecodeBitOne(Counter[ctx4].p1());
    Counter[ctx4].update(bit);
    buf[i]=bit;
  }
  memcpy(linebuf,buf,m_width*sizeof(uint8_t));
  rowcounter++;
}
