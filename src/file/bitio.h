#ifndef BITIO
#define BITIO

#if defined(__GNUC__) && defined(__x86_64__)
// fast log2 on x86
static inline uint32_t ilog2(uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}
#else
static inline uint32_t ilog2(uint32_t x) {
  uint32_t y=0;
  while (x>>=1) ++y;
  return y;
}
#endif

// 64-bit bitbuffer, fills bitbuf frow Low to High
class BitBuffer64LH {
  public:
      BitBuffer64LH(uint8_t *buffer)
      :pbuf((uint64_t*)buffer)
      {
        bitbuf=bitcnt=bytes_processed=0;
      }
      void PutBits(int val,int n)
      {
        bitbuf=bitbuf|(uint64_t(val) << bitcnt);
        bitcnt+=n;
        if (bitcnt>=64) {
            *pbuf++=bitbuf; // flush 4-bytes
            bitcnt-=64;
            bitbuf=val>>(n-bitcnt);
            bytes_processed+=8;
        }
      }
      void Flush() { // we waste at most 7-bytes
        if (bitcnt) {*pbuf++=bitbuf;bytes_processed+=8;};
      }
      int GetBits(int n) { //GetBits will read over eof
        uint64_t val=bitbuf>>(64-bitcnt);
        if (bitcnt<=n) {
          bitbuf=*pbuf++; // refill the buffer
          val|=bitbuf<<bitcnt; // add remaining bits
          bitcnt+=64;
          bytes_processed+=8;
        }
        bitcnt-=n;
        return val&((1<<n)-1); // return masked val*/
      }
      void PutEliasGamma(int val)
      {
         int q=ilog2(val);
         PutBits(0,q);
         PutBits(1,1);
         if (q) PutBits(val&((1<<q)-1),q);
      }
      int GetEliasGamma()
      {
         int q=0,r=0;
         while (GetBits(1)==0) q++;
         if (q) r=GetBits(q);
         return (1<<q)+r;
      }
      int GetBytesProcessed(){return bytes_processed;};
  private:
    uint64_t *pbuf;
    uint64_t bitbuf;
    int bitcnt,bytes_processed;
};

// 8-bit bitbuffer, fills bitbuf from high to low
class BitBuffer {
  public:
      BitBuffer(uint8_t *buffer)
      :pbuf(buffer)
      {
        bitbuf=bitcnt=bytes_processed=0;
      }
      inline void PutBit(int bit) {
        bitcnt++;
        bitbuf=bitbuf|(bit<<(8-bitcnt));
        if (bitcnt==8) {
           *pbuf++=bitbuf;
           bitcnt=bitbuf=0;
           bytes_processed++;
        }
      }
      inline int GetBit() {
        if (bitcnt==0) {
          bitbuf=*pbuf++;
          bitcnt=8;
          bytes_processed++;
        }
        bitcnt--;
        return (bitbuf>>bitcnt)&1;
      }
      void PutBits(int64_t val,int nbits) { // put nbits going from low to high
        for (int k=0;k<nbits;k++) PutBit((val>>k)&1);
      }
      int64_t GetBits(int nbits) { // get nbits going from low to high
        int64_t val=0;
        for (int k=0;k<nbits;k++) val|=(int64_t)GetBit()<<k;
        return val;
      }
      static int EstimateK(int N,int A) { // doesn't handle overflows correctly
        int k;
        for (k=0;(N<<k)<A;k++);
        return k;
      }
      void PutEliasGamma(int val)
      {
         int q=ilog2(val);
         for (int i=0;i<q;i++) PutBit(1);
         PutBit(0);
         if (q) PutBits(val&((1<<q)-1),q);
      }
      void PutEliasGamma(int64_t val)
      {
         int64_t x=val;
         int q=0;
         while (x>>=1) ++q;

         for (int i=0;i<q;i++) PutBit(1);
         PutBit(0);
         if (q) PutBits(val&(((int64_t)1<<q)-1),q);
      }
      int64_t GetEliasGamma()
      {
         int64_t q=0;
         while (GetBit()) q++;
         return ((int64_t)1<<q)+GetBits(q);
      }
      void PutRice(int val,int k) { // write varible length rice code with param k<32
         int q=val>>k;
         for (int i=0;i<q;i++) PutBit(1);
         PutBit(0);
         if (k) PutBits(val&((1<<k)-1),k);
      }
      int GetRice(int k)
      {
         int q=0,r=0;
         while (GetBit()) q++;
         if (k) r=GetBits(k);
         return (q<<k)+r;
      }
      void PutExpGolomb(int val,int k) {
        while (val >= (1<<k)) {
            PutBit(1);
            val-=(1<<k);
            k++;
        }
        PutBit(0);
        if (k) PutBits(val &((1<<k)-1),k);
      }
      int GetBytesProcessed(){return bytes_processed;};
      void Flush() {
        if (bitcnt) {*pbuf++=bitbuf;bytes_processed++;};
      }
  private:
    int bitbuf,bitcnt,bytes_processed;
    uint8_t *pbuf;
};

// 8-bit bitbuffer, fills bitbuf from high to low
class BitBufferSafe {
  public:
      BitBufferSafe(vector <uint8_t>&buffer)
      :buffer_(buffer)
      {
        bitbuf=bitcnt=bytes_processed=0;
      }
      inline void StoreByte(uint8_t val)
      {
        if (buffer_.size()<=bytes_processed)  {
          cerr << "bitbuf: resize " << ((buffer_.size()+1)*2) << endl;
          buffer_.resize((buffer_.size()+1)*2);
        }
        buffer_[bytes_processed++]=val;
      }
      inline void PutBit(int bit) {
        bitcnt++;
        bitbuf=bitbuf|(bit<<(8-bitcnt));
        if (bitcnt==8) {
           StoreByte(bitbuf);
           bitcnt=bitbuf=0;
        }
      }
      inline int GetBit() {
        if (bitcnt==0) {
          if (buffer_.size()<=bytes_processed) {
            cerr << "bitbuf: attempt to read over end of input buffer\n";
            bitbuf=0;
          }
          else bitbuf=buffer_[bytes_processed++];
          bitcnt=8;
        }
        bitcnt--;
        return (bitbuf>>bitcnt)&1;
      }
      void PutBits(int64_t val,int nbits) { // put nbits going from low to high
        for (int k=0;k<nbits;k++) PutBit((val>>k)&1);
      }
      int64_t GetBits(int nbits) { // get nbits going from low to high
        int64_t val=0;
        for (int k=0;k<nbits;k++) val|=(int64_t)GetBit()<<k;
        return val;
      }
      static int EstimateK(int N,int A) { // doesn't handle overflows correctly
        int k;
        for (k=0;(N<<k)<A;k++);
        return k;
      }
      void PutEliasGamma(int val)
      {
         int q=ilog2(val);
         for (int i=0;i<q;i++) PutBit(1);
         PutBit(0);
         if (q) PutBits(val&((1<<q)-1),q);
      }
      void PutEliasGamma(int64_t val)
      {
         int64_t x=val;
         int q=0;
         while (x>>=1) ++q;

         for (int i=0;i<q;i++) PutBit(1);
         PutBit(0);
         if (q) PutBits(val&(((int64_t)1<<q)-1),q);
      }
      int64_t GetEliasGamma()
      {
         int64_t q=0;
         while (GetBit()) q++;
         return ((int64_t)1<<q)+GetBits(q);
      }
      void PutRice(int val,int k) { // write varible length rice code with param k<32
         int q=val>>k;
         for (int i=0;i<q;i++) PutBit(1);
         PutBit(0);
         if (k) PutBits(val&((1<<k)-1),k);
      }
      int GetRice(int k)
      {
         int q=0,r=0;
         while (GetBit()) q++;
         if (k) r=GetBits(k);
         return (q<<k)+r;
      }
      void PutExpGolomb(int val,int k) {
        while (val >= (1<<k)) {
            PutBit(1);
            val-=(1<<k);
            k++;
        }
        PutBit(0);
        if (k) PutBits(val &((1<<k)-1),k);
      }
      int GetBytesProcessed(){return bytes_processed;};
      void Flush() {
        if (bitcnt) StoreByte(bitbuf);
      }
  private:
    int bitbuf,bitcnt;
    uint32_t bytes_processed;
    vector <uint8_t> &buffer_;
};


#endif // BITIO
