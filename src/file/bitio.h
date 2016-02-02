#ifndef BITIO
#define BITIO

// fast log2 on x86
static inline uint32_t ilog2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}

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
      void PutBits(int val,int nbits) { // put nbits going from low to high
        for (int k=0;k<nbits;k++) PutBit( (val>>k)&1);
      }
      int GetBits(int nbits) { // get nbits going from low to high
        int val=0;
        for (int k=0;k<nbits;k++) val|=(GetBit()<<k);
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
      int GetEliasGamma()
      {
         int q=0;
         while (GetBit()) q++;
         return (1<<q)+GetBits(q);
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

#endif // BITIO
