#include "tiff.h"

IMGTIFF::IMGTIFF()
:tif((TIFF*)0),
scanline_size(0),tiles_size(0),tilewidth(0),tileheight(0),nrow(0),
comp(0),bps(0),spp(0),isTiled(false)
{
  TIFFSetWarningHandler(&myTIFFWarningHandler);
}

int IMGTIFF::ReadHeader(const char *fname)
{
  tif=TIFFOpen(fname,"r");
  if (tif) {
    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&width);
    TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height);
    TIFFGetField(tif,TIFFTAG_BITSPERSAMPLE,&bps);
    TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&spp);
    TIFFGetField(tif,TIFFTAG_TILEWIDTH,&tilewidth);
    TIFFGetField(tif,TIFFTAG_TILELENGTH,&tileheight);
    TIFFGetField(tif,TIFFTAG_COMPRESSION,&comp);

    scanline_size=TIFFScanlineSize(tif);
    tiles_size=TIFFTileSize(tif);

    /*uint32_t number_of_stripes=TIFFNumberOfStrips(tif);
    uint32_t number_of_tiles=TIFFNumberOfTiles(tif);*/

    isTiled=TIFFIsTiled(tif);

    return 0;
  } else return 1;
}

void IMGTIFF::PrintInfo()
{
  cout << "TIFF: " << width << "x" << height;
  if (isTiled) cout << ", Tiled (" << tilewidth << "x" << tileheight << ")";
  switch (comp) {
    case 5: cout << ", LZW";break;
    case 8: cout << ", Deflate";break;
    case 32773: cout << ", Packbits";break;
    default: cout << ", Unknown compression";break;
  }
  cout << endl;
  cout << bps << " Bits per Sample, " << spp << " Samples per Pixel" << endl;
}

void IMGTIFF::Close()
{
  TIFFClose(tif);
}

void IMGTIFF::StartReader()
{
  rowbuffer=new uint8_t[width];
  if (isTiled) {
    tiffbuf=_TIFFmalloc(tiles_size);

    tile_buf=new uint8_t*[tileheight];
    for (uint32_t i=0;i<tileheight;i++) tile_buf[i]=new uint8_t[scanline_size];

  } else tiffbuf=_TIFFmalloc(scanline_size);

  nrow=0;
  nrow_mod=0;
}

void IMGTIFF::StopReader()
{
  if (rowbuffer) delete []rowbuffer;
  if (isTiled) {
    for (uint32_t i=0;i<tileheight;i++) delete []tile_buf[i];
    delete []tile_buf;
  }
  _TIFFfree(tiffbuf);
}

int IMGTIFF::ReadRow(uint8_t*buf)
{
  if (!isTiled) {
    int error=TIFFReadScanline(tif,tiffbuf,nrow);
    for (int i=0;i<width;i++) buf[i]=((uint8_t*)tiffbuf)[i];
    nrow++;
    if (error==1) return scanline_size;
    else return 0;
  } else {
    if (nrow_mod==0) {
      for (int ncol=0;ncol<width;ncol+=tilewidth) {
        TIFFReadTile(tif,tiffbuf,ncol,nrow,0,0);

        for (uint32_t j=0;j<tileheight;j++) {
          for (uint32_t i=0;i<tilewidth;i++) {
              if (ncol+i<(uint32_t)width) tile_buf[j][ncol+i]=((uint8_t*)tiffbuf)[j*tilewidth+i];
            }
        }
      }
      nrow+=tileheight;
    }
    for (int i=0;i<width;i++) buf[i]=tile_buf[nrow_mod][i];
    nrow_mod=(nrow_mod+1)%tileheight;
    return 0;
  }
}
