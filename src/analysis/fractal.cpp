#include "fractal.h"
#include "..\utils.h"

FractalDimension::FractalDimension(IMG &image)
:myIMG(image),box_count(0)
{

}

void FractalDimension::PrintProgress(int row,int nrows)
{
  std::cout << "Reading: " << Utils::ConvertFixed(row/float(nrows)*100.,1) << "%\r";
}

bool FractalDimension::AnalyseBlock(int nblock,int box_size)
{
  for (int y=0;y<box_size;y++)
    for (int x=0;x<box_size;x++) {
         int pixel=row_block[y][(nblock*box_size)+x];
         if (pixel) return true;
    }
  return false;
}

void FractalDimension::AnalyseRowBlock(int nblocksx,int box_size)
{
  for (int nblock=0;nblock<nblocksx;nblock++) {
    if (AnalyseBlock(nblock,box_size)) box_count++;
  }
}

void FractalDimension::AnalyseBoxSize(int box_size)
{
  box_count=0;

  myIMG.StartReader();

  row_block.resize(box_size,std::vector<uint8_t>(myIMG.GetWidth(),0));
  std::cout << "mem usage: " << ((uint64_t(box_size)*uint64_t(myIMG.GetWidth()))>>20) << " mb" << std::endl;

  int nblocksy=myIMG.GetHeight()/box_size;
  int nblocksx=myIMG.GetWidth()/box_size;
  int64_t num_boxes=(int64_t)nblocksx*(int64_t)nblocksy;
  std::cout << "boxsize: " << box_size << " pixels (" << num_boxes  << " boxes)" << std::endl;
  for (int nblock=0;nblock<nblocksy;nblock++)
  {
     for (int blockrow=0;blockrow<box_size;blockrow++) {
      int nbytes=myIMG.ReadRow();
      for (int i=0;i<myIMG.GetWidth();i++) row_block[blockrow][i]=myIMG.rowbuffer[i];
     }
     AnalyseRowBlock(nblocksx,box_size);
     PrintProgress(nblock+1,nblocksy);
  }
   PrintProgress(nblocksy,nblocksy);std::cout << std::endl;

  std::cout << "box count:    " << box_count << std::endl;

  myIMG.StopReader();
}


