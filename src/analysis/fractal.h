#ifndef FRACTAL_H
#define FRACTAL_H

#include "..\file\img.h"
#include "..\model\sic.h"

class FractalDimension {
  public:
    FractalDimension(IMG &image);
    void Analyse(const std::string &fname);
    void AnalyseBoxSize(int box_size);
    int64_t GetBoxCount(){return box_count;};
  private:
    void PrintProgress(int row,int nrows);
    bool AnalyseBlock(int nblock,int box_size);
    void AnalyseRowBlock(int nblocksx,int box_size);
    std::vector<std::vector<uint8_t>> row_block;
    IMG &myIMG;
    int64_t box_count;
};
#endif // FRACTAL_H
