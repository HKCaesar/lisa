#include "shape.h"

void ShapeFile::GeoRef(int rx,int ry,int &xpos,int &ypos)
{
  xpos=ypos=0;
  double ref_long=ref_left+(double)rx*ref_cellsize;
  double ref_lat=ref_top-(double)ry*ref_cellsize;

  xpos=round((ref_long-left)/cellsize);
  ypos=round((top-ref_lat)/cellsize);
}

bool ShapeFile::IsCLASS(int rx,int ry)
{
  if (referenced && vec_classes.size()) {
    int x,y;
    GeoRef(rx,ry,x,y);
    if (x < 0 || x>=width) return false;
    if (y < 0 || y>=height) return false;
    int val=data[y*width+x];
    for (auto c : vec_classes)
      if (c==val) return true;
    return false;
  } else return true; // accept every pixel if we have no shape
}

void ShapeFile::ReadFile()
{
    // stats
    int64_t num_pixel_valid,num_pixel_total;
    num_pixel_valid=num_pixel_total=0;
    int64_t val_total=0.;
    int min_value,max_value;
    max_value=std::numeric_limits<int>::min();
    min_value=std::numeric_limits<int>::max();

    std::vector <int>vecline;
    vecline.resize(width);

    linenum=0;
    while (linenum<height)
    {
        if (!getLine(line)) cout << "error reading file at line: "<<linenum << endl;
        else if (!SplitTokens<int>(line,vecline)) cout << "error parsing file at line: " << linenum << endl;
        else {
          for (size_t i=0;i<vecline.size();i++) {
            num_pixel_total++;
            int val=vecline[i];
            data[linenum*width+i]=val;

            if (val!=nodata_value) {
              if (val<min_value) min_value=val;
              else if (val>max_value) max_value=val;

              num_pixel_valid++;val_total+=val;
            }
          }
        }
        linenum++;
        if (linenum%100==0) cout << "Reading: " << Utils::ConvertFixed(linenum*100/(double)height,1) << "%\r";
    }
    cout << "Reading: " << Utils::ConvertFixed(linenum*100/(double)height,1) << "%\n";
    cout << "total pixels: " << num_pixel_total << " ("<<width<<"x"<<height<<"="<<width*height<<")\n";
    cout << "valid pixels: " << num_pixel_valid << " ("<<Utils::ConvertFixed(num_pixel_valid*100./(double)num_pixel_total,1)<<"%)\n";
    cout << "mean: " << (val_total/(double)num_pixel_valid) << " t/ha\n";
    cout << "min : " << min_value << "\n";
    cout << "max : " << max_value << "\n";
}

bool ShapeFile::ReadShapeFile(const std::string &fname)
{
  if (OpenRead(fname)) {
    cout << "reading rasterfile: '" << fname << "'" << endl;
    if (ReadHeader()==0) {
      PrintInfo();
      data=new int[width*height];
      ReadFile();
      fclose(file);
      return true;
    } else return false;
  } else return false;
}
