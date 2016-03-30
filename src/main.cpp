#include "global.h"
#include "common/cmdoptions.h"
#include "cmdline.h"
#include "analysis/map.h"

void PrintVersion(int mode=0)
{
  cout << LISA_VERSION << endl;
  if (mode) {
  cout << "compiled ";
  #ifdef __GNUC__
    cout << "with GNU " << __VERSION__;
  #endif
  #if __x86_64__
    cout << " (64-bit)";
  #else
    cout << " (32-bit)";
  #endif
  cout << " on " << __DATE__ << endl;
  }
  #ifdef TIFF_SUPPORT
    cout << endl << TIFFGetVersion() << endl;
  #endif
}

const std::string LISA_USAGE={
"lisa [-/--options]\n\n"
"-a,--analyze      analyze connected components of a raster file (#-area cells)\n"
"  -f,--fragment   statistics using minimum fragment size in [ha] (default: 0)\n"
"  -s,--save       save results to .csv file, 1=small, 2=large\n"
"  -p,--pixel      minimum pixel length for edge detection, default: 1\n"
"  -w,--write      write clusterlabel data (bin/lab), 1=clusters+labels, 2=labels\n"
"  -b              mean biomass\n"
"  --nrows         number of rows to process\n"
"  --agb-file      saatchi agb biomass file [t/ha]\n"
"  --bthres        biomass threshold [t/ha] (default: 0 t/ha)\n"
"  --check         check consistency of component analysis\n"
"  --rloss         relative carbon loss in edge areas, default: 0.5\n"
"  --flush         flush clusters to use fixed amount of memory\n"
"-c,--convert      convert raster file into [bri] file (g=globcover)\n"
"-m,--map          produce a density map, out of bin/lab file\n"
"  --map-scale     number of output classes 1...256\n"
"  --map-type      0=closs, 1=core/area\n"
"  --map-reduction reduction factor, default: 500\n"
"  -e,--extend     top,left,right,bottom\n"
"-t,--test         test consistence of lisa\n"
"-v,--verbose      verbosity level [0-2] (default: 1)\n"
"-d,--dept         edge effect dept of d [m], default: 100\n"
"--info            info about raster file\n"
"--input           inputfile\n"
"--output          outputfile\n"
"--version         print version info\n"
"--threshold       forest cover threshold for forest/nonforest map\n"
"--force           force overwrite of files\n"
"supported raster file formats: asc, pgm, tiff, bri\n"
};

int main(int argc,char *argv[])
{
    AnalyzeOptions AnalyzeOptions;

    int verbosity_level=1;
    ComandLine::METHOD cmode=ComandLine::ANALYZE;

    bool globcover=false;
    int reduction_factor=500;
    int map_scale=256;
    int map_type=0;
    bool force_overwrite=false;
    std::string str_ifile,str_ofile;
    std::string str_bfile;
    geoExtend myGeoExtend;

    if (argc < 2) {
      cout << LISA_USAGE << endl;
      return 1;
    }

    CmdOptions myCmdOpt(argc,argv);

    if (myCmdOpt.searchOption("-e","--extend")) {
      std::vector <double>vextend;
      myCmdOpt.getOpt(vextend);
      if (vextend.size()!=4) cout << "warning: unexpected size of extend=" << vextend.size() << endl;
      else {
        myGeoExtend.top=vextend[0];
        myGeoExtend.left=vextend[1];
        myGeoExtend.right=vextend[2];
        myGeoExtend.bottom=vextend[3];
      }
    }
    if (myCmdOpt.searchOption("-a","--analyze")) cmode=ComandLine::ANALYZE;
    if (myCmdOpt.searchOption("","--info")) cmode=ComandLine::INFO;
    if (myCmdOpt.searchOption("-c","--convert")) cmode=ComandLine::CONVERT;
    if (myCmdOpt.searchOption("-t","--test")) cmode=ComandLine::TEST;
    if (myCmdOpt.searchOption("-m","--map")) cmode=ComandLine::MAP;
    if (myCmdOpt.searchOption("","--version")) cmode=ComandLine::VERSION;
    if (myCmdOpt.searchOption("","--force")) force_overwrite=true;
    if (myCmdOpt.searchOption("-v","--verbose")) myCmdOpt.getOpt(verbosity_level);
    if (myCmdOpt.searchOption("-d","--dept")) {myCmdOpt.getOpt(AnalyzeOptions.edge_dept);};
    if (myCmdOpt.searchOption("","--threshold")) {myCmdOpt.getOpt(AnalyzeOptions.forest_cover_threshold);};
    if (myCmdOpt.searchOption("","--nrows")) {myCmdOpt.getOpt(AnalyzeOptions.nrows);};
    if (myCmdOpt.searchOption("","--flush")) {AnalyzeOptions.flush_clusters=true;};
    if (myCmdOpt.searchOption("-f","--fragment")) {myCmdOpt.getOpt(AnalyzeOptions.min_fragment_size);};
    if (myCmdOpt.searchOption("-w","--write")) myCmdOpt.getOpt(AnalyzeOptions.write_mode);
    if (myCmdOpt.searchOption("-s","--save")) myCmdOpt.getOpt(AnalyzeOptions.save_mode);
    if (myCmdOpt.searchOption("-b","")) myCmdOpt.getOpt(AnalyzeOptions.mean_biomass);
    if (myCmdOpt.searchOption("-p","--pixel")) myCmdOpt.getOpt(AnalyzeOptions.pixel_len);
    if (myCmdOpt.searchOption("-i","--input")) myCmdOpt.getOpt(str_ifile);
    if (myCmdOpt.searchOption("-o","--output")) myCmdOpt.getOpt(str_ofile);
    if (myCmdOpt.searchOption("","--agb-file")) myCmdOpt.getOpt(str_bfile);
    if (myCmdOpt.searchOption("","--map-scale")) {myCmdOpt.getOpt(map_scale);map_scale=std::max(std::min(map_scale,256),1);};
    if (myCmdOpt.searchOption("","--map-type")) myCmdOpt.getOpt(map_type);
    if (myCmdOpt.searchOption("","--map-reduction")) myCmdOpt.getOpt(reduction_factor);
    if (myCmdOpt.searchOption("","--bthres")) myCmdOpt.getOpt(AnalyzeOptions.bthres);
    if (myCmdOpt.searchOption("","--rloss")) myCmdOpt.getOpt(AnalyzeOptions.relative_carbon_loss);
    if (myCmdOpt.searchOption("","--check")) AnalyzeOptions.check_consistency=true;

    if (verbosity_level>1) {
       cout << "mode:          " << cmode << endl;
       cout << "edge dept:     " << AnalyzeOptions.edge_dept << endl;
       cout << "mean biosmass: " << AnalyzeOptions.mean_biomass << endl;
       cout << "min fragment:  " << AnalyzeOptions.min_fragment_size << endl;
       cout << "pixel len:     " << AnalyzeOptions.pixel_len << endl;
       cout << "rel. c-loss:   " << AnalyzeOptions.relative_carbon_loss << endl;
       cout << "savemode:      " << AnalyzeOptions.save_mode << endl;
       cout << "writemode:     " << AnalyzeOptions.write_mode << endl;
       cout << "bthres:        " << AnalyzeOptions.bthres << endl;
       cout << "extend:        " << myGeoExtend.top<<","<<myGeoExtend.left<<","<<myGeoExtend.right<<","<<myGeoExtend.bottom<<endl;
       cout << "infile:        '" << str_ifile << "'" << endl;
       cout << "outfile:       '" << str_ofile << "'" << endl;
       cout << "agb-file:      '" << str_bfile << "'" << endl;
    }

    ComandLine myCmdLine;
    switch (cmode) {
      case ComandLine::ANALYZE:myCmdLine.Analyze(str_ifile,str_bfile,AnalyzeOptions);break;
      case ComandLine::CONVERT:
      case ComandLine::INFO:myCmdLine.Convert(str_ifile,str_ofile,cmode,globcover,force_overwrite,myGeoExtend);break;
      case ComandLine::TEST:myCmdLine.TestConsistency();break;
      case ComandLine::MAP:myCmdLine.CreateMap(str_ifile,str_ofile,reduction_factor,map_scale,map_type,AnalyzeOptions.edge_dept,myGeoExtend);break;
      case ComandLine::VERSION:PrintVersion(1);break;
      default: cout << "unknown mode: " << cmode << endl;break;
    }
    return 0;
}
