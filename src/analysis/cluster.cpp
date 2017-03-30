#include "cluster.h"
#include "..\model\sic.h"

Cluster::Cluster(BRIOptions &Options)
:opt(Options),bri_width(opt.myIMG.GetWidth()),mask_rows(2,std::vector<char>(bri_width)),bufrows((2*(opt.analyze_opt.pixel_len+1))+1),
vborder(4),vcorridor(4)
{
  //rowbuffer=new uint8_t[bri_width];
  lookahead_rows=Options.analyze_opt.pixel_len+1;
  wrows=new int64_t*[bufrows];
  for (int i=0;i<bufrows;i++) wrows[i]=new int64_t[bri_width];

  minlabel=0;
  num_1pixel=0;
}

Cluster::~Cluster()
{
  //if (rowbuffer) delete []rowbuffer,rowbuffer=0;
  for (int i=0;i<bufrows;i++) {if (wrows[i]) delete []wrows[i],wrows[i]=0;};
  delete []wrows,wrows=0;
}

int64_t Cluster::FindRoot(int64_t x)
{
  while (cdata[x]>0) x=cdata[x];
  return x;
}

int64_t Cluster::FindCollapse(int64_t x)
{
  int64_t y=FindRoot(x);

  while (cdata[x]>0) { // collapse the tree
    int z=cdata[x];
    cdata[x]=y;
    x=z;
  }
  return y;
}

int64_t Cluster::GetNumRoots()
{
  int64_t num_roots=0;
  for (int64_t i=1;i<=max_cluster_label;i++)
    if (cdata[i]<0) num_roots++;
  return num_roots;
}

// calculate edge-effected area
double Cluster::CalculateEdgeAreaDE(double area,double edge_len,double edge_effect_dept)
{
  double edge_area=0.;
  if (edge_len*edge_effect_dept>4.*area) edge_area=area;
  else  // calculate size-index from didham & ewers core area model
  {
    const double si=edge_len/(2.*sqrt(M_PI*area));
    edge_area=edge_effect_dept*(edge_len - (si*si)*M_PI*edge_effect_dept);
  }
  return edge_area;
}

// calculate edge-effected area
double Cluster::CalculateEdgeAreaDE(double area,double edge_len)
{
  return CalculateEdgeAreaDE(area,edge_len,opt.analyze_opt.edge_dept);
}

double Cluster::CalculateEdgeAreaCircle(double area)
{
  double r=sqrt(area/M_PI);
  if (r < opt.analyze_opt.edge_dept) return area;
  else {
    double edge_area=area-((r-opt.analyze_opt.edge_dept)*(r-opt.analyze_opt.edge_dept)*M_PI);
    return edge_area;
  }
}


// Calculate C-Loss for a forest fragment in Gt
// input: biomass in t, area in m^2, edge_area in m^2
double Cluster::CalculateCLoss(double biomass,double area_m2,double edge_area)
{
  double biomass_ha=(biomass*10000.)/area_m2; // [t/ha]
  double carbon_ha=0.5*biomass_ha;
  double c_loss=(opt.analyze_opt.relative_carbon_loss*edge_area/10000.*carbon_ha)/1000000000.; //[Gt]
  return c_loss;
}

void Cluster::AddClusterStats(int64_t parea,const tcelldata &cell)
{
   double areaha=cell.area/10000.;
   if (areaha>=opt.analyze_opt.min_fragment_size) {
      if (areaha < 10.) myStats.num_clusters10ha++;
      if (areaha < 50.) myStats.num_clusters50ha++;
      myStats.cell_area+=parea;
      myStats.total_area+=cell.area;
      myStats.total_border_len+=cell.border;
      double edge_area_de=CalculateEdgeAreaDE(cell.area,cell.border);
      myStats.total_edge_area_de+=edge_area_de;
      myStats.total_edge_area_circle+=CalculateEdgeAreaCircle(cell.area);

      double rel_area=((cell.area-edge_area_de)*100.)/cell.area; // each fragment has one of four state conditions
      if (rel_area>75) {myStats.fragment_state_total[3]++;myStats.fragment_state_area[3]+=cell.area;}
      else if (rel_area>50) {myStats.fragment_state_total[2]++;myStats.fragment_state_area[2]+=cell.area;}
      else if (rel_area>25) {myStats.fragment_state_total[1]++;myStats.fragment_state_area[1]+=cell.area;}
      else  {myStats.fragment_state_total[0]++;myStats.fragment_state_area[0]+=cell.area;}

      myStats.total_biomass+=(cell.biomass/1000000000.); //Gt

      myStats.total_closs+=CalculateCLoss(cell.biomass,cell.area,edge_area_de);

      if (cell.area>myStats.max_area) myStats.max_area=cell.area;
      if (cell.area<myStats.min_area) myStats.min_area=cell.area;
      myStats.num_clusters++;
    }
}

void Cluster::CalculateStats()
{
  if (opt.analyze_opt.flush_clusters) { // read clusters from file
    ofs_clusterfile.open(opt.str_clusterflushfile,ios::binary|ios::in);
    if (!ofs_clusterfile.is_open()) {
      cout << "warning: could not open: '" << opt.str_clusterflushfile << "'\n";
      return;
    }
    uint8_t buffer[8*4];
    for (int64_t i=0;i<total_roots_written;i++)
    {
       ofs_clusterfile.read((char*)buffer,8*4);
       tcelldata cell;
       int64_t parea=Utils::Get64LH(buffer);
       cell.area=Utils::GetDouble(buffer+8);
       cell.border=Utils::GetDouble(buffer+16);
       cell.biomass=Utils::GetDouble(buffer+24);
       AddClusterStats(parea,cell);
    }
    ofs_clusterfile.close();
  } else {
    for (int64_t i=1;i<=max_cluster_label;i++) {
      if (cdata[i]<0) AddClusterStats(-cdata[i],clusterdata[i]);
    }
  }
  if (myStats.num_clusters)  myStats.mean_area=myStats.total_area/((double)myStats.num_clusters*10000.);
}

double Cluster::CalculateBorder(inter_cell &icell,bool left,bool right,bool top,bool bottom)
{
  double border_len=0.;
  if (left) {border_len+=icell.pixel_height;};
  if (right) {border_len+=icell.pixel_height;};
  if (top) {border_len+=icell.pixel_width_top;};
  if (bottom) {border_len+=icell.pixel_width_bottom;};
  return border_len;
}

void Cluster::CalculateBorder(inter_cell &icell,double &border_len,double &corridor_len)
{
  border_len=corridor_len=0.;
  if (vborder[0]) {border_len+=icell.pixel_height;if (vcorridor[0]) corridor_len+=icell.pixel_height;};
  if (vborder[1]) {border_len+=icell.pixel_height;if (vcorridor[1]) corridor_len+=icell.pixel_height;};
  if (vborder[2]) {border_len+=icell.pixel_width_top;if (vcorridor[2]) corridor_len+=icell.pixel_width_top;};
  if (vborder[3]) {border_len+=icell.pixel_width_bottom;if (vcorridor[3]) corridor_len+=icell.pixel_width_bottom;};
}

void Cluster::DetectBorders(int row,int cur_row,int i)
{
  int j;
  const int pixel_len=opt.analyze_opt.pixel_len;
  const int pixel_len1=pixel_len+1;

  std::fill(begin(vborder),end(vborder),true);
  std::fill(begin(vcorridor),end(vcorridor),true);
  if (i>=pixel_len) {
    for (j=1;j<=pixel_len;j++) {if (wrows[cur_row][i-j]!=0) {vborder[0]=false;break;}};
    if (vborder[0]) {if (i>=pixel_len1 && wrows[cur_row][i-pixel_len1]==0) vcorridor[0]=false;};
  }
  if (i<bri_width-pixel_len) {
    for (j=1;j<=pixel_len;j++) {if (wrows[cur_row][i+j]!=0) {vborder[1]=false;break;}};
    if (vborder[1]) {if (i<bri_width-pixel_len1 && wrows[cur_row][i+pixel_len1]==0) vcorridor[1]=false;};
  }
  if (row>=pixel_len) {
    int trow=cur_row;
    for (j=1;j<=pixel_len;j++) {
      if (--trow<0) trow=bufrows-1;
      if (wrows[trow][i]!=0) {vborder[2]=false;break;};
    };
    if (vborder[2]) {
      if (--trow<0) trow=bufrows-1;
      if (row>=pixel_len1 && wrows[trow][i]==0) vcorridor[2]=false;
    }
  }
  if (row<num_rows-pixel_len) {
    int trow=cur_row;
    for (j=1;j<=pixel_len;j++) {
      if (++trow>=bufrows) trow=0;
      if (wrows[trow][i]!=0) {vborder[3]=false;break;}
    };
    if (vborder[3]) {
      if (++trow>=bufrows) trow=0;
      if (row<num_rows-pixel_len1 && wrows[trow][i]==0) vcorridor[3]=false;
    }
    //for (j=1;j<=pixel_len;j++) {if (wrows[Utils::SMod(cur_row+j,bufrows)][i]!=0) {vborders[3]=false;break;}};
    //if (row<num_rows-pixel_len1 && wrows[Utils::SMod(cur_row+pixel_len1,bufrows)][i]==0) vcorridor[3]=false;
  }
}

void Cluster::DetectBorders(int row,int cur_row,int i,bool &bleft,bool &bright,bool &btop,bool &bbottom)
{
  int j;
  const int pixel_len=opt.analyze_opt.pixel_len;

  bleft=btop=bright=bbottom=true;
  if (i>=pixel_len) {
    for (j=1;j<=pixel_len;j++) {if (wrows[cur_row][i-j]!=0) {bleft=false;break;}};
  }
  if (i<bri_width-pixel_len) {
    for (j=1;j<=pixel_len;j++) {if (wrows[cur_row][i+j]!=0) {bright=false;break;}};
  }
  if (row>=pixel_len) {
    for (j=1;j<=pixel_len;j++) {if (wrows[Utils::SMod(cur_row-j,bufrows)][i]!=0) {btop=false;break;}};
  }
  if (row<num_rows-pixel_len) {
    for (j=1;j<=pixel_len;j++) {if (wrows[Utils::SMod(cur_row+j,bufrows)][i]!=0) {bbottom=false;break;}};
  }
}

// only for testing, should be obsolet
void Cluster::WriteMarkedRow(int64_t *clusterrow,uint32_t width,FILE *file)
{
  memset(rowtmp,0,width*sizeof(int64_t));
  uint32_t outsize=RLEPack::PackRow(clusterrow,width,clusterrowdata_);
  if (outsize>maxrowdatasize_) cout << "cluster: outsize>maxrowdatasize\n";
  RLEPack::UnpackRow(clusterrowdata_,width,rowtmp);
  int nerr=0;
  for (uint32_t i=0;i<width;i++) {
      if (rowtmp[i]!=clusterrow[i]) {cout << "error at pos " << i << " " << rowtmp[i] << "!=" << clusterrow[i] << endl;nerr++;};
  }
  if (nerr) cout << "cluster: " << nerr << " errors in decompression ("<<outsize<<")\n";
  //else cout << "decompression ok\n";
  uint8_t tbuf[4];Utils::Put32LH(tbuf,outsize);fwrite(tbuf,1,4,file);
  fwrite(&clusterrowdata_[0],1,outsize,file);
}

// we spend most processing time here
void Cluster::ProcessRow(int irow,int row_offset,int cur_row,int mask_ptr,int corner_left,int corner_right)
{
  ShapeFile &myShapeFile=opt.SF;
  int numrow=irow+row_offset;

  inter_cell icell=opt.Proj.GetCellDim(numrow);

  int64_t *currow=wrows[cur_row];
  for (int i=corner_left;i<corner_right;i++)
  {
    if (currow[i] && (mask_rows[mask_ptr][i]=myShapeFile.IsCLASS(i,numrow))) {
      if (opt.analyze_opt.calc_surface_area) myStats.surface_area+=icell.pixel_area;
      num_1pixel++;

      int64_t pleft=0,ptop=0; // load the left & top pixels for the hoshen-kopelman algorithm
      if (i>0) pleft=mask_rows[mask_ptr][i-1]?currow[i-1]:0; // beware: mask_rows is invalid if currow[i-1]=0, but pleft=0, independent of the state of mask_rows
      if (irow>0) ptop=mask_rows[mask_ptr?0:1][i]?wrows[(cur_row-1)<0?(bufrows-1):(cur_row-1)][i]:0; // same applies for prev row

      double border_len=0,corridor_len=0;
      DetectBorders(irow,cur_row,i);
      CalculateBorder(icell,border_len,corridor_len);
      myStats.total_corridor_len+=corridor_len;

      max_border_pixel+=4;

      double biomass_m2=opt.BMass.getBiomassRef(i,numrow);

      if (pleft==0 && ptop==0) // new label
      {
         max_cluster_label++;
         currow[i]=max_cluster_label;

         tcelldata tcell;tcell.area=icell.pixel_area;tcell.border=border_len;tcell.biomass=icell.pixel_area*biomass_m2;
         cdata.push_back(-1);
         clusterdata.push_back(tcell);
      } else if (ptop==0 || pleft==0 || (ptop==pleft)) //cluster already seen, but no conflict
      {
        int64_t root=FindCollapse(std::max(pleft,ptop));
        currow[i]=root;
        cdata[root]--;
        clusterdata[root].area+=icell.pixel_area;
        clusterdata[root].border+=border_len;
        clusterdata[root].biomass+=(icell.pixel_area*biomass_m2);
      } else // conflict: map the larger root to the smaller one, collapsing the individual trees
      {
        int64_t root1=FindCollapse(pleft);
        int64_t root2=FindCollapse(ptop);
        int64_t lmin=std::min(root1,root2);
        int64_t lmax=std::max(root1,root2);
        if (lmin!=lmax)
        {
          cdata[lmin]+=cdata[lmax];
          cdata[lmax]=lmin;
          clusterdata[lmin].area+=clusterdata[lmax].area;
          clusterdata[lmin].border+=clusterdata[lmax].border;
          clusterdata[lmin].biomass+=clusterdata[lmax].biomass;
        }
        currow[i]=lmin;
        cdata[lmin]--;
        clusterdata[lmin].area+=icell.pixel_area;
        clusterdata[lmin].border+=border_len;
        clusterdata[lmin].biomass+=(icell.pixel_area*biomass_m2);
      }
    } else {
      if (opt.analyze_opt.calc_surface_area && myShapeFile.IsCLASS(i,numrow)) myStats.surface_area+=icell.pixel_area;
    }
  }
  if (opt.analyze_opt.write_mode==1)WriteMarkedRow(currow,bri_width,clusterfile1);
}

void Cluster::PrintHist(std::vector <int64_t> &hist,std::string header)
{
  cout << header << endl;
  for (size_t i=0;i<hist.size();i++) {
     if (i==0) cout << "[0-1):" << hist[i] << endl;
     else cout << "[10^" << (i-1) << "-10^"<<(i)<<"):" << hist[i] << endl;
  }
}

void Cluster::PrintHist(std::vector <double> &hist,std::string header,std::string unit)
{
  cout << header << endl;
  for (size_t i=0;i<hist.size();i++) {
     if (i==0) cout << "[0-1):" << std::fixed << std::setprecision(4) << hist[i] << unit << endl;
     else cout << "[10^" << (i-1) << "-10^"<<(i)<<"):" << std::fixed << std::setprecision(4) << hist[i] << unit << endl;
  }
}

void Cluster::WriteHist(ofstream &file,std::vector <int64_t> &hist,std::string header)
{
  file << header << "\n";
  for (size_t i=0;i<hist.size();i++)  {
    if (i==0) file << hist[i];
    else file << "," << hist[i];
  }
  file<<"\n";
}


void Cluster::WriteHist(ofstream &file,std::vector <double> &hist,std::string header)
{
  file << header << "\n";
  for (size_t i=0;i<hist.size();i++)  {
    if (i==0) file << std::fixed << std::setprecision(4) << hist[i];
    else file << "," << std::fixed << std::setprecision(4) << hist[i];
  }
  file<<"\n";
}

void Cluster::AddClusterSmallStats(const tcelldata &cell,vector <int64_t>&hist_area,vector <double>&hist_totalarea,vector <double>&hist_totaledge,vector <double>&hist_biomass,vector <double>&hist_totalloss)
{
    if (cell.area/10000.>opt.analyze_opt.min_fragment_size) {
        int dclass=floor(log10(cell.area/10000.));
        if (dclass<0) dclass=0;
        else if (dclass<9) dclass++;
        if (dclass>=0 && dclass<=9)
        {
          hist_area[dclass]++;
          hist_totalarea[dclass]+=Utils::SqMetre_To_MillHa(cell.area);

          double edge_area_de=CalculateEdgeAreaDE(cell.area,cell.border);

          hist_totaledge[dclass]+=Utils::SqMetre_To_MillHa(edge_area_de);

          hist_biomass[dclass]+=cell.biomass/1000000000.;
          hist_totalloss[dclass]+=CalculateCLoss(cell.biomass,cell.area,edge_area_de);

        } else cout << "warning: too large fragment detected: " << cell.area/10000. << " ha" << endl;
    }
}

void Cluster::SaveSmallClusterData(std::string &fname)
{
  // calculate histograms
  std::vector <int64_t>hist_area(10);
  std::vector <double>hist_totalarea(10);
  std::vector <double>hist_totaledge(10);
  std::vector <double>hist_biomass(10);
  std::vector <double>hist_totalloss(10);

  if (opt.analyze_opt.flush_clusters) { // read clusters from file
    ofs_clusterfile.open(opt.str_clusterflushfile,ios::binary|ios::in);
    if (!ofs_clusterfile.is_open()) {
      cout << "warning: could not open: '" << opt.str_clusterflushfile << "'\n";
      return;
    }
    uint8_t buffer[8*4];
    for (int64_t i=0;i<total_roots_written;i++)
    {
       ofs_clusterfile.read((char*)buffer,8*4);
       tcelldata cell;
       //int64_t parea=Utils::Get64LH(buffer);
       cell.area=Utils::GetDouble(buffer+8);
       cell.border=Utils::GetDouble(buffer+16);
       cell.biomass=Utils::GetDouble(buffer+24);
       AddClusterSmallStats(cell,hist_area,hist_totalarea,hist_totaledge,hist_biomass,hist_totalloss);
    }
    ofs_clusterfile.close();
  } else {
    for (int64_t i=1;i<=max_cluster_label;i++)
      if (cdata[i]<0) {
        AddClusterSmallStats(clusterdata[i],hist_area,hist_totalarea,hist_totaledge,hist_biomass,hist_totalloss);
      }
  }

  PrintHist(hist_area,"fragment distribution");
  PrintHist(hist_totalarea,"area distribution"," 10^6 ha");
  PrintHist(hist_totaledge,"edge area distribution"," 10^6 ha");
  PrintHist(hist_biomass,"biomass distribution"," Gt");
  PrintHist(hist_totalloss,"c-loss distribution"," Gt");

  ofstream myfile (fname);
  if (myfile.is_open()) {
    WriteHist(myfile,hist_area,"fragment distribution");
    WriteHist(myfile,hist_totalarea,"area distribution (10^6 ha)");
    WriteHist(myfile,hist_totaledge,"edge area distribution (10^6 ha)");
    WriteHist(myfile,hist_biomass,"biomass distribution (Gt)");
    WriteHist(myfile,hist_totalloss,"c-loss distribution (Gt)");
    myfile.close();
  }
}

// save full clusters
void Cluster::SaveFullClusterData(std::string &fname)
{
  char stmp[256];
  /*std::stringstream sout;
  sout << std::fixed << std::setprecision(2);*/

  ofstream myfile (fname);
  if (myfile.is_open())
  {
    for (size_t i=0;i<cdata.size();i++)
    {
       if (cdata[i]<0 && clusterdata[i].area/10000.>opt.analyze_opt.min_fragment_size)
       {
          snprintf(stmp,256,"%jd,%0.2f,%0.2f\n",(intmax_t)-cdata[i],clusterdata[i].area,clusterdata[i].border);
          myfile << stmp;

          /* too slow
          sout.clear();
          sout << (-cdata[i])<<","<<clusterdata[i].area<<","<<clusterdata[i].border<<"\n";*/

       }
    }
    myfile.close();
  }
  else cout << "Unable to open file";
}

void Cluster::PrintProgress(int y,int height)
{
  cout << y << "/" << height << "\r";
}

// calculate carbon-loss for a given cluster label [t/ha]
double Cluster::CalculateCLossPerHA(int64_t label)
{
  double area_m2=clusterdata[label].area;
  double edge_area_de=CalculateEdgeAreaDE(area_m2,clusterdata[label].border);
  if (edge_area_de < 0) cout << "warning: edge_area < 0\n";
  double c_loss=(0.5*opt.analyze_opt.relative_carbon_loss*edge_area_de*clusterdata[label].biomass*10000.)/(area_m2*area_m2);
  return c_loss;
}

//somehow overlaps with "savefullclusterdata"
void Cluster::WriteLabelFile()
{
  char stmp[256];
  ofstream labelfile(opt.str_labelfile);
  for (int64_t i=1;i<=max_cluster_label;i++) {
     if (cdata[i]<0) {
        snprintf(stmp,256,"%jd,%jd,%0.4f,%0.4f,%0.4f\n",(intmax_t)i,(intmax_t)-cdata[i],clusterdata[i].area,clusterdata[i].border,CalculateCLossPerHA(i));
        labelfile << stmp;
     }
  }
  labelfile.close();
}

void Cluster::WriteClusterfile()
{
  clusterfile1=NULL;
  FILE *clusterfile2=NULL;
  clusterfile1=fopen(opt.str_clusterfile1.c_str(),"rb");
  cout << "reading pass1-clusterfile: '" << opt.str_clusterfile1 << "'" << endl;

  cout << "writing pass2-clusterfile: '" << opt.str_clusterfile2 << "'" << endl;
  clusterfile2=Utils::OpenWriteFILE(opt.str_clusterfile2);
  if (clusterfile2==NULL)  {
    cout << "warning: could not create '" << opt.str_clusterfile2 << "'\n";
    return;
  }
  uint8_t tbuf[8];
  Utils::Put32LH(tbuf,bri_width);
  Utils::Put32LH(tbuf+4,num_rows);
  fwrite(tbuf,1,8,clusterfile2);

  double total_biomass=0.0;
  int64_t total_cells=0;
  for (int row=0;row<num_rows;row++)
  {
    uint32_t nread=fread(tbuf,1,4,clusterfile1);
    if (nread!=4) cout << " error reading clusterfile1 at line: " << row << endl;
    else {
       uint32_t compsize=Utils::Get32LH(tbuf);
       if (clusterrowdata_.size()<compsize) clusterrowdata_.resize(compsize);
       nread=fread(&clusterrowdata_[0],1,compsize,clusterfile1);

       if (nread!=compsize) cout << " error reading clusterfile1 at line: " << row << endl;
       else {
         RLEPack::UnpackRow(clusterrowdata_,bri_width,labelrow);
         for (int i=0;i<bri_width;i++) {
          if (labelrow[i]) {
            int64_t root=FindRoot(labelrow[i]);
            double biomass_per_cell=clusterdata[root].biomass/(-cdata[root]);
            total_cells++;
            total_biomass+=biomass_per_cell;
            labelrow[i]=root;
           }
         }
         WriteMarkedRow(labelrow,bri_width,clusterfile2);
       }
    }
    if (opt.analyze_opt.verbose && (row+1)%100==0) PrintProgress(row+1,num_rows);
  }
  cout << "total cells:   " << total_cells << endl;
  cout << "total biomass: " << total_biomass/1000000000.0 << endl;

  fclose(clusterfile1);
  fclose(clusterfile2);
  cout << "deleting '" << opt.str_clusterfile1 << "'\n";
  std::remove(opt.str_clusterfile1.c_str());
}

// check connected components for consistency
void Cluster::CheckClusters()
{
  cout << "checking clusters..." << endl;

  // test if the sum of all root-sizes matches "num_1pixel"
  int64_t num_1test=0;
  for (size_t i=0;i<cdata.size();i++) {
     if (cdata[i]<0) num_1test+=-cdata[i];
  }
  cout << "test sum of root-sizes: ";
  if (num_1test==num_1pixel) cout << "passed" << endl;
  else cout << "failed" << endl;
}

int Cluster::UnpackRow(int64_t *dstrow,uint8_t *srcrow,int len)
{
  for (int i=0;i<len;i++) {
    uint8_t val=srcrow[i];
    if (opt.analyze_opt.threshold_fc) { // thresholding activated?
      if (val>opt.analyze_opt.forest_cover_threshold) val=1;
      else val=0;
    }
    if (val!=0 && val!=1) {Utils::PrintWarning("invalid input value: "+std::to_string(val));return 1;}
    dstrow[i]=val;
  }
  return 0;
}

// attempts to compress the cdata tree
// unfortunately uses temporary twice as much memory, needs to be fixed
void Cluster::CompressTree(int cur_row)
{
   std::vector <int64_t> cdata_new;
   std::vector <int64_t> root_map;
   std::vector <tcelldata> clusterdata_new;

   root_map.resize(cdata.size());

   cdata_new.push_back(cdata[0]);
   clusterdata_new.push_back(clusterdata[0]);

   cout << "compressing tree: ";
   uint64_t max_cluster=0;
   for (size_t i=1;i<cdata.size();i++) {
     if (cdata[i]<0) {
        cdata_new.push_back(cdata[i]);
        clusterdata_new.push_back(clusterdata[i]);
        root_map[i]=(max_cluster+1);
        max_cluster++;
     }
   }
   max_cluster_label=max_cluster;

   int64_t *currow=wrows[cur_row];
   for (int i=0;i<bri_width;i++)
   {
     if (currow[i]) currow[i]=root_map[FindRoot(currow[i])];
   }
   double dratio=(max_cluster+1)*100./(double)cdata.size();
   cout << Utils::ConvertFixed(dratio,2) << "%" << endl;

   cdata=cdata_new;
   clusterdata=clusterdata_new;
}


// write all unused roots to cluster file
void Cluster::FlushClusters(int cur_row)
{
  std::vector <int64_t>root_used(cdata.size());
  std::vector <int64_t> cdata_new;
  std::vector <tcelldata> clusterdata_new;

  cdata_new.push_back(cdata[0]);
  clusterdata_new.push_back(clusterdata[0]);

  int64_t nroots=0;
  if (cur_row!=-1) {
    int64_t *currow=wrows[cur_row];
    for (int i=0;i<bri_width;i++) { // save the roots which are used in 'cur_row'
      if (currow[i]) {
        int64_t root=FindRoot(currow[i]);
        if (!root_used[root]) {
          root_used[root]=++nroots;
          cdata_new.push_back(cdata[root]);
          clusterdata_new.push_back(clusterdata[root]);
        }
        currow[i]=root_used[root];
      }
    }
    //cout << "numroots: " << nroots << ", " << (cdata.size()-nroots-1) << " unused labels" << endl;
  }

  int64_t roots_written=0;
  uint8_t buffer[8*4];
  for (size_t i=1;i<cdata.size();i++) {
    if (cdata[i]<0 && !root_used[i]) {
        roots_written++;
        Utils::Put64LH(buffer,-cdata[i]);
        Utils::PutDouble(buffer+8,clusterdata[i].area);
        Utils::PutDouble(buffer+16,clusterdata[i].border);
        Utils::PutDouble(buffer+24,clusterdata[i].biomass);
        ofs_clusterfile.write((char*)buffer,8*4);
    };
  }
  //cout << "roots written: " << roots_written << endl;
  total_roots_written+=roots_written;

  max_cluster_label=nroots;
  cdata=cdata_new;
  clusterdata=clusterdata_new;

  root_used.resize(0);
  cdata_new.resize(0);
  clusterdata_new.resize(0);
}

void Cluster::DeleteTempFiles()
{
 if (opt.analyze_opt.flush_clusters && std::remove(opt.str_clusterflushfile.c_str())!=0)  {
    cout << "warning: could not delete '"  << opt.str_clusterflushfile << "'\n";
 }
}

void Cluster::ClusterAnalyzation(const geoExtend &myextend)
{
  myStats.Reset(); // in case CalculateStats is called a second time

  int pleft,pright,ptop,pbottom;
  Frame::SetExtend(opt.Proj.getLeft(),opt.Proj.getTop(),opt.Proj.getCellsize(),myextend,opt.myIMG.GetWidth(),opt.myIMG.GetHeight(),pleft,ptop,pright,pbottom);
  std::cerr << "warning: SetExtend not fully implemented" << std::endl;
  num_rows=pbottom-ptop;

  //num_rows=1000;
  //std::cout <<
  //std::cout << myGeoExtend.top << " " << myGeoExtend.left << " " << myGeoExtend.right << " " << myGeoExtend.bottom << std::endl;

  if (opt.analyze_opt.write_mode==1) {
    cout << "writing pass1-clusterfile: '" << opt.str_clusterfile1 << "'" << endl;
    clusterfile1=Utils::OpenWriteFILE(opt.str_clusterfile1);
    if (clusterfile1==NULL) {
      cerr << "warning: could not create '" << opt.str_clusterfile1 << "'\n";
      return;
    }
    labelrow=new int64_t[bri_width];
    rowtmp=new int64_t[bri_width];
    maxrowdatasize_=bri_width*10;
  }
  if (opt.analyze_opt.flush_clusters) {
    ofs_clusterfile.open(opt.str_clusterflushfile,ios::binary|ios::out);
    if (!ofs_clusterfile.is_open()) {
       cout << "warning: could not create '" << opt.str_clusterflushfile << "'\n";
       return;
    }
  }
  max_border_pixel=0;
  cdata.resize(0);
  clusterdata.resize(0);
  opt.myIMG.StartReader();

  /*if (opt.analyze_opt.nrows>0) endrow=opt.analyze_opt.nrows;
  else endrow=opt.myIMG.GetHeight();*/

  max_cluster_label=0;
  total_roots_written=0;
  cdata.push_back(max_cluster_label); // dummy label
  tcelldata tcell;tcell.area=tcell.border=0.;tcell.biomass=0.;
  clusterdata.push_back(tcell);

  for (int i=0;i<ptop-1;i++) {
    opt.myIMG.ReadRow();
    if (i%10==0) cout << "skipping " << (ptop-1) << " lines: " << Utils::ConvertFixed(i*100/(double)(ptop-1),1) << "%\r";
  }
  std::cout << std::endl;
  //cout << "skipping " << (ptop-1) << " lines: " << Utils::ConvertFixed(100,1) << "%\n";

  int lookahead=0;
  for (lookahead=0;lookahead<lookahead_rows;lookahead++) // buffer lookahead_rows
  {
    opt.myIMG.ReadRow();
    UnpackRow(wrows[lookahead],opt.myIMG.rowbuffer,bri_width);
  }
  int row_ptr=lookahead;
  int cur_ptr=0;
  int mask_ptr=0;

  for (int row=0;row<num_rows;row++)
  {
      if (row < num_rows-lookahead_rows) { // read in the next row, process cur_row
        opt.myIMG.ReadRow();
        UnpackRow(wrows[row_ptr],opt.myIMG.rowbuffer,bri_width);
      }
      ProcessRow(row,ptop,cur_ptr,mask_ptr,pleft,pright);

      //if (row%10000==0) CompressTree(cur_ptr);
      if (opt.analyze_opt.flush_clusters && (row+1)%100==0) FlushClusters(cur_ptr);

      cur_ptr=(cur_ptr+1)%bufrows;
      row_ptr=(row_ptr+1)%bufrows;
      mask_ptr=(mask_ptr+1)%2;

      if (opt.analyze_opt.verbose && (row+1)%100==0) PrintProgress(row+1,num_rows);
  }
  if (opt.analyze_opt.flush_clusters) {
     FlushClusters(-1);
     ofs_clusterfile.close();
     cout << "flush_clusters: " << total_roots_written << " roots written\n";
  }

  if (opt.analyze_opt.verbose) {PrintProgress(num_rows,num_rows);std::cout << endl;};
  opt.myIMG.StopReader();

  if (opt.analyze_opt.write_mode==1) {
    fclose(clusterfile1);
    WriteClusterfile();
    WriteLabelFile();

    delete []labelrow;
    delete []rowtmp;
  } else if (opt.analyze_opt.write_mode==2) {
    WriteLabelFile();
  }

  CalculateStats();

  if (opt.analyze_opt.verbose) {
  cout << "cdata size: " << ((cdata.size()*sizeof(int64_t))>>20) << " MiB";
  cout << ", clusterdata size: " << ((clusterdata.size()*sizeof(tcelldata))>>20) << " MiB" << endl;
  cout << endl;

  cout << "number of clusters:  " << myStats.num_clusters << " (min size: " << opt.analyze_opt.min_fragment_size << " ha, forest cover: " << opt.analyze_opt.forest_cover_threshold << "%)" << endl;
  if (myStats.num_clusters) {
  double r1=(double)myStats.num_clusters10ha/(double)myStats.num_clusters*100.0;
  double r2=(double)myStats.num_clusters50ha/(double)myStats.num_clusters*100.0;
  cout << "fraction < 10ha   :  " << std::fixed << std::setprecision(2) << r1 << " %" << endl;
  cout << "fraction < 50ha   :  " << std::fixed << std::setprecision(2) << r2 << " %" << endl;
  }
  cout << "surface area:        " << std::fixed << std::setprecision(2) << Utils::SqMetre_To_MillHa(myStats.surface_area) << " 10^6 ha" << endl;
  cout << "total area:          " << myStats.cell_area << " pixels = " << std::fixed << std::setprecision(2) << Utils::SqMetre_To_MillHa(myStats.total_area) << " 10^6 ha" << endl;
  cout << "mean area:           " << std::fixed << std::setprecision(4) << myStats.mean_area << " ha" << endl;
  cout << "max area:            " << std::fixed << std::setprecision(4) << Utils::SqMetre_To_MillHa(myStats.max_area) << " 10^6 ha" << endl;
  cout << "min area:            " << std::fixed << std::setprecision(4) << (myStats.min_area/10000) << " ha" << endl;
  cout << "edge len:            " << std::fixed << std::setprecision(4) << Utils::Metre_To_MillKm(myStats.total_border_len) << " 10^6 km (pixel-len: " << opt.analyze_opt.pixel_len << ")" << endl;
  cout << "portion corridor:    " << std::fixed << std::setprecision(4) << Utils::Metre_To_MillKm(myStats.total_corridor_len) << " 10^6 km" << endl;
  //cout << "max edge len:        " << std::fixed << std::setprecision(4) << Utils::Metre_To_MillKm((double)max_border_pixel*sqrt(opt.Proj.GetMeanPixelArea())) << " 10^6 km" << endl;
  cout << "edge area (DE):      " << std::fixed << std::setprecision(4) << Utils::SqMetre_To_MillHa(myStats.total_edge_area_de) << " 10^6 ha" << ", edge effect dept: " << std::fixed << std::setprecision(1) << opt.analyze_opt.edge_dept << " m" << endl;
  if (myStats.total_area)
  cout << "edge area/area:      " << std::fixed << std::setprecision(2) << (myStats.total_edge_area_de/myStats.total_area*100) << " %" << endl;
  cout << "edge area (Circle):  " << std::fixed << std::setprecision(4) << Utils::SqMetre_To_MillHa(myStats.total_edge_area_circle) << " 10^6 ha" << endl;
  if (myStats.total_area)
  cout << "edge area/area:      " << std::fixed << std::setprecision(2) << (myStats.total_edge_area_circle/myStats.total_area*100) << " %" << endl;
  cout << endl;
  cout << "total biomass:       " << std::fixed << std::setprecision(2) << myStats.total_biomass << " Gt" << endl;
  cout << "mean biomass:        " << std::fixed << std::setprecision(2) << (myStats.total_biomass*1000)/(Utils::SqMetre_To_MillHa(myStats.total_area)) << " t/ha" << endl;
  cout << "total C-stock:       " << std::fixed << std::setprecision(2) << (0.5*myStats.total_biomass) << " Gt" << endl;
  cout << "total C-loss:        " << std::fixed << std::setprecision(2) << myStats.total_closs<< " Gt (rel. loss e=" <<opt.analyze_opt.relative_carbon_loss<<")"<< endl;
  std::cout << endl << "fragment state (core/total):" << std::endl;
  //if (myStats.num_clusters && (myStats.total_area>0.)) {
    for (int i=0;i<4;i++) {
      switch (i) {
        case 0:std::cout << "<=25%: ";break;
        case 1:std::cout << "<=50%: ";break;
        case 2:std::cout << "<=75%: ";break;
        case 3:std::cout << "else : ";break;
      }
      std::cout << std::setw(5) << std::setprecision(1) << ((double)myStats.fragment_state_total[i]*100.)/(double)myStats.num_clusters << "% clusters, ";
      std::cout << std::setw(5) << std::setprecision(1) << ((double)myStats.fragment_state_area[i]*100.)/(double)myStats.total_area << "% area" << std::endl;
    //}
  }
  }
}
