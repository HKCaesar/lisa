#include "cluster.h"

ClusterBRI::ClusterBRI(BRIOptions &Options)
:opt(Options),bri_width(opt.myIMG.GetWidth())
{
  //rowbuffer=new uint8_t[bri_width];
  lookahead_rows=Options.pixel_len;
  bufrows=(2*lookahead_rows)+1;
  wrows=new int64_t*[bufrows];
  for (int i=0;i<bufrows;i++) wrows[i]=new int64_t[bri_width];
  minlabel=0;
  num_1pixel=0;
}

ClusterBRI::~ClusterBRI()
{
  //if (rowbuffer) delete []rowbuffer,rowbuffer=0;
  for (int i=0;i<bufrows;i++) {if (wrows[i]) delete []wrows[i],wrows[i]=0;};
  delete []wrows,wrows=0;
}

int64_t ClusterBRI::FindRoot(int64_t x)
{
  while (cdata[x]>0) x=cdata[x];
  return x;
}

int64_t ClusterBRI::FindCollapse(int64_t x)
{
  int64_t y=FindRoot(x);

  while (cdata[x]>0) { // collapse the tree
    int z=cdata[x];
    cdata[x]=y;
    x=z;
  }
  return y;
}

int64_t ClusterBRI::GetNumRoots()
{
  int64_t num_roots=0;
  for (int64_t i=1;i<=max_cluster_label;i++)
    if (cdata[i]<0) num_roots++;
  return num_roots;
}

// calculate edge-effected area
double ClusterBRI::CalculateEdgeAreaDE(double area,double edge_len)
{
  double edge_area=0.;
  if (edge_len*opt.edge_dept>4.*area) edge_area=area;
  else  // calculate size-index from didham & ewers core area model
  {
    const double si=edge_len/(2.*sqrt(M_PI*area));
    edge_area=opt.edge_dept*(edge_len - (si*si)*M_PI*opt.edge_dept);
  }
  return edge_area;
}

double ClusterBRI::CalculateEdgeAreaCircle(double area)
{
  double r=sqrt(area/M_PI);
  if (r < opt.edge_dept) return area;
  else {
    double edge_area=area-((r-opt.edge_dept)*(r-opt.edge_dept)*M_PI);
    return edge_area;
  }
}


// Calculate C-Loss for a forest fragment in Gt
// input: biomass in t, area in m^2, edge_area in m^2
double ClusterBRI::CalculateCLoss(double biomass,double area_m2,double edge_area)
{
  double biomass_ha=(biomass*10000.)/area_m2; // [t/ha]
  double carbon_ha=0.5*biomass_ha;
  double c_loss=(opt.relative_carbon_loss*edge_area/10000.*carbon_ha)/1000000000.; //[Gt]
  return c_loss;
}

void ClusterBRI::AddClusterStats(int64_t parea,double carea,double cborder,double cbiomass)
{
   double areaha=carea/10000.;
   if (areaha>=opt.min_fragment_size) {
      if (areaha < 10.) myStats.num_clusters10ha++;
      if (areaha < 50.) myStats.num_clusters50ha++;
      myStats.cell_area+=parea;
      myStats.total_area+=carea;
      myStats.total_border_len+=cborder;
      double edge_area_de=CalculateEdgeAreaDE(carea,cborder);
      myStats.total_edge_area_de+=edge_area_de;
      myStats.total_edge_area_circle+=CalculateEdgeAreaCircle(carea);

      myStats.total_biomass+=(cbiomass/1000000000.); //Gt

      myStats.total_closs+=CalculateCLoss(cbiomass,carea,edge_area_de);
      if (carea>myStats.max_area) myStats.max_area=carea;
      myStats.num_clusters++;
    }
}

void ClusterBRI::CalculateStats()
{
  myStats.Reset(); // in case CalculateStats is called a second time
  if (opt.flush_clusters) { // read clusters from file
    ofs_clusterfile.open(opt.str_clusterflushfile,ios::binary|ios::in);
    if (!ofs_clusterfile.is_open()) {
      cout << "warning: could not open: '" << opt.str_clusterflushfile << "'\n";
      return;
    }
    uint8_t buffer[8*4];
    for (int64_t i=0;i<total_roots_written;i++)
    {
       ofs_clusterfile.read((char*)buffer,8*4);
       int64_t parea=Utils::Get64LH(buffer);
       double  carea=Utils::GetDouble(buffer+8);
       double  cborder=Utils::GetDouble(buffer+16);
       double  cbiomass=Utils::GetDouble(buffer+24);
       AddClusterStats(parea,carea,cborder,cbiomass);
    }
    ofs_clusterfile.close();
  } else {
    for (int64_t i=1;i<=max_cluster_label;i++) {
      if (cdata[i]<0) AddClusterStats(-cdata[i],clusterdata[i].area,clusterdata[i].border,clusterdata[i].biomass);
    }
  }
  if (myStats.num_clusters)  myStats.mean_area=myStats.total_area/((double)myStats.num_clusters*10000.);
}

double ClusterBRI::CalculateBorder(inter_cell &icell,bool left,bool right,bool top,bool bottom,int &border_pixel)
{
  border_pixel=0;
  double border_len=0.;
  if (left) {border_len+=icell.pixel_height;border_pixel++;};
  if (right) {border_len+=icell.pixel_height;border_pixel++;};
  if (top) {border_len+=icell.pixel_width;border_pixel++;};
  if (bottom) {border_len+=icell.pixel_width;border_pixel++;};
  return border_len;
}

void ClusterBRI::DetectBorders(int row,int cur_row,int i,bool &bleft,bool &bright,bool &btop,bool &bbottom)
{
  int j;
  bleft=btop=bright=bbottom=true;
  if (i>=opt.pixel_len) {
    for (j=1;j<=opt.pixel_len;j++) {if (wrows[cur_row][i-j]!=0) {bleft=false;break;}};
  }
  if (i<bri_width-opt.pixel_len) {
    for (j=1;j<=opt.pixel_len;j++) {if (wrows[cur_row][i+j]!=0) {bright=false;break;}};
  }
  if (row>=opt.pixel_len) {
    for (j=1;j<=opt.pixel_len;j++) {if (wrows[Utils::SMod(cur_row-j,bufrows)][i]!=0) {btop=false;break;}};
  }
  if (row<endrow-opt.pixel_len) {
    for (j=1;j<=opt.pixel_len;j++) {if (wrows[Utils::SMod(cur_row+j,bufrows)][i]!=0) {bbottom=false;break;}};
  }
}

// only for testing, should be obsolet
void ClusterBRI::WriteMarkedRow(int64_t *clusterrow,uint32_t width,FILE *file)
{
  int outsize=RLEPack::PackRow(clusterrow,width,clusterrowdata);
  RLEPack::UnpackRow(clusterrowdata,width,rowtmp);
  for (uint32_t i=0;i<width;i++) {
    if (rowtmp[i]!=clusterrow[i]) cout << "error in decompression\n" << endl;
  }
  uint8_t tbuf[4];Utils::Put32LH(tbuf,outsize);fwrite(tbuf,1,4,file);
  fwrite(clusterrowdata,1,outsize,file);
}

void ClusterBRI::ProcessRow(int row,int cur_row)
{
  inter_cell icell=opt.Proj.GetCellDim(row);
  int64_t *currow=wrows[cur_row];
  for (int i=0;i<bri_width;i++)
  {
    if (currow[i])
    {
      num_1pixel++;

      int64_t pleft,ptop; // load the left & top pixels for the hoshen-kopelman algorithm
      pleft=ptop=0;
      if (i>0) pleft=currow[i-1];
      if (row>0) ptop=wrows[Utils::SMod((cur_row-1),bufrows)][i];

      // detect borders
      bool bleft,bright,btop,bbottom;
      DetectBorders(row,cur_row,i,bleft,bright,btop,bbottom);

      int border_pixel;
      double border_len=CalculateBorder(icell,bleft,bright,btop,bbottom,border_pixel);
      max_border_pixel+=4;

      double biomass_m2=opt.BMass.getBiomassRef(i,row);

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
    }
  }
  if (opt.write_clusterlabel==1)WriteMarkedRow(currow,bri_width,clusterfile1);
}

void ClusterBRI::PrintHist(std::vector <int64_t> &hist,std::string header)
{
  cout << header << endl;
  for (size_t i=0;i<hist.size();i++) {
     if (i==0) cout << "[0-1):" << hist[i] << endl;
     else cout << "[10^" << (i-1) << "-10^"<<(i)<<"):" << hist[i] << endl;
  }
}

void ClusterBRI::PrintHist(std::vector <double> &hist,std::string header,std::string unit)
{
  cout << header << endl;
  for (size_t i=0;i<hist.size();i++) {
     if (i==0) cout << "[0-1):" << std::fixed << std::setprecision(4) << hist[i] << unit << endl;
     else cout << "[10^" << (i-1) << "-10^"<<(i)<<"):" << std::fixed << std::setprecision(4) << hist[i] << unit << endl;
  }
}

void ClusterBRI::WriteHist(ofstream &file,std::vector <int64_t> &hist,std::string header)
{
  file << header << "\n";
  for (size_t i=0;i<hist.size();i++)  {
    if (i==0) file << hist[i];
    else file << "," << hist[i];
  }
  file<<"\n";
}


void ClusterBRI::WriteHist(ofstream &file,std::vector <double> &hist,std::string header)
{
  file << header << "\n";
  for (size_t i=0;i<hist.size();i++)  {
    if (i==0) file << std::fixed << std::setprecision(4) << hist[i];
    else file << "," << std::fixed << std::setprecision(4) << hist[i];
  }
  file<<"\n";
}

void ClusterBRI::SaveSmallClusterData(std::string &fname)
{
  // calculate histograms
  std::vector <int64_t>hist_area(10);
  std::vector <double>hist_totalarea(10);
  std::vector <double>hist_totaledge(10);
  std::vector <double>hist_biomass(10);
  std::vector <double>hist_totalloss(10);
  for (int64_t i=1;i<=max_cluster_label;i++)
    if (cdata[i]<0 && clusterdata[i].area/10000.>opt.min_fragment_size) {
        const double area_m2=clusterdata[i].area;
        int dclass=floor(log10(area_m2/10000.));
        if (dclass<0) dclass=0;
        else if (dclass<9) dclass++;
        if (dclass>=0 && dclass<=9)
        {
          hist_area[dclass]++;
          hist_totalarea[dclass]+=Utils::SqMetre_To_MillHa(area_m2);

          double edge_area_de=CalculateEdgeAreaDE(area_m2,clusterdata[i].border);

          hist_totaledge[dclass]+=Utils::SqMetre_To_MillHa(edge_area_de);

          double biomass=(clusterdata[i].biomass);
          hist_biomass[dclass]+=biomass/1000000000.;
          hist_totalloss[dclass]+=CalculateCLoss(clusterdata[i].biomass,clusterdata[i].area,edge_area_de);
        } else cout << "warning: too large fragment detected: " << area_m2/10000. << " ha" << endl;
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
void ClusterBRI::SaveFullClusterData(std::string &fname)
{
  char stmp[256];
  /*std::stringstream sout;
  sout << std::fixed << std::setprecision(2);*/

  ofstream myfile (fname);
  if (myfile.is_open())
  {
    for (size_t i=0;i<cdata.size();i++)
    {
       if (cdata[i]<0 && clusterdata[i].area/10000.>opt.min_fragment_size)
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

void ClusterBRI::PrintProgress(int y,int height)
{
  cout << y << "/" << height << "\r";
}

// calculate carbon-loss for a given cluster label [t/ha]
double ClusterBRI::CalculateCLossPerHA(int64_t label)
{
  double area_m2=clusterdata[label].area;
  double edge_area_de=CalculateEdgeAreaDE(area_m2,clusterdata[label].border);
  if (edge_area_de < 0) cout << "warning: edge_area < 0\n";
  double c_loss=(0.5*opt.relative_carbon_loss*edge_area_de*clusterdata[label].biomass*10000.)/(area_m2*area_m2);
  return c_loss;
}

void ClusterBRI::WriteLabelFile()
{
  ofstream labelfile(opt.str_labelfile);
  for (int64_t i=1;i<=max_cluster_label;i++) {
     if (cdata[i]<0) {
        double closs=CalculateCLossPerHA(i);
        labelfile << i << "," << Utils::ConvertFixed(closs,4) << "\n";
     }
  }
  labelfile.close();
}

void ClusterBRI::WriteClusterfile()
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
  Utils::Put32LH(tbuf+4,endrow);
  fwrite(tbuf,1,8,clusterfile2);

  double total_biomass=0.0;
  int64_t total_cells=0;
  for (row=0;row<endrow;row++)
  {
    uint32_t nread=fread(tbuf,1,4,clusterfile1);
    if (nread!=4) cout << " error reading clusterfile1 at line: " << row << endl;
    else {
       uint32_t outsize=Utils::Get32LH(tbuf);
       nread=fread(clusterrowdata,1,outsize,clusterfile1);
       if (nread!=outsize) cout << " error reading clusterfile1 at line: " << row << endl;
       else {
         RLEPack::UnpackRow(clusterrowdata,bri_width,labelrow);
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
    if (opt.verbose && (row+1)%100==0) PrintProgress(row+1,endrow);
  }
  cout << "total cells:   " << total_cells << endl;
  cout << "total biomass: " << total_biomass/1000000000.0 << endl;

  fclose(clusterfile1);
  fclose(clusterfile2);
  cout << "deleting '" << opt.str_clusterfile1 << "'\n";
  std::remove(opt.str_clusterfile1.c_str());
}

// check connected components for consistency
void ClusterBRI::CheckClusters()
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

int ClusterBRI::UnpackRow(int64_t *dstrow,uint8_t *srcrow,int len)
{
  for (int i=0;i<len;i++) {
    uint8_t val=srcrow[i];
    if (opt.forest_cover_threshold>0) {
      if (val>opt.forest_cover_threshold) val=1;
      else val=0;
    }
    if (val!=0 && val!=1) {Utils::PrintWarning("invalid input value: "+std::to_string(val));return 1;}
    dstrow[i]=val;
  }
  return 0;
}

// attempts to compress the cdata tree
// unfortunately uses temporary twice as much memory, needs to be fixed
void ClusterBRI::CompressTree(int cur_row)
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
void ClusterBRI::FlushClusters(int cur_row)
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
}

void ClusterBRI::ClusterAnalyzation()
{
  if (opt.write_clusterlabel==1) {
    cout << "writing pass1-clusterfile: '" << opt.str_clusterfile1 << "'" << endl;
    clusterfile1=Utils::OpenWriteFILE(opt.str_clusterfile1);
    if (clusterfile1==NULL) {
      cout << "warning: could not create '" << opt.str_clusterfile1 << "'\n";
      return;
    }
    labelrow=new int64_t[bri_width];
    rowtmp=new int64_t[bri_width];
    clusterrowdata=new uint8_t[bri_width*8];
  }
  if (opt.flush_clusters) {
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
  endrow=opt.myIMG.GetHeight();
  //endrow=20000;

  max_cluster_label=0;
  total_roots_written=0;
  cdata.push_back(max_cluster_label); // dummy label
  tcelldata tcell;tcell.area=tcell.border=0.;tcell.biomass=0.;
  clusterdata.push_back(tcell);

  int lookahead=0;
  for (lookahead=0;lookahead<lookahead_rows;lookahead++) // buffer lookahead_rows
  {
    opt.myIMG.ReadRow();
    UnpackRow(wrows[lookahead],opt.myIMG.rowbuffer,bri_width);
  }
  int row_ptr=lookahead;
  int cur_ptr=0;

  for (row=0;row<endrow;row++)
  {
      if (row < endrow-lookahead_rows) { // read in the next row, process cur_row
        opt.myIMG.ReadRow();
        UnpackRow(wrows[row_ptr],opt.myIMG.rowbuffer,bri_width);
      }
      ProcessRow(row,cur_ptr);

      //if (row%10000==0) CompressTree(cur_ptr);
      if (opt.flush_clusters && (row+1)%100==0) FlushClusters(cur_ptr);

      cur_ptr=(cur_ptr+1)%bufrows;
      row_ptr=(row_ptr+1)%bufrows;

      if (opt.verbose && (row+1)%100==0) PrintProgress(row+1,endrow);
  }
  if (opt.flush_clusters) {
     FlushClusters(-1);
     ofs_clusterfile.close();
     cout << "flush_clusters: " << total_roots_written << " roots written\n";
  }

  if (opt.verbose) {PrintProgress(endrow,endrow);cout << endl;};
  opt.myIMG.StopReader();

  if (opt.write_clusterlabel==1) {
    fclose(clusterfile1);
    WriteClusterfile();
    WriteLabelFile();

    delete []labelrow;
    delete []clusterrowdata;
    delete []rowtmp;
  } else if (opt.write_clusterlabel==2) {
    WriteLabelFile();
  }

  CalculateStats();
  if (std::remove(opt.str_clusterflushfile.c_str())!=0)  {
     cout << "warning: could not delete '"  << opt.str_clusterflushfile << "'\n";
  }

  if (opt.verbose) {
  cout << "cdata size: " << ((cdata.size()*sizeof(int64_t))>>20) << " MiB";
  cout << ", clusterdata size: " << ((clusterdata.size()*sizeof(tcelldata))>>20) << " MiB" << endl;
  cout << endl;

  cout << "number of clusters:  " << myStats.num_clusters << " (min size: " << opt.min_fragment_size << " ha, forest cover: " << opt.forest_cover_threshold << "%)" << endl;
  if (myStats.num_clusters) {
  double r1=(double)myStats.num_clusters10ha/(double)myStats.num_clusters*100.0;
  double r2=(double)myStats.num_clusters50ha/(double)myStats.num_clusters*100.0;
  cout << "fraction < 10ha   :  " << std::fixed << std::setprecision(2) << r1 << " %" << endl;
  cout << "fraction < 50ha   :  " << std::fixed << std::setprecision(2) << r2 << " %" << endl;
  }
  cout << "total area:          " << myStats.cell_area << " pixels = " << std::fixed << std::setprecision(2) << Utils::SqMetre_To_MillHa(myStats.total_area) << " 10^6 ha" << endl;
  cout << "mean area:           " << std::fixed << std::setprecision(4) << myStats.mean_area << " ha" << endl;
  cout << "max area:            " << std::fixed << std::setprecision(4) << Utils::SqMetre_To_MillHa(myStats.max_area) << " 10^6 ha" << endl;
  cout << "edge len:            " << std::fixed << std::setprecision(4) << Utils::Metre_To_MillKm(myStats.total_border_len) << " 10^6 km (pixel-len: " << opt.pixel_len << ")" << endl;
  //cout << "max edge len:        " << std::fixed << std::setprecision(4) << Utils::Metre_To_MillKm((double)max_border_pixel*sqrt(opt.Proj.GetMeanPixelArea())) << " 10^6 km" << endl;
  cout << "edge area (DE):      " << std::fixed << std::setprecision(4) << Utils::SqMetre_To_MillHa(myStats.total_edge_area_de) << " 10^6 ha" << ", edge effect dept: " << std::fixed << std::setprecision(1) << opt.edge_dept << " m" << endl;
  if (myStats.total_area)
  cout << "edge area/area:      " << std::fixed << std::setprecision(2) << (myStats.total_edge_area_de/myStats.total_area*100) << " %" << endl;
  cout << "edge area (Circle):  " << std::fixed << std::setprecision(4) << Utils::SqMetre_To_MillHa(myStats.total_edge_area_circle) << " 10^6 ha" << endl;
  if (myStats.total_area)
  cout << "edge area/area:      " << std::fixed << std::setprecision(2) << (myStats.total_edge_area_circle/myStats.total_area*100) << " %" << endl;
  cout << endl;
  cout << "total biomass:       " << std::fixed << std::setprecision(2) << myStats.total_biomass << " Gt" << endl;
  cout << "mean biomass:        " << std::fixed << std::setprecision(2) << (myStats.total_biomass*1000)/(Utils::SqMetre_To_MillHa(myStats.total_area)) << " t/ha" << endl;
  cout << "total C-stock:       " << std::fixed << std::setprecision(2) << (0.5*myStats.total_biomass) << " Gt" << endl;
  cout << "total C-loss:        " << std::fixed << std::setprecision(2) << myStats.total_closs<< " Gt (rel. loss e=" <<opt.relative_carbon_loss<<")"<< endl;
  }
}
