#include "proj.h"

int Projection::ReadCoordinateFile(const std::string &fname)
{
  ifstream myfile (fname);
  if (myfile.is_open())
  {
    std::string line;
    while ( getline (myfile,line) )
    {
      std::string rline=Utils::RemoveWhite(line);
      if (rline.length())
      {
        size_t pos=rline.find("=");
        if (pos!=string::npos) {
          std::string key=rline.substr(0,pos);
          std::string val=rline.substr(pos+1);
          Utils::ToUpper(key);
          if (key.compare("TOP")==0) top=std::stod(val);
          else if (key.compare("LEFT")==0) left=std::stod(val);
          else if (key.compare("RIGHT")==0) right=std::stod(val);
          else if (key.compare("BOTTOM")==0) bottom=std::stod(val);
        } else cout << "unknown option: '" << rline << "'" << endl;
      }
    }
    myfile.close();

    #if 1
    cout << endl;
      double height=CalcDeltaLat(bottom,top);
      double lat_mid=bottom+(top-bottom)/2;
      double width=CalcDeltaLong(lat_mid,right-left);
      cout << "proj: size at "<<Utils::ConvertFixed(lat_mid,2) << " degree: " << Utils::ConvertFixed(width/1000.0,2) << "km x " << (height/1000.0) << "km "<<endl;

    /*double lat_mid=top+(bottom-top)/2;
    double dwidth=CalcDist_Haversine(lat_mid,left,lat_mid,left+1);
    double dheight=CalcDist_Haversine(top,left,bottom,left);
    cout << "haversine: size "<<Utils::ConvertFixed(dwidth*(right-left)/1000.0,2) << "km x " << (dheight/1000.0) << "km" << endl;

    dwidth=CalcDist_Vincenty(lat_mid,left,lat_mid,left+1);
    dheight=CalcDist_Vincenty(top,left,bottom,left);
    cout << "vincenty:  size "<<Utils::ConvertFixed(dwidth*(right-left)/1000.0,2) << "km x " << (dheight/1000.0) << "km "<<endl;
    cout << endl;*/

    #endif
    return 0;
  } else return 1;
}

void Projection::SetProjection(double ttop,double tleft,double tbottom,double tright)
{
  top=ttop;left=tleft;bottom=tbottom;right=tright;
}

// calculate shortest distance on a spehre using haversine formula
double Projection::CalcDist_Haversine(double lat1,double long1,double lat2,double long2,double precision)
{
  const double r=6371000.785; // mean earth radius
  double phi1=Utils::ToRadians(lat1);
  double phi2=Utils::ToRadians(lat2);
  double dphi=Utils::ToRadians(lat2-lat1);
  double dlambda=Utils::ToRadians(long2-long1);

  double t1=sin(dphi/2.);
  double t2=sin(dlambda/2.);
  double a=t1*t1+cos(phi1)*cos(phi2)*t2*t2;
  double c=2.*atan2(sqrt(a),sqrt(1.-a));
  double s=r*c;
  if (precision>0.)  return roundf(s * precision) / precision;
  else return s;
}

// calculate shortest distance on an oblate ellipsoid using vencenty formula
double Projection::CalcDist_Vincenty(double lat1,double long1,double lat2,double long2,double precision)
{
  double phi1=Utils::ToRadians(lat1);
  double phi2=Utils::ToRadians(lat2);
  double U1=atan((1.0-WGS84_f)*tan(phi1));
  double U2=atan((1.0-WGS84_f)*tan(phi2));
  double L=Utils::ToRadians(long2-long1);
  double sinU1=sin(U1);
  double sinU2=sin(U2);
  double cosU1=cos(U1);
  double cosU2=cos(U2);


  int iter_limit=256;
  double t1,t2;
  double sinLambda,cosLambda,sinAlpha,sinSigma,cosSigma,cosSqAlpha,lambdaP,cos2SigmaM,sigma,C;
  double lambda=L;
  do {
    sinLambda=sin(lambda);
    cosLambda=cos(lambda);

    t1=cosU2*sinLambda;
    t2=cosU1*sinU2-sinU1*cosU2*cosLambda;
    sinSigma=sqrt(t1*t1+t2*t2);

    if (fabs(sinSigma)<std::numeric_limits<double>::epsilon()) { // co-incident points
       return 0;
    }

    cosSigma=sinU1*sinU2+cosU1*cosU2*cosLambda;

    sigma=atan2(sinSigma,cosSigma);
    sinAlpha=cosU1*cosU2*sinLambda/sinSigma;

    cosSqAlpha=1.0-sinAlpha*sinAlpha;

    cos2SigmaM=cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

    if (std::isnan(cos2SigmaM)) cos2SigmaM=0.0; // equatorial line, cosSqAlpha=0

    C = WGS84_f / 16.0 * cosSqAlpha * (4.0 + WGS84_f * (4.0 - 3.0 * cosSqAlpha));
    lambdaP = lambda;
    lambda  = L + (1 - C) * WGS84_f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM)));
  } while (fabs(lambda-lambdaP)>EPS && --iter_limit>0);

  if (iter_limit==0) return NAN;  // failed to converge

  double uSq = cosSqAlpha * (WGS84_a * WGS84_a - WGS84_b * WGS84_b) / (WGS84_b * WGS84_b);
  double A = 1. + uSq / 16384. * (4096. + uSq * (-768. + uSq * (320. - 175. * uSq)));
  double B = uSq / 1024. * (256. + uSq * (-128. + uSq * (74. - 47. * uSq)));
  double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4. * (cosSigma * (-1. + 2. * cos2SigmaM * cos2SigmaM) - B / 6. * cos2SigmaM * (-3. + 4. * sinSigma * sinSigma) * (-3. + 4. * cos2SigmaM * cos2SigmaM)));
  double s = WGS84_b * A * (sigma - deltaSigma);
  if (precision>0.)  return roundf(s * precision) / precision;
  else return s;
}

// get the geo-reference for a certain pixel
void Projection::GetLatLong(int w,int h,double &geo_width,double &geo_height)
{
  geo_width=left+(double)w*cellsize; // longitude
  geo_height=top-(double)h*cellsize; // latitude
}

double Projection::GetLong(int w)
{
  return left+(double)w*cellsize;
}

void Projection::SetDummyInterpolation(int len)
{
  inter_matrix.resize(1);
  //dy=height;
  inter_cell icell;
  icell.pixel_width_top=icell.pixel_width_bottom=(double)len;
  icell.pixel_height=(double)len;
  icell.pixel_area=(double)(len*len) ;
  inter_matrix[0]=icell;
}

int Projection::GetPPosLong(double dLong)
{
  return round(dLong/cellsize);
}

// calculate the distance in [m] of dLong at a certain latitude
// Peter Osborne 2013, "The Mercator Projections"
double Projection::CalcDeltaLong(double lat,double dLong=1.0)
{
  double phi=Utils::ToRadians(lat);
  double sinphi=sin(phi);
  double v_phi=WGS84_a/sqrt(1.-WGS84_e2*sinphi*sinphi); // radius of curvature along a parallel
  return v_phi*cos(phi)*Utils::ToRadians(dLong);

  //double w=atan(WGS84_b/WGS84_a*tan(phi)); // alternative formula
  //return WGS84_a*cos(w)*Utils::ToRadians(1);
}

// calculate the distance in [m] on the meridian for a certain latitude
// Peter Osborne 2013, "The Mercator Projections"
double Projection::CalcDeltaLat(double lat)
{
  const double phi=Utils::ToRadians(lat);
  const double A0=6367449.156;
  const double A2=16038.509;
  const double A4=16.833;
  const double A6=0.022;
  const double A8=0.00003;
  return A0*phi-A2*sin(2*phi)+A4*sin(4*phi)-A6*sin(6*phi)+A8*sin(8*phi);
}

double Projection::CalcDeltaLat(double lat1,double lat2)
{
  return CalcDeltaLat(lat2)-CalcDeltaLat(lat1);
}

/* calculate the surface area of a latitude-longitude quadrangle
// vincenty formula (shortest distance), approximate the area with a trapezoid
// assumes epsg:4326 projection
// useful for small areas: vincenty goes along a geodesic not along a circle of latitudes for the width
inter_cell Projection::CalcArea_geodesic(double top,double left,double bottom,double right)
{
  double sheight=CalcDist_Vincenty(top,left,bottom,left); // calculate vertical distance of left edge
  double top_width=CalcDist_Vincenty(top,left,top,right); // calculate horizontal distances
  double bottom_width=CalcDist_Vincenty(bottom,left,bottom,right);

  double t=(bottom_width-top_width)/2;
  double h=sqrt(sheight*sheight-t*t);
  inter_cell cell;
  cell.pixel_width_top=top_width; // calculate lengths in pixels
  cell.pixel_width_bottom=bottom_width;
  cell.pixel_height=sheight;
  cell.pixel_area=(top_width+bottom_width)/2.*h;
  return cell;
}*/

// calculate the surface area of latitude-longitude rectangle using trapezoid approximation
// assumes epsg:4326 projection
// usefull for small areas
void Projection::CalcArea_trapezoid_parallel(double top,double left,double bottom,double right,double &sheight,double &top_width,double &bottom_width,double &area)
{
  sheight=CalcDeltaLat(bottom,top);
  top_width=CalcDeltaLong(top,right-left);
  bottom_width=CalcDeltaLong(bottom,right-left);

  double t=(bottom_width-top_width);
  double h=sqrt(sheight*sheight-t*t/4);
  area=(top_width+bottom_width)/2.*h;
}

// calculate the surface area of a pixel at top latitude top
inter_cell Projection::CalcPixelArea_trapezoid(double top)
{
  inter_cell cell;
  CalcArea_trapezoid_parallel(top,0,top-cellsize,cellsize,cell.pixel_height,cell.pixel_width_top,cell.pixel_width_bottom,cell.pixel_area);
  return cell;
}

// calculate the surface area of latitude-longitude rectangle from top to zero using exact formula
double Projection::CalcArea_rectangle(double top,double left,double right)
{
  double phi=Utils::ToRadians(top);
  double sinphi=sin(phi);

  double t=WGS84_b/WGS84_a;
  double e=sqrt(1.0-t*t);
  double zm=1-e*sinphi;
  double zp=1+e*sinphi;

  double area=M_PI*WGS84_b*WGS84_b*(atanh(e*sinphi)/e+sinphi/(zp*zm));
  return area*((right-left)/360);
}

inter_cell Projection::CalcPixelArea_exact(double top)
{
  double a1=CalcArea_rectangle(top,0,cellsize);
  double a2=CalcArea_rectangle(top-cellsize,0,cellsize);

  inter_cell cell;
  cell.pixel_height=CalcDeltaLat(top-cellsize,top);
  cell.pixel_width_top=CalcDeltaLong(top,cellsize);
  cell.pixel_width_bottom=CalcDeltaLong(top-cellsize,cellsize);
  cell.pixel_area=a1-a2;
  return cell;
}

void Projection::GenerateInterpolation()
{
  inter_matrix.resize(height);
  //double len_test=CalcDist_Vincenty(top,left,bottom,left); // calculate len from top to bottom
  double len_test=CalcDeltaLat(bottom,top); // calculate len from top to bottom
  double l1,t1;
  //int dy=ceil( (double)height/(double)ncells);
  cout << "proj: calculating pixel-area for " << height << " latitude(s)\n";

  #if 0
  cout << "proj: 1degree of longitude at 45 degree lat (geodesic) = " << Utils::ConvertFixed(CalcDist_Vincenty(45,0,45,1)/1000,4) << "km" << endl;
  cout << "proj: 1degree of latitude  at 45 degree lat (geodesic) = " << Utils::ConvertFixed(CalcDist_Vincenty(45,0,46,0)/1000,4) << "km" << endl;

  cout << "proj: 1degree of longitude at 45 degree lat (parallel) = " << Utils::ConvertFixed(CalcDeltaLong(45)/1000,4) << "km" << endl;
  cout << "proj: 1degree of latitude  at 45 degree lat (parallel) = " << Utils::ConvertFixed(CalcDeltaLat(45,46)/1000,4) << "km" << endl;
  #endif

  double len_calc=0.;
  double mean_pixel_area=0;
  double latitude=top;

  double min_width=std::numeric_limits<double>::max();
  double min_height=std::numeric_limits<double>::max();
  double max_width=0.;
  double max_height=0.;

  double edge_distance=50.0; // length of non-forest distance to declare a forest-pixel as edge
  for (int i=0;i<height;i++)
  {
    GetLatLong(0,i,l1,t1); // slightly better precicion than using (latitude,latitude-cellsize)


    inter_cell cell1=CalcPixelArea_trapezoid(t1); // calc using trapezoid
    inter_cell cell2=CalcPixelArea_exact(t1);

    mean_pixel_area+=cell1.pixel_area;

    if (fabs(cell2.pixel_area-cell1.pixel_area)>1E-5) std::cerr << "proj: pixel areas mismatch" << endl;

    double mean_pixel_width=(cell2.pixel_width_top+cell2.pixel_width_top)/2.; // doesn't really matter
    double mean_pixel_height=cell2.pixel_height;
    cell2.npixel_horiz=ceil(edge_distance/mean_pixel_width); // calculate number of pixels needed to match "edge_distance"
    cell2.npixel_vert=ceil(edge_distance/mean_pixel_height);
    if (cell2.npixel_vert>max_npixel_vert) max_npixel_vert=cell2.npixel_vert;

    if (fabs(t1)<1E-10 || (i==0) || (i==height-1)) {
      std::cout << "proj: pixels for edge distance of " << Utils::ConvertFixed(edge_distance,2) << "m at lat " << std::setw(5) << Utils::ConvertFixed(t1,1) << "rad: " << cell2.npixel_horiz << "x" << cell2.npixel_vert << std::endl;
    }

    max_width=std::max(max_width,cell2.pixel_width_top);
    max_width=std::max(max_width,cell2.pixel_width_bottom);
    min_width=std::min(min_width,cell2.pixel_width_bottom);
    min_width=std::min(min_width,cell2.pixel_width_top);
    max_height=std::max(max_height,cell2.pixel_height);
    min_height=std::min(min_height,cell2.pixel_height);

    inter_matrix[i]=cell2;
    len_calc+=cell2.pixel_height;
    latitude-=cellsize;
  }
  //cout << latitude << " " << bottom << endl;
  double diff=fabs(len_test-len_calc);
  mean_pixel_area/=(double)height;
  // total difference should not be larger than 0.5m
  std::cout << "proj: mean-pixel area: " << Utils::ConvertFixed(mean_pixel_area,2) << " m^2" << std::endl;
  std::cout << "min width: " << Utils::ConvertFixed(min_width,2) << "m, height: " << Utils::ConvertFixed(min_height,2) << "m" << std::endl;
  std::cout << "max width: " << Utils::ConvertFixed(max_width,2) << "m, height: " << Utils::ConvertFixed(max_height,2) << "m" << std::endl;
  if (diff>0.) std::cout << "proj: calculated distances do not match (" << diff << " m)" << endl;
  std::cout << endl;
}
