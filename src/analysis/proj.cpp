#include "proj.h"

int Projection::ReadCoordinateFile(std::string &fname)
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
    cout << endl;
    CalculateGeoDelta();
    double dwidth=CalcDist_Haversine(top+(bottom-top)/2,left,top+(bottom-top)/2,right);
    double dheight=CalcDist_Haversine(top,left,bottom,left);
    double pixel_area=(dwidth/(double)width)*(dheight/(double)height);
    cout << "haversine: size "<<(dwidth/1000.0) << "km x " << (dheight/1000.0) << "km, " << " pixel area: " << pixel_area << " m^2" << endl;

    dwidth=CalcDist_Vincenty(top+(bottom-top)/2,left,top+(bottom-top)/2,right);
    dheight=CalcDist_Vincenty(top,left,bottom,left);
    pixel_area=(dwidth/(double)width)*(dheight/(double)height);
    mean_pixelarea=pixel_area;
    cout << "vincenty:  size "<<(dwidth/1000.0) << "km x " << (dheight/1000.0) << "km, " << " pixel area: " << pixel_area << " m^2" << endl;
    cout << endl;
    return 0;
  } else return 1;
}

double Projection::CalcDist_Haversine(double lat1,double long1,double lat2,double long2)
{
  double r=6371000.785; // mean earth radius
  double phi1=Utils::ToRadians(lat1);
  double phi2=Utils::ToRadians(lat2);
  double dphi=Utils::ToRadians(lat2-lat1);
  double dlambda=Utils::ToRadians(long2-long1);
  double t1=sin(dphi/2.);
  double t2=sin(dlambda/2.);
  double a=t1*t1+cos(phi1)*cos(phi2)*t2*t2;
  double c=2.*atan2(sqrt(a),sqrt(1.-a));
  return r*c;
}

double Projection::CalcDist_Vincenty(double lat1,double long1,double lat2,double long2)
{
  double a=6378137.0; //length of semi-major axis
  double f=1.0/298.257223563; //flattening of the ellipsoid
  double b=(1.0-f)*a; //length of semi-minor axis of the ellipsoid

  double phi1=Utils::ToRadians(lat1);
  double phi2=Utils::ToRadians(lat2);
  double U1=atan((1.0-f)*tan(phi1));
  double U2=atan((1.0-f)*tan(phi2));
  double L=Utils::ToRadians(long2-long1);
  double sinU1=sin(U1);
  double sinU2=sin(U2);
  double cosU1=cos(U1);
  double cosU2=cos(U2);


  int iter_limit=100;
  double t1,t2;
  double sinLambda,cosLambda,sinAlpha,sinSigma,cosSigma,cosSqAlpha,lambdaP,cos2SigmaM,sigma,C;
  double lambda=L;
  do {
    sinLambda=sin(lambda);
    cosLambda=cos(lambda);

    t1=cosU2*sinLambda;
    t2=cosU1*sinU2-sinU1*cosU2*cosLambda;
    sinSigma=sqrt(t1*t1+t2*t2);
    cosSigma=sinU1*sinU2+cosU1*cosU2*cosLambda;

    sigma=atan2(sinSigma,cosSigma);
    sinAlpha=cosU1*cosU2*sinLambda/sinSigma;

    cosSqAlpha=1.0-sinAlpha*sinAlpha;

    if (cosSqAlpha>EPS) cos2SigmaM=cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
    else cos2SigmaM=0.0;

    C = f / 16.0 * cosSqAlpha * (4.0 + f * (4.0 - 3.0 * cosSqAlpha));
    lambdaP = lambda;
    lambda  = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM)));
  } while (fabs(lambda-lambdaP)>EPS && --iter_limit>0);

  if (iter_limit==0) return -1;

  double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
  double A = 1. + uSq / 16384. * (4096. + uSq * (-768. + uSq * (320. - 175. * uSq)));
  double B = uSq / 1024. * (256. + uSq * (-128. + uSq * (74. - 47. * uSq)));
  double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4. * (cosSigma * (-1. + 2. * cos2SigmaM * cos2SigmaM) - B / 6. * cos2SigmaM * (-3. + 4. * sinSigma * sinSigma) * (-3. + 4. * cos2SigmaM * cos2SigmaM)));
  double s = b * A * (sigma - deltaSigma);
  return s;
  //return roundf(s * 100) / 100;
}

void Projection::GetLatLong(int w,int h,double &geo_width,double &geo_height)
{
  geo_width=left+(double)w*cellsize; // longitude
  geo_height=top-(double)h*cellsize; // latitude
}

void Projection::SetDummyInterpolation(int len)
{
  inter_matrix.resize(1);
  dy=height;
  inter_cell icell;
  icell.pixel_width=(double)len;
  icell.pixel_height=(double)len;
  icell.pixel_area=(double)(len*len) ;
  inter_matrix[0]=icell;
}

void Projection::GenerateInterpolation(int ncells)
{
  inter_matrix.resize(ncells);
  double len_test=CalcDist_Vincenty(top,left,bottom,left); // calculate len from top to bottom
  double l1,t1,l2,t2;
  dy=ceil( (double)height/(double)ncells);
  cout << "calculating pixel-area for " << ncells << " latitude(s), ";
  cout << "clipping " << ncells*dy-height << " pixels" << endl;
  double len_calc=0.;
  for (int i=0;i<ncells;i++)
  {
    int y1=i*dy;
    int y2=(i+1)*dy;
    if (y2>height) {y2=height;};
    GetLatLong(0,y1,l1,t1);
    GetLatLong(0,y2,l2,t2);
    double len_verti=CalcDist_Vincenty(t1,left,t2,left);
    double len_horiz=CalcDist_Vincenty(t1+(t2-t1)/2.,left,t1+(t2-t1)/2.,right);

    double len_pixel_verti=len_verti/double(y2-y1);
    double len_pixel_horiz=len_horiz/double(width);
    double pixel_area=len_pixel_horiz*len_pixel_verti;

    inter_cell icell;
    icell.pixel_width=len_pixel_horiz;
    icell.pixel_height=len_pixel_verti;
    icell.pixel_area=pixel_area;
    inter_matrix[i]=icell;
    //cout << fixed << setprecision(2) << t1+(t2-t1)/2. << " degree: " << icell.pixel_width << " m x" << icell.pixel_height << " m=" << icell.pixel_area << " m^2" << endl;

    len_calc+=len_verti;
  }
  double diff=fabs(len_test-len_calc);
  // total difference should not be larger than 0.5m
  if (diff>0.) cout << "warning: calculated distances do not match (" << diff << " m)" << endl;
  cout << endl;
}
