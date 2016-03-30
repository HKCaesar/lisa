#include "cmdoptions.h"
#include "..\utils.h"

void CmdOptions::getOpt(std::string &val)
{
  if (optfound && optvalue.length()) {
    val=optvalue;
  };
}

void CmdOptions::getOpt(int &val)
{
  if (optfound && optvalue.length()) {
     try {
       val=std::stoi(optvalue);
     } catch (std::invalid_argument&) {
     }
  };
}

void CmdOptions::getOpt(double &val)
{
  if (optfound && optvalue.length()) {
     try {
       val=std::stod(optvalue);
     } catch (std::invalid_argument&) {
     }
  };
}

void CmdOptions::getOpt(std::vector <double>&dtokens)
{
  std::vector <std::string> stokens;
  StringUtils::Tokenize(optvalue,stokens,",");
  for (size_t i=0;i<stokens.size();i++) {
      double val=std::stod(stokens[i]);
      dtokens.push_back(val);
  }
}

bool CmdOptions::searchOption(const std::string &sshort,const std::string &&sslong)
{
  std::string usshort=StringUtils::toupper(sshort);
  std::string usslong=StringUtils::toupper(sslong);
  optvalue="";
  optfound=false;
  int i=1;
  while (i<argc) {
    std::string sarg=StringUtils::toupper(std::string(argv[i]));
    std::string arglong(sarg);
    std::string argshort(sarg.substr(0,2));

    if (usshort.length() && argshort.compare(usshort)==0) {
      if (sarg.length()>2) optvalue=sarg.substr(2);
      optfound=true;
      return optfound;
    } else if (usslong.length() && arglong.compare(usslong)==0) {
      if (i<argc-1) optvalue=std::string(argv[i+1]);
      optfound=true;
      return optfound;
    }
    i++;
  }
  return optfound;
}

