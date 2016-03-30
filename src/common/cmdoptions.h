#ifndef CMDOPTIONS_H
#define CMDOPTIONS_H

#include "..\global.h"

class CmdOptions {
  public:
    CmdOptions(int argc,char *argv[]):argc(argc),argv(argv) {};
    bool searchOption(const std::string &sshort,const std::string &&sslong); // check if option exists
    void getOpt(std::string &val);
    void getOpt(int &val);
    void getOpt(double &val);
    void getOpt(std::vector <double>&dtokens);
private:
    int argc;
    char **argv;
    bool optfound;
    std::string optvalue;
};

#endif // CMDOPTIONS_
