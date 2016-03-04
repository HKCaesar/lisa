#ifndef CMDOPTIONS_H
#define CMDOPTIONS_H

#include "global.h"

class CmdOptions {
  public:
    CmdOptions(int argc,char *argv[]):argc(argc),argv(argv) {};
    bool SearchOption(const std::string &sshort,const std::string &&sslong); // check if option exists
    void getopt(std::string &val);
    void getopt(int &val);
    void getopt(double &val);
    void getopt(std::vector <double>&dtokens);
private:
    int argc;
    char **argv;
    bool optfound;
    std::string optvalue;
};

#endif // CMDOPTIONS_
