#ifndef GLOBAL_H
#define GLOBAL_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <list>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <sys/time.h>

using namespace std;

// enable to compile against libtiff
#define TIFF_SUPPORT

#define LISA_VERSION "Large Image Spatial Analysis v1.00b (c) 2014-2016 - Sebastian Lehmann"

#define BOOL(x) (!(!(x)))

#define BitSet(arg,posn) ((arg) | (1L << (posn)))
#define BitClear(arg,posn) ((arg) & ~(1L << (posn)))
#define BitTest(arg,posn) BOOL((arg) & (1L << (posn)))
#define BitFlip(arg,posn) ((arg) ^ (1L << (posn)))

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define EPS (1E-12)


#endif // GLOBAL_H
