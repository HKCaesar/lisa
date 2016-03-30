#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include "../global.h"

// time measuring via C++ 11 chrono
class Timer {
  typedef std::chrono::high_resolution_clock Clock;
  //typedef std::chrono::milliseconds TimeT;
  typedef std::chrono::duration<double, std::ratio<1>> durationS;
  typedef std::chrono::duration<double, std::milli> durationMS;
  typedef Clock::time_point Timepoint;
  public:
    void Start() {tstart=Clock::now();};
    void Stop() {tstop=Clock::now();};
    double ElapsedMS() {
      durationMS elapsed=tstop-tstart;
      return elapsed.count();
      //cout << std::chrono::duration_cast<double,TimeT>(tstop - tstart).count() << endl;
      //return std::chrono::duration_cast<TimeT>(tstop - tstart).count();
    };
    double ElapsedS() {
      durationS elapsed=tstop-tstart;
      return elapsed.count();
    };
  private:
    Timepoint tstart,tstop;
};
#endif // TIMER_H
