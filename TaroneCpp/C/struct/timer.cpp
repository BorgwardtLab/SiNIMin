/*
* @Author: Anja Gumpinger
* @Date:   2018-11-16 13:27:55
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-10 15:25:16
*/

#ifndef _timer_cpp_
#define _timer_cpp_

#include <fstream>
#include <time.h>
#include <sys/resource.h>


class Timer
{
  public:
    double io;
    double meta;
    double proc_edge;
    double total;
    int n_cores;

    // Constructor.
    Timer(int n);

    void write_profiling(std::string filename);
};


Timer::Timer(int n)
{
  io = 0;
  meta = 0;
  proc_edge = 0;
  total = 0;
  n_cores = n;
}


/*
  Function to write profiling file.
*/
void Timer::write_profiling(std::string filename)
{

  std::ofstream file(filename);

  file << "Total time: " << total << "secs." << std::endl;
  file << "\t * Initialization: " << io << "secs." << std::endl;
  file << "\t * Processing edges: " << proc_edge << "secs.";
  file << std::endl << std::endl;

  file << "Number of cores: " << n_cores << std::endl;

  file.close();
}


/*
  Measure running time (Felipe Llinares).
*/
double measureTime(){
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}


#endif
