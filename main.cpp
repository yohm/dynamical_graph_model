//#define DEBUG
//#define BASIC
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/timer.hpp>
#include "dynamical_graph_model.cpp"

int main (int argc, char* argv[]) {

  if(argc != 4) {
    std::cerr << "Usage : ./network-full.out <connectance> <tmax> <seed>" << std::endl;
    exit(1);
  }

  const double connectance = boost::lexical_cast<double>(argv[1]);
  const uint32_t tmax = boost::lexical_cast<uint32_t>(argv[2]);
  const uint64_t seed = boost::lexical_cast<uint64_t>(argv[3]);
  
  std::cerr << "Lists of given parameters are as follows:" << std::endl
            << "connectance:\t" << connectance << std::endl
            << "tmax:\t" << tmax << std::endl
            << "seed:\t" << seed << std::endl;

  //ofstreams
  DynamicalGraph eco( seed, connectance);
  boost::timer t;
  eco.Run(tmax);
  eco.LifetimeHistoOutput("lifetime.dat");
  eco.DiversityHistoOutput("diversity_histo.dat");
  eco.ExtinctionSizeHistoOutput("extinction_histo.dat");
  std::cerr << "elapsed time : " << t.elapsed() << std::endl;

  return 0;
}


//----------------------end of the program-------------------------------
