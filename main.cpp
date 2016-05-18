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

  if(argc != 5) {
    std::cerr << "Usage : ./network-full.out <connectance> <t_init> <t_measure> <seed>" << std::endl;
    exit(1);
  }

  const double connectance = boost::lexical_cast<double>(argv[1]);
  const uint32_t t_init = boost::lexical_cast<uint32_t>(argv[2]);
  const uint32_t t_measure = boost::lexical_cast<uint32_t>(argv[3]);
  const uint64_t seed = boost::lexical_cast<uint64_t>(argv[4]);
  
  std::cerr << "Lists of given parameters are as follows:" << std::endl
            << "connectance:\t" << connectance << std::endl
            << "t_init:\t" << t_init << std::endl
            << "t_measure:\t" << t_measure << std::endl
            << "seed:\t" << seed << std::endl;

  //ofstreams
  DynamicalGraph eco( seed, connectance);
  boost::timer t;
  eco.Run(t_init, t_measure);
  eco.LifetimeHistoOutput("lifetime.dat");
  eco.DiversityHistoOutput("diversity_histo.dat");
  eco.ExtinctionSizeHistoOutput("extinction_histo.dat");
  std::ofstream json("_output.json");
  json << "{ \"diversity\": " << eco.AverageDiversity()
       << ", \"CC\": " << eco.AverageCC() << " }" << std::endl;
  json.close();
  std::cerr << "elapsed time : " << t.elapsed() << std::endl;

  return 0;
}


//----------------------end of the program-------------------------------
