#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>
#include "dynamical_graph_model.cpp"

int main (int argc, char* argv[]) {

  if(argc != 5) {
    std::cerr << "Usage : ./main.out <connectance> <t_init> <t_measure> <seed>" << std::endl;
    exit(1);
  }

  const double connectance = std::stod(argv[1]);
  const uint32_t t_init = std::stoul(argv[2]);
  const uint32_t t_measure = std::stoul(argv[3]);
  const uint64_t seed = std::stoull(argv[4]);
  
  std::cerr << "Lists of given parameters are as follows:" << std::endl
            << "connectance:\t" << connectance << std::endl
            << "t_init:\t" << t_init << std::endl
            << "t_measure:\t" << t_measure << std::endl
            << "seed:\t" << seed << std::endl;

  //ofstreams
  DynamicalGraph dg( seed, connectance);
  std::clock_t start = std::clock();
  dg.Run(t_init, t_measure);
  dg.LifetimeHistoOutput("lifetime.dat");
  dg.DiversityHistoOutput("diversity_histo.dat");
  dg.ExtinctionSizeHistoOutput("extinction_histo.dat");
  std::ofstream json("_output.json");
  json << "{ \"diversity\": " << dg.AverageDiversity()
       << ", \"link_density\": " << dg.AverageLinkDensity()
       << ", \"CC\": " << dg.AverageCC() << " }" << std::endl;
  json.close();
  std::clock_t end = std::clock();
  std::cerr << "elapsed time : " << static_cast<double>(end-start)/CLOCKS_PER_SEC << std::endl;

  return 0;
}


//----------------------end of the program-------------------------------
