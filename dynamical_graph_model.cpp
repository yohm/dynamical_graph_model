#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <set>
#include <map>
#include <boost/cstdint.hpp>
#include <boost/timer.hpp>
#include <boost/random.hpp>

//================================================
class Species {
public:
  Species(uint32_t immigrationTime) : m_immigrationTime(immigrationTime), m_fitness(0.0) {};
  void MakeInteractionWith(Species* other, double coefficient);
  void DeleteInteractions( std::set<Species*> &dyingSpecies); // delete all the interactions
  // insert species that can go extinct into dyingSpecies
  bool IsGoingExtinct();
  std::string ToString();
  double Fitness() { return m_fitness; }
  uint32_t ImmigrationTime() { return m_immigrationTime; }
private:
  const uint32_t m_immigrationTime;
  std::map<Species*, double> m_incomingInteractions;
  std::set<Species*> m_outgoingInteractions;
  double m_fitness;
};

void Species::MakeInteractionWith(Species* other, double coefficient) {
  m_incomingInteractions.insert( std::make_pair(other, coefficient) );
  m_fitness += coefficient;
  other->m_outgoingInteractions.insert(this);
}

void Species::DeleteInteractions( std::set<Species*>& dyingSpecies) {
  for(std::map<Species*, double>::iterator it = m_incomingInteractions.begin();
      it != m_incomingInteractions.end(); ++it) {
    it->first->m_outgoingInteractions.erase(this);
  }
  for(std::set<Species*>::iterator it = m_outgoingInteractions.begin();
      it != m_outgoingInteractions.end(); ++it) {
    Species* other = *it;
    std::map<Species*, double>::iterator found = other->m_incomingInteractions.find(this);
    other->m_fitness -= found->second;
    other->m_incomingInteractions.erase(found);

    if( other->IsGoingExtinct() ) {
      dyingSpecies.insert(other);
    }
  }
}

bool Species::IsGoingExtinct() {
  return (m_fitness < 0.0 && m_incomingInteractions.size() > 0);
}

std::string Species::ToString() {
  std::ostringstream oss;
  oss << "immigrationTime: " << m_immigrationTime << std::endl;
  oss << "incoming:" << std::endl;
  for( std::map<Species*, double>::iterator it = m_incomingInteractions.begin();
    it != m_incomingInteractions.end(); ++it) {
    oss << "  " << it->first << ", " << it->second << std::endl;
  }
  return oss.str();
}

//==============================================

class DynamicalGraph {
public:
  DynamicalGraph( uint64_t seed, double connectance);
  ~DynamicalGraph() {};
  void Run( uint32_t tmax);
  void LifetimeHistoOutput( const char* filename);
  void DiversityHistoOutput( const char* filename);
  void ExtinctionSizeHistoOutput( const char* filename);
private:
  const double m_connectance;

  std::set<Species*> m_species;
  std::set<Species*> m_dyingSpecies;
  uint32_t m_currentTime;

  boost::mt19937 * pRnd;
  void Update();

  std::map<uint32_t, uint32_t> extinction_histo;
  std::map<uint32_t, uint32_t> diversity_histo;
  std::map<uint32_t, uint32_t> lifetime_histo;

  void AddOneSpecies();
  void AddRandomInteractions(Species* new_species);
  void Selection();
  bool RemoveNegativeFitnessSpecies();
  void Extinct(Species* s);
  int Diversity();
};

//================================================
DynamicalGraph::DynamicalGraph(uint64_t seed, double t_connectance)
: m_connectance(t_connectance), m_currentTime(0) {
  pRnd = new boost::mt19937(seed);
}
//================================================
void DynamicalGraph::Run(uint32_t tmax) {
  std::ofstream fout("timeseries.dat");
  for(m_currentTime=0; m_currentTime<tmax; m_currentTime++) {
    Update();
    if( m_currentTime%1024 == 0 ) {
      std::cerr << "t : " << m_currentTime << std::endl;
      fout << m_currentTime << ' ' << Diversity() << std::endl;
    }
  }
  fout.close();
}
//================================================
void DynamicalGraph::Update() {
  int previous_diversity = Diversity();

  AddOneSpecies();
  Selection();

  int current_diversity = Diversity();
  int extinction_size = previous_diversity + 1 - current_diversity;

  // calculate histograms
  if( extinction_histo.find(extinction_size) != extinction_histo.end() ) {
    extinction_histo[extinction_size] += 1;
  } else {
    extinction_histo[extinction_size] = 1;
  }

  if( diversity_histo.find(current_diversity) != diversity_histo.end() ) {
    diversity_histo[current_diversity] += 1;
  } else {
    diversity_histo[current_diversity] = 1;
  }
}
//================================================
void DynamicalGraph::AddOneSpecies() {
  Species* new_species = new Species(m_currentTime);
  AddRandomInteractions(new_species);
  m_species.insert(new_species);
}

//=================================================
void DynamicalGraph::AddRandomInteractions(Species* s) {
  boost::random::uniform_01<> uniform;
  boost::random::normal_distribution<> nd;
  for(std::set<Species*>::iterator it = m_species.begin(); it != m_species.end(); ++it) {
    if( uniform(*pRnd) < m_connectance ) {
      s->MakeInteractionWith( *it, nd(*pRnd) );
    }
    if( uniform(*pRnd) < m_connectance ) {
      (*it)->MakeInteractionWith( s, nd(*pRnd) );
      if( (*it)->IsGoingExtinct() ) { m_dyingSpecies.insert(*it); }
    }
  }

  if( s->IsGoingExtinct() ) { m_dyingSpecies.insert(s); }
}

//================================================
void DynamicalGraph::Selection() {
  bool bNegativeExists = true;
  while(bNegativeExists) {
    bNegativeExists = RemoveNegativeFitnessSpecies();
  }

  m_dyingSpecies.clear();
}

//================================================
bool DynamicalGraph::RemoveNegativeFitnessSpecies() {
  // remove negative fitness species
  bool negative_exists = false;
  // Find minimum
  typedef std::set<Species*>::iterator SETIT;
  SETIT it_min = m_dyingSpecies.begin();
  double min = 0.0;
  for(SETIT it = m_dyingSpecies.begin(); it != m_dyingSpecies.end(); ++it) {
    double f = (*it)->Fitness();
    if( min > f ) {
      min = f;
      it_min = it;
      negative_exists = true;
    }
  }
  if( negative_exists ) {
    Extinct(*it_min);
    m_dyingSpecies.erase(it_min);
  }
  return negative_exists;
}
//================================================
void DynamicalGraph::Extinct(Species* s) {
  s->DeleteInteractions(m_dyingSpecies);
  m_species.erase(s);

  int32_t lifetime = m_currentTime - s->ImmigrationTime();

  if( lifetime_histo.find(lifetime) != lifetime_histo.end() ) {
    lifetime_histo[lifetime] += 1;
  } else {
    lifetime_histo[lifetime] = 1;
  }

  delete s;
}

//================================================
void DynamicalGraph::LifetimeHistoOutput( const char* filename) {
  std::ofstream fout(filename);

  for(std::map<uint32_t, uint32_t>::iterator it = lifetime_histo.begin();
      it != lifetime_histo.end(); ++it) {
    fout << it->first << ' ' << it->second << std::endl;
  }
  fout.close();
}
//================================================
void DynamicalGraph::ExtinctionSizeHistoOutput( const char* filename) {
  std::ofstream fout(filename);

  for(std::map<uint32_t, uint32_t>::iterator it = extinction_histo.begin();
      it != extinction_histo.end(); ++it) {
    fout << it->first << ' ' << it->second << std::endl;
  }
  fout.close();
}
//================================================
void DynamicalGraph::DiversityHistoOutput( const char* filename) {
  std::ofstream fout(filename);

  for(std::map<uint32_t, uint32_t>::iterator it = diversity_histo.begin();
      it != diversity_histo.end(); ++it) {
    fout << it->first << ' ' << it->second << std::endl;
  }
  fout.close();
}
//================================================
int DynamicalGraph::Diversity() {
  return m_species.size();
}
