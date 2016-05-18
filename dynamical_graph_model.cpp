#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <algorithm>
#include <boost/cstdint.hpp>
#include <boost/random.hpp>

//================================================
class Species {
public:
  Species(uint32_t immigrationTime) : m_immigrationTime(immigrationTime), m_fitness(0.0) {};

  struct comp {
    bool operator() (const Species* lhs, const Species* rhs) const {
      return lhs->m_immigrationTime < rhs->m_immigrationTime;
    }
  };
  typedef std::map<Species*, double, comp> SpeciesMap;
  typedef std::set<Species*, comp> SpeciesSet;

  void MakeInteractionWith(Species* other, double coefficient);
  void DeleteInteractions( SpeciesSet& dyingSpecies); // delete all the interactions
  // insert species that can go extinct into dyingSpecies
  bool IsGoingExtinct();
  std::string ToString();
  double Fitness() { return m_fitness; }
  uint32_t ImmigrationTime() { return m_immigrationTime; }
  double LocalCC() const;
private:
  const uint32_t m_immigrationTime;
  SpeciesMap m_incomingInteractions;
  SpeciesSet m_outgoingInteractions;
  double m_fitness;
};

void Species::MakeInteractionWith(Species* other, double coefficient) {
  m_incomingInteractions.insert( std::make_pair(other, coefficient) );
  m_fitness += coefficient;
  other->m_outgoingInteractions.insert(this);
}

void Species::DeleteInteractions( SpeciesSet& dyingSpecies) {
  for( auto pair : m_incomingInteractions ) {
    pair.first->m_outgoingInteractions.erase(this);
  }
  for( auto other : m_outgoingInteractions ) {
    auto found = other->m_incomingInteractions.find(this);
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
  for( auto pair : m_incomingInteractions ) {
    oss << "  " << pair.first << ", " << pair.second << std::endl;
  }
  return oss.str();
}

double Species::LocalCC() const {
  SpeciesSet neighbors = m_outgoingInteractions;
  for( auto pair : m_incomingInteractions ) {
    neighbors.insert( pair.first );
  }
  size_t k = neighbors.size();
  if( k <= 1 ) { return 0.0; }

  size_t num_triads = 0;
  for( auto pSpecies : neighbors ) {
    SpeciesSet nn = pSpecies->m_outgoingInteractions;
    std::vector<Species*> common;
    std::set_intersection( neighbors.begin(), neighbors.end(), nn.begin(), nn.end(), std::back_inserter(common), Species::comp() );
    num_triads += common.size();
  }

  return ( static_cast<double>(num_triads) / (k*(k-1)) );
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

  Species::SpeciesSet m_species;
  Species::SpeciesSet m_dyingSpecies;
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
  double CC();
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
      fout << m_currentTime << ' ' << Diversity() << ' ' << CC() << std::endl;
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
  for( auto pSpecies : m_species ) {
    if( uniform(*pRnd) < m_connectance ) {
      s->MakeInteractionWith( pSpecies, nd(*pRnd) );
    }
    if( uniform(*pRnd) < m_connectance ) {
      pSpecies->MakeInteractionWith( s, nd(*pRnd) );
      if( pSpecies->IsGoingExtinct() ) { m_dyingSpecies.insert(pSpecies); }
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
  // Find minimum
  Species* pMinSpecies = NULL;
  double min = 0.0;
  for( auto pSpecies : m_dyingSpecies ) {
    double f = pSpecies->Fitness();
    if( min > f ) {
      min = f;
      pMinSpecies = pSpecies;
    }
  }
  if( pMinSpecies != NULL ) {
    Extinct(pMinSpecies);
    m_dyingSpecies.erase(pMinSpecies);
  }
  return (pMinSpecies != NULL);
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

  for( auto pair : lifetime_histo ) {
    fout << pair.first << ' ' << pair.second << std::endl;
  }
  fout.close();
}
//================================================
void DynamicalGraph::ExtinctionSizeHistoOutput( const char* filename) {
  std::ofstream fout(filename);

  for( auto pair : extinction_histo ) {
    fout << pair.first << ' ' << pair.second << std::endl;
  }
  fout.close();
}
//================================================
void DynamicalGraph::DiversityHistoOutput( const char* filename) {
  std::ofstream fout(filename);

  for( auto pair : diversity_histo ) {
    fout << pair.first << ' ' << pair.second << std::endl;
  }
  fout.close();
}
//================================================
int DynamicalGraph::Diversity() {
  return m_species.size();
}

//================================================
double DynamicalGraph::CC() {
  if( m_species.empty() ) { return 0.0; }
  double total = 0.0;
  for( auto s : m_species ) {
    double d = s->LocalCC();
    total += d;
  }
  return total / m_species.size();
}
