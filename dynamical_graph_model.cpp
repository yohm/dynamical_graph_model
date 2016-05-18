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
class Histogram {
public:
  Histogram() {};
  std::map<uint32_t, uint32_t> histo;

  void Add( uint32_t x ) {
    if( histo.find(x) != histo.end() ) {
      histo[x] += 1;
    }
    else {
      histo[x] = 1;
    }
  }
  void Print( std::ostream& out ) const {
    for( auto pair : histo ) {
      out << pair.first << ' ' << pair.second << std::endl;
    }
  }
  void Clear() { histo.clear(); }
};

//==============================================

class DynamicalGraph {
public:
  DynamicalGraph( uint64_t seed, double connectance);
  ~DynamicalGraph() {};
  void Run( uint32_t t_init, uint32_t t_measure);
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

  Histogram extinction_histo;
  Histogram diversity_histo;
  Histogram lifetime_histo;
  void ClearHisto() {
    extinction_histo.Clear();
    diversity_histo.Clear();
    lifetime_histo.Clear();
  };

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
void DynamicalGraph::Run(uint32_t t_init, uint32_t t_measure) {
  std::ofstream fout("timeseries.dat");

  for( ; m_currentTime<t_init; m_currentTime++) {
    Update();
    if( m_currentTime%1024 == 0 ) {
      std::cerr << "t : " << m_currentTime << std::endl;
      fout << m_currentTime << ' ' << Diversity() << ' ' << CC() << std::endl;
    }
  }

  ClearHisto();

  for( ; m_currentTime<t_init+t_measure; m_currentTime++) {
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
  extinction_histo.Add( extinction_size );

  diversity_histo.Add( current_diversity );
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

  lifetime_histo.Add( lifetime );

  delete s;
}

//================================================
void DynamicalGraph::LifetimeHistoOutput( const char* filename) {
  std::ofstream fout(filename);
  lifetime_histo.Print( fout );
  fout.close();
}
//================================================
void DynamicalGraph::ExtinctionSizeHistoOutput( const char* filename) {
  std::ofstream fout(filename);
  extinction_histo.Print( fout );
  fout.close();
}
//================================================
void DynamicalGraph::DiversityHistoOutput( const char* filename) {
  std::ofstream fout(filename);
  diversity_histo.Print( fout );
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
