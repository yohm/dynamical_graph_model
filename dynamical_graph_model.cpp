#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <algorithm>
#include <boost/cstdint.hpp>
#include <boost/random.hpp>
#include <boost/function_output_iterator.hpp>

//================================================
class Species {
public:
  Species(uint32_t immigrationTime) :
    m_immigrationTime(immigrationTime), m_fitness(0.0), m_cacheReady(false) {};

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
  size_t InDegree() const;
  size_t OutDegree() const;
private:
  const uint32_t m_immigrationTime;
  SpeciesMap m_incomingInteractions;
  SpeciesSet m_outgoingInteractions;
  double m_fitness;
  std::vector<uint32_t> m_outLinkCache;
  bool m_cacheReady;
  void MakeOutLinkCache();
};

void Species::MakeInteractionWith(Species* other, double coefficient) {
  m_incomingInteractions.insert( std::make_pair(other, coefficient) );
  m_fitness += coefficient;
  other->m_outgoingInteractions.insert(this);
  other->m_cacheReady = false;
}

void Species::DeleteInteractions( SpeciesSet& dyingSpecies) {
  for( auto pair : m_incomingInteractions ) {
    pair.first->m_outgoingInteractions.erase(this);
    pair.first->m_cacheReady = false;
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

  std::vector<uint32_t> neighbor_ids;
  for( auto n: neighbors ) {
    n->MakeOutLinkCache();
    neighbor_ids.push_back( n->m_immigrationTime );
  }

  size_t num_triads = 0;
  auto counter = [&num_triads](uint32_t){ num_triads++; };
  for( auto pSpecies : neighbors ) {
    std::set_intersection( neighbor_ids.begin(), neighbor_ids.end(),
        pSpecies->m_outLinkCache.begin(), pSpecies->m_outLinkCache.end(),
        boost::make_function_output_iterator(counter) );
  }

  return ( static_cast<double>(num_triads) / (k*(k-1)) );
}

size_t Species::InDegree() const {
  return m_incomingInteractions.size();
}

size_t Species::OutDegree() const {
  return m_outgoingInteractions.size();
}

void Species::MakeOutLinkCache() {
  if( m_cacheReady == false ) {
    m_outLinkCache.clear();
    for( auto s: m_outgoingInteractions ) {
      uint32_t t = s->m_immigrationTime;
      m_outLinkCache.push_back(t);
    }
    m_cacheReady = true;
  }
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
  double AverageDiversity() const;
  double AverageLinkDensity() const;
  double AverageCC() const;
private:
  const double m_connectance;
  Species::SpeciesSet m_species;
  Species::SpeciesSet m_dyingSpecies;
  uint32_t m_currentTime;
  boost::mt19937 * pRnd;
  Histogram extinction_histo;
  Histogram diversity_histo;
  Histogram lifetime_histo;
  double m_CCsum;
  size_t m_CCcount;
  double m_densitySum;
  size_t m_densityCount;

  void Update();
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
  double LinkDensity();
  double CC();
};

//================================================
DynamicalGraph::DynamicalGraph(uint64_t seed, double connectance)
: m_connectance(connectance), m_currentTime(0), m_CCsum(0.0), m_CCcount(0),
  m_densitySum(0.0), m_densityCount(0) {
  pRnd = new boost::mt19937(seed);
}
//================================================
void DynamicalGraph::Run(uint32_t t_init, uint32_t t_measure) {
  std::ofstream fout("timeseries.dat");

  for( ; m_currentTime<t_init; m_currentTime++) {
    Update();
    if( m_currentTime%1024 == 0 ) {
      std::cerr << "t : " << m_currentTime << std::endl;
      fout << m_currentTime << ' ' << Diversity()
        << ' ' << LinkDensity()
        << ' ' << CC() << std::endl;
    }
  }

  ClearHisto();

  for( ; m_currentTime<t_init+t_measure; m_currentTime++) {
    Update();
    m_densitySum += LinkDensity();
    m_densityCount++;
    if( m_currentTime % 128 == 0 ) {
      m_CCsum += CC();
      m_CCcount++;
    }
    if( m_currentTime%1024 == 0 ) {
      std::cerr << "t : " << m_currentTime << std::endl;
      fout << m_currentTime << ' ' << Diversity()
        << ' ' << LinkDensity()
        << ' ' << CC() << std::endl;
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
  }
  return (pMinSpecies != NULL);
}
//================================================
void DynamicalGraph::Extinct(Species* s) {
  s->DeleteInteractions(m_dyingSpecies);
  m_species.erase(s);

  int32_t lifetime = m_currentTime - s->ImmigrationTime();

  lifetime_histo.Add( lifetime );
  m_dyingSpecies.erase(s);

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
//================================================
double DynamicalGraph::LinkDensity() {
  size_t total = 0;
  for( auto s : m_species ) {
    total += s->OutDegree();
  }
  size_t n = m_species.size();
  return static_cast<double>(total) / (n*(n-1));
}
//================================================
double DynamicalGraph::AverageDiversity() const {
  uint64_t sum = 0;
  uint64_t count = 0;
  for( auto key_freq : diversity_histo.histo ) {
    sum += key_freq.first * key_freq.second;
    count += key_freq.second;
  }
  return static_cast<double>(sum)/count;
}
//================================================
double DynamicalGraph::AverageLinkDensity() const {
  return m_densitySum / m_densityCount;
}
//================================================
double DynamicalGraph::AverageCC() const {
  return m_CCsum / m_CCcount;
}

