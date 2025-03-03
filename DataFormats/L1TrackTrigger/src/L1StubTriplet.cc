#include "DataFormats/L1TrackTrigger/interface/L1StubTriplet.h"

/// Default Constructor
/// NOTE: to be used with setSomething(...) methods
L1StubTriplet::L1StubTriplet() :
  rapprox_(3,0),
  zapprox_(3,0),
  bend_(3,0),
  rzbin_(3,0),
  index_(3,0),
  layerdisk_(3,0),
  sector_(0),
  region_(0),
  tpdunit_(0)
{
  /// Set default data members
//   rapprox_.clear();
//   zapprox_.clear();
//   bend_.clear();
//   rzbin_.clear();
//   index_.clear();
//   layerdisk_.clear();
  
//   rapprox_.reserve(3); 
//   zapprox_.reserve(3); 
//   bend_.reserve(3); 
//   rzbin_.reserve(3); 
//   index_.reserve(3); 
//   layerdisk_.reserve(3); 

}

// Destructor definition
L1StubTriplet::~L1StubTriplet() {
  // Cleanup logic (if any) goes here
}