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
  rvalue_(3,0),
  sector_(0),
  region_(0),
  tpdunit_(0),
  firstbin_out_(0),
  firstbin_in_(0),
  diffmax_out_(0),
  diffmax_in_(0),
  rzeff_out_(0),
  rzeff_in_(0)
{
}

// Destructor definition
L1StubTriplet::~L1StubTriplet() {
  // Cleanup logic (if any) goes here
}