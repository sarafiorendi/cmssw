#include "DataFormats/L1TrackTrigger/interface/L1StubTriplet.h"

/// Default Constructor
/// NOTE: to be used with setSomething(...) methods
L1StubTriplet::L1StubTriplet() :
  rapprox_(3,0),
  zapprox_(3,0),
  bend_(3,0),
  phi_(3,0),
  rzbin_(3,0),
  index_(3,0),
  layerdisk_(3,0),
  rvalue_(3,0),
  sector_(0),
  region_(0),
  tpdunit_(0),
  firstbin_out_(0),
  firstbin_in_(0),
  firstbin_out_pair_(0),
  diffmax_out_(0),
  diffmax_in_(0),
  diffmax_out_pair_(0),
  rzeff_out_(0),
  rzeff_in_(0),
  start_in_(0),
  start_out_pair_(0)
{
}

// L1StubTriplet::L1StubTriplet(const trklet::Stub* inner,
//                              const trklet::Stub* middle,
//                              const trklet::Stub* outer,
//                              unsigned int sector,
//                              int region,
//                              int tpdunit)
//     : rapprox_{inner ->rapprox(), 
//                middle->rapprox(), 
//                outer ->rapprox()},
//       zapprox_{inner ->zapprox(), 
//                middle->zapprox(), 
//                outer ->zapprox()},
//       bend_   {inner ->bend().value(), 
//                middle->bend().value(), 
//                outer ->bend().value()},
//       phi_    {inner ->phiapprox(0., 0), 
//                middle->phiapprox(0., 0), 
//                outer ->phiapprox(0., 0)},
//       index_  {inner ->stubindex().value(), 
//                middle->stubindex().value(), 
//                outer ->stubindex().value()},
//       layerdisk_{inner ->layerdisk(), 
//                  middle->layerdisk(), 
//                  outer ->layerdisk()},
//       sector_(sector),
//       region_(region),
//       tpdunit_(tpdunit) 
// {}

// Destructor definition
L1StubTriplet::~L1StubTriplet() {
  // Cleanup logic (if any) goes here
}