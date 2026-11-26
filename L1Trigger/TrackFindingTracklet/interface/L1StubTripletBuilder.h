#ifndef L1TRIGGER_TRACKFINDINGTRACKLET_L1STUBTRIPLETBUILDER_H
#define L1TRIGGER_TRACKFINDINGTRACKLET_L1STUBTRIPLETBUILDER_H

#include "DataFormats/L1TrackTrigger/interface/L1StubTriplet.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"

namespace trklet {

  L1StubTriplet makeL1StubTriplet(const Stub* inner,
                                 const Stub* middle,
                                 const Stub* outer,
                                 unsigned int sector,
                                 int region,
                                 int tpdunit);


  L1StubTriplet makeL1StubTriplet(const Stub* inner,
                                 const Stub* middle,
                                 const Stub* outer,
                                 unsigned int sector,
                                 int region,
                                 int tpdunit,
                                 int rzbinfirst_out,
                                 int rzbinfirst_in,
                                 int rzbinfirst_out_pair,
                                 int rzdiffmax_out,
                                 int rzdiffmax_in,
                                 int rzdiffmax_out_pair,
                                 int rzbin_out,
                                 int rzbin_in,
                                 int ibin_in,
                                 int ibin_out,
                                 int inner_rzbin,
                                 int outer_rzbin                
                                 );
}

#endif