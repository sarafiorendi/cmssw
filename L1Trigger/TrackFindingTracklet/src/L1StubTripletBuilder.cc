#include "L1Trigger/TrackFindingTracklet/interface/L1StubTripletBuilder.h"

namespace trklet {

  L1StubTriplet makeL1StubTriplet(const Stub* inner,
                                 const Stub* middle,
                                 const Stub* outer,
                                 unsigned int sector,
                                 int region,
                                 int tpdunit) {
  
    L1StubTriplet triplet;
  
    triplet.setStubRapprox(0, inner->rapprox());
    triplet.setStubRapprox(1, middle->rapprox());
    triplet.setStubRapprox(2, outer->rapprox());
  
    triplet.setStubZapprox(0, inner->zapprox());
    triplet.setStubZapprox(1, middle->zapprox());
    triplet.setStubZapprox(2, outer->zapprox());
  
    triplet.setStubBend(0, inner->bend().value());
    triplet.setStubBend(1, middle->bend().value());
    triplet.setStubBend(2, outer->bend().value());
  
    triplet.setStubPhi(0, inner->phiapprox(0., 0));
    triplet.setStubPhi(1, middle->phiapprox(0., 0));
    triplet.setStubPhi(2, outer->phiapprox(0., 0));
  
    triplet.setStubIndex(0, inner->stubindex().value());
    triplet.setStubIndex(1, middle->stubindex().value());
    triplet.setStubIndex(2, outer->stubindex().value());
  
    triplet.setStubLayerdisk(0, inner->layerdisk());
    triplet.setStubLayerdisk(1, middle->layerdisk());
    triplet.setStubLayerdisk(2, outer->layerdisk());
  
    triplet.setSector(sector);
    triplet.setRegion(region);
    triplet.setTPDUnit(tpdunit);
  
    return triplet;
  }


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
                                 int ibin_out,
                                 int ibin_in,
                                 int inner_rzbin,
                                 int outer_rzbin
                                 ) {  

    L1StubTriplet triplet;
  
    triplet.setStubRapprox(0, inner->rapprox());
    triplet.setStubRapprox(1, middle->rapprox());
    triplet.setStubRapprox(2, outer->rapprox());
  
    triplet.setStubZapprox(0, inner->zapprox());
    triplet.setStubZapprox(1, middle->zapprox());
    triplet.setStubZapprox(2, outer->zapprox());
  
    triplet.setStubBend(0, inner->bend().value());
    triplet.setStubBend(1, middle->bend().value());
    triplet.setStubBend(2, outer->bend().value());
  
    triplet.setStubPhi(0, inner->phiapprox(0., 0));
    triplet.setStubPhi(1, middle->phiapprox(0., 0));
    triplet.setStubPhi(2, outer->phiapprox(0., 0));
  
    triplet.setStubIndex(0, inner->stubindex().value());
    triplet.setStubIndex(1, middle->stubindex().value());
    triplet.setStubIndex(2, outer->stubindex().value());
  
    triplet.setStubLayerdisk(0, inner->layerdisk());
    triplet.setStubLayerdisk(1, middle->layerdisk());
    triplet.setStubLayerdisk(2, outer->layerdisk());

    triplet.setStubRValue(0, inner->r().value());
    triplet.setStubRValue(1, middle->r().value());
    triplet.setStubRValue(2, outer->r().value());

    triplet.setStubRZbin(0,  inner_rzbin); //(innervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1)));
    triplet.setStubRZbin(1, 0); // dummy fill
    triplet.setStubRZbin(2, outer_rzbin);//(outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1)));
  
    triplet.setSector(sector);
    triplet.setRegion(region);
    triplet.setTPDUnit(tpdunit);

    triplet.setFirstBinOut(rzbinfirst_out);
    triplet.setFirstBinIn(rzbinfirst_in);
    triplet.setFirstBinOutPair(rzbinfirst_out_pair);
    triplet.setDiffMaxOut(rzdiffmax_out);
    triplet.setDiffMaxIn(rzdiffmax_in);
    triplet.setDiffMaxOutPair(rzdiffmax_out_pair);
    triplet.setRZEffOut(rzbin_out);
    triplet.setRZEffIn(rzbin_in);
    triplet.setLargeBinIn(ibin_in);
    triplet.setLargeBinOutPair(ibin_out);
            

    return triplet;
  }

}
