#include "L1Trigger/TrackFindingTracklet/interface/TripletEngineUnit.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"

using namespace trklet;

// TripletEngineUnit
//
// This script sets a processing unit for the TrackletProcessorDisplaced
// based on information from the middle stub and its projections in and out,
// and a new triplet seed is checked and created for each processing step
//
// Author: Claire Savard, Nov. 2024

TripletEngineUnit::TripletEngineUnit(const Settings* const settings,
                                     unsigned int layerdisk1,
                                     unsigned int layerdisk2,
                                     unsigned int layerdisk3,
                                     unsigned int iSeed,
                                     std::vector<VMStubsTEMemory*> innervmstubs,
                                     std::vector<VMStubsTEMemory*> outervmstubs)
    : settings_(settings), candtriplets_(3) {
  idle_ = true;
  layerdisk1_ = layerdisk1;
  layerdisk2_ = layerdisk2;
  layerdisk3_ = layerdisk3;
  iSeed_ = iSeed;
  innervmstubs_ = innervmstubs;
  outervmstubs_ = outervmstubs;
}

void TripletEngineUnit::init(const TrpEData& trpdata) {
  trpdata_ = trpdata;
  istub_out_ = 0;
  istub_in_ = 0;
  nproj_out_ = 0;
  nproj_in_ = 0;
  idle_ = false;

  assert(!trpdata_.projbin_out_.empty() && !trpdata_.projbin_in_.empty());
  std::tie(next_out_, outmem_, nstub_out_) = trpdata_.projbin_out_[0];
  std::tie(next_in_, inmem_, nstub_in_) = trpdata_.projbin_in_[0];
}

void TripletEngineUnit::reset() {
  idle_ = true;
  goodtriplet_ = false;
  goodtriplet__ = false;
  candtriplets_.reset();
}

// void TripletEngineUnit::step() {
void TripletEngineUnit::step(std::vector<L1StubTriplet>& foundtriplets, unsigned int iSector, int iTC, int count_trpunits) {
  if (goodtriplet__) {
    candtriplets_.store(candtriplet__);
  }

  goodtriplet__ = goodtriplet_;
  candtriplet__ = candtriplet_;

  goodtriplet_ = false;

  if (idle_ || nearfull_) {
    return;
  }

  // get inner and outer projected stub for certain next value
  int ibin_out = trpdata_.start_out_ + next_out_;
  int ibin_in = trpdata_.start_in_ + next_in_;
  const VMStubTE& outervmstub = outervmstubs_[outmem_]->getVMStubTEBinned(ibin_out, istub_out_);
  const VMStubTE& innervmstub = innervmstubs_[inmem_]->getVMStubTEBinned(ibin_in, istub_in_);

  // check if r/z of outer stub is within projection range
  int rzbin = (outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1));
//   std::cout << "rzbin: " << rzbin << std::endl;

  if (trpdata_.start_out_ != ibin_out)
    rzbin += 8;
  if (rzbin < trpdata_.rzbinfirst_out_ || rzbin - trpdata_.rzbinfirst_out_ > trpdata_.rzdiffmax_out_) {
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "Outer stub rejected because of wrong r/z bin";
    }
  } else {
    candtriplet_ =
        std::tuple<const Stub*, const Stub*, const Stub*>(innervmstub.stub(), trpdata_.stub_, outervmstub.stub());
    goodtriplet_ = true;
//     std::cout << "goodtriplet: \n\t r:" << innervmstub.stub()->rapprox() << " \t z:" << innervmstub.stub()->zapprox() << " \tbend" << innervmstub.stub()->bend().value() << " \t" << innervmstub.stub()->stubindex().value() 
//                                         << " -> layerdisk: " << innervmstub.stub()-> layerdisk() // << " \t" << innervmstub.stub()->phiapprox() 
//                                         <<" \n\t r:" << trpdata_.stub_->rapprox() << " \t z:" << trpdata_.stub_->zapprox()  << " \tbend" << trpdata_.stub_->bend().value() << " \t" << trpdata_.stub_->stubindex().value()// << " \t" << trpdata_.stub_->phiapprox() 
//                                         <<" \n\t r:" << outervmstub.stub()->rapprox()  << " \t z:" << outervmstub.stub()->zapprox()  << " \tbend" << outervmstub.stub()->bend().value() << " \t" << outervmstub.stub()->stubindex().value()
//                                         << std::endl;
                                        
    L1StubTriplet myTriplet;
    myTriplet.setStubRapprox(0, innervmstub.stub()->rapprox());
    myTriplet.setStubRapprox(1, trpdata_.stub_->rapprox());
    myTriplet.setStubRapprox(2, outervmstub.stub()->rapprox());

    myTriplet.setStubZapprox(0, innervmstub.stub()->zapprox());
    myTriplet.setStubZapprox(1, trpdata_.stub_->zapprox());
    myTriplet.setStubZapprox(2, outervmstub.stub()->zapprox());

    myTriplet.setStubBend(0, innervmstub.stub()->bend().value());
    myTriplet.setStubBend(1, trpdata_.stub_->bend().value());
    myTriplet.setStubBend(2, outervmstub.stub()->bend().value());

    myTriplet.setStubRZbin(0, (innervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1)));
    myTriplet.setStubRZbin(1, (trpdata_.rzbinfirst_out_));
    myTriplet.setStubRZbin(2, (outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1)));
    
    myTriplet.setStubIndex(0, innervmstub.stub()->stubindex().value());
    myTriplet.setStubIndex(1, trpdata_.stub_->stubindex().value());
    myTriplet.setStubIndex(2, outervmstub.stub()->stubindex().value());

    myTriplet.setStubLayerdisk(0, innervmstub.stub()->layerdisk());
    myTriplet.setStubLayerdisk(1, trpdata_.stub_->layerdisk());
    myTriplet.setStubLayerdisk(2, outervmstub.stub()->layerdisk());
    
    myTriplet.setSector(iSector);
    myTriplet.setRegion(iTC);
    myTriplet.setTPDUnit(count_trpunits);
    foundtriplets.push_back(myTriplet);
  }

  // go to next projection (looping through all inner stubs for each outer stub)
  istub_in_++;
  if (istub_in_ >= nstub_in_) {  // if gone through all in stubs, move to next in proj bin
    nproj_in_++;
    istub_in_ = 0;
    if (nproj_in_ >= trpdata_.projbin_in_.size()) {  // if gone through all in proj bins, move to next out stub
      istub_out_++;
      nproj_in_ = 0;
      if (istub_out_ >= nstub_out_) {  // if gone through all out stubs, move to next out proj bin
        nproj_out_++;
        istub_out_ = 0;
        if (nproj_out_ >=
            trpdata_.projbin_out_.size()) {  // if gone through all out proj bins, reset everything and stop engine unit
          istub_in_ = 0;
          istub_out_ = 0;
          nproj_in_ = 0;
          nproj_out_ = 0;
          idle_ = true;
          return;
        }
        // get next out proj bin
        std::tie(next_out_, outmem_, nstub_out_) = trpdata_.projbin_out_[nproj_out_];
      }
    }
    // get next in proj bin
    std::tie(next_in_, inmem_, nstub_in_) = trpdata_.projbin_in_[nproj_in_];
  }
}
