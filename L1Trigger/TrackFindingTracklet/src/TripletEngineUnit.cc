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
  int ibin_in = trpdata_.start_in_ + next_in_; // big vm bin to get the stub memory
  const VMStubTE& outervmstub = outervmstubs_[outmem_]->getVMStubTEBinned(ibin_out, istub_out_);
  const VMStubTE& innervmstub = innervmstubs_[inmem_]->getVMStubTEBinned(ibin_in, istub_in_);

  // check if r/z of outer stub is within projection range
  int rzbin_out = (outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1));
  int rzbin_in  = (innervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1)); // 3 bits 100//sara
  
//   if (!innervmstub.stub()->isPSmodule())
//     std::cout << "inner one is 2S " << std::endl;
//   if (!outervmstub.stub()->isPSmodule())
//     std::cout << "outer one is 2S " << std::endl;

  if (trpdata_.start_out_ != ibin_out)
    rzbin_out += 8;
    
  if (rzbin_out < trpdata_.rzbinfirst_out_ || rzbin_out - trpdata_.rzbinfirst_out_ > trpdata_.rzdiffmax_out_) {
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "Outer stub rejected because of wrong r/z bin";
    }
  } else {
    if ((trpdata_.start_in_ != ibin_in))
      rzbin_in += (1 << NFINERZBITS);
    //  do not passcondition on inner
    bool applyCut = true;
    if ( applyCut && (rzbin_in < trpdata_.rzbinfirst_in_ || rzbin_in - trpdata_.rzbinfirst_in_ > trpdata_.rzdiffmax_in_) ){
      if (settings_->debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "Inner stub rejected because of wrong r/z bin";
      }
    } else {  // here it should have passed both inner and outer stub conditions
      candtriplet_ =
          std::tuple<const Stub*, const Stub*, const Stub*>(innervmstub.stub(), trpdata_.stub_, outervmstub.stub());
      goodtriplet_ = true;
                                        
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
      myTriplet.setStubRZbin(1, 0); // dummy fill
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
      
      myTriplet.setFirstBinOut(trpdata_.rzbinfirst_out_);
      myTriplet.setFirstBinIn(trpdata_.rzbinfirst_in_);
      myTriplet.setDiffMaxOut(trpdata_.rzdiffmax_out_);
      myTriplet.setDiffMaxIn(trpdata_.rzdiffmax_in_);
      myTriplet.setRZEffOut(rzbin_out);
      myTriplet.setRZEffIn(rzbin_in);
      
      
      foundtriplets.push_back(myTriplet);
    }
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
