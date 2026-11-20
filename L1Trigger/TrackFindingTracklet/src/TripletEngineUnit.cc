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
                                     unsigned int iAllStub,
                                     unsigned int nbitsfinephi,
                                     unsigned int nbitsfinephiouterdiff,
                                     unsigned int nbitsfinephiinnerdiff,
                                     const TrackletLUT* pttablemiddle,
                                     const TrackletLUT* pttableouter,
                                     const TrackletLUT* pttablemiddlein,
                                     const TrackletLUT* pttableinner,
                                     std::vector<VMStubsTEMemory*> innervmstubs,
                                     std::vector<VMStubsTEMemory*> outervmstubs)
    : innervmstubs_(std::move(innervmstubs)),
      outervmstubs_(std::move(outervmstubs)),
      settings_(settings),
      layerdisk1_(layerdisk1),
      layerdisk2_(layerdisk2),
      layerdisk3_(layerdisk3),
      iSeed_(iSeed),
      iAllStub_(iAllStub),
      nbitsfinephi_(nbitsfinephi),
      nbitsfinephiouterdiff_(nbitsfinephiouterdiff),
      nbitsfinephiinnerdiff_(nbitsfinephiinnerdiff),
      nearfull_(false),
      nmem_out_(0),
      nmem_in_(0),
      istub_out_(0),
      istub_in_(0),
      next_out_(0),
      next_in_(0),
      nstub_out_(0),
      nstub_in_(0),
      outmem_(0),
      inmem_(0),
      nproj_out_(0),
      nproj_in_(0),
      idle_(true),
      pttablemiddle_(pttablemiddle),
      pttableouter_(pttableouter),
      pttablemiddlein_(pttablemiddlein),
      pttableinner_(pttableinner),    
      candtriplets_(3) {}

void TripletEngineUnit::init(const TrpEData& trpdata) {
  trpdata_ = trpdata;
  istub_out_ = 0;
  istub_in_ = 0;
  nproj_out_ = 0;
  nproj_in_ = 0;
  idle_ = false;

  assert(!trpdata_.projbin_out_.empty() && !trpdata_.projbin_in_.empty());
  std::tie(next_out_, outmem_, nstub_out_, phi_out_) = trpdata_.projbin_out_[0];
  std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[0];
}

void TripletEngineUnit::reset() {
  idle_ = true;
  goodtriplet_ = false;
  goodtriplet__ = false;
  candtriplets_.reset();
}

void TripletEngineUnit::step() {
  if (goodtriplet__) {
    candtriplets_.store(candtriplet__);
  }

  goodtriplet__ = goodtriplet_;
  candtriplet__ = candtriplet_;

  goodtriplet_ = false;

  if (idle_ || nearfull_) {
    return;
  }

  if (!nextGoodInner()) {      // nothing usable left for this outer stub
    advanceOuterStub();
    if (idle_) return;
    if (!nextGoodInner()) {    // every inner stub is bad: no outer stub can make a triplet
      idle_ = true;            // this middle stub is done
      return;
    }
  }
  
  // get inner and outer projected stub for certain next value
  int ibin_out = trpdata_.start_out_ + next_out_; // this is the coarse bin
  const VMStubTE& outervmstub = outervmstubs_[outmem_]->getVMStubTEBinned(ibin_out, istub_out_);

  int rzbin_out = (outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1));
  if (trpdata_.start_out_ != ibin_out)  // if looking at the "next" bin
    rzbin_out += (1 << NFINERZBITS);

  bool rz_out_ok = !(rzbin_out < trpdata_.rzbinfirst_out_ ||
                     rzbin_out - trpdata_.rzbinfirst_out_ > trpdata_.rzdiffmax_out_);

  // phi / bend cuts on outer
  FPGAWord ifinephiouter = outervmstub.finephi();
  assert(ifinephiouter == outervmstub.finephi());
  // retrieve lut value: first calculate the phi value, as done in TrackletLUT
  int outerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(outervmstub.stub()->layerdisk()))) +
                     phi_out_ * (1 << settings_->nfinephi(1, iSeed_)) + ifinephiouter.value();
                     
  if (phi_out_ >= 2 * settings_->nvmte(1, iSeed_))                     
      outerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(outervmstub.stub()->layerdisk()))) +
                     (phi_out_ - 24) * (1 << settings_->nfinephi(1, iSeed_)) + ifinephiouter.value(); // 24 = 3 * settings_.nvmte(1, iSeed) 

  int idphi_out = outerfinephi - trpdata_.middlefinephi_;
  bool inrange_out = true;
  if (iSeed_ == 8 || iSeed_ == 9) {
    inrange_out = (idphi_out < (1 << (nbitsfinephiouterdiff_ - 1))) && (idphi_out >= -(1 << (nbitsfinephiouterdiff_ - 1)));
  }  
  int idphi_out_for_index = idphi_out & ((1 << nbitsfinephiouterdiff_) - 1);

  bool pass_pt_cut_out = true;
  FPGAWord outerbend = outervmstub.bend();
  if (iSeed_ == 8 || iSeed_ == 9) {
    int ptouterindex = (idphi_out_for_index << outerbend.nbits()) + outerbend.value();
    int ptmiddleindex = (idphi_out_for_index << trpdata_.middlebend_.nbits()) + trpdata_.middlebend_.value();
    pass_pt_cut_out = pttablemiddle_->lookup(ptmiddleindex) && pttableouter_->lookup(ptouterindex);
  }

  if (!(rz_out_ok && inrange_out && pass_pt_cut_out)) {
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "Outer stub rejected";
    }
    advanceOuterStub();   
    return;  
  }

  // retrieve inner stub
  int ibin_in = trpdata_.start_in_ + next_in_;
  const VMStubTE& innervmstub = innervmstubs_[inmem_]->getVMStubTEBinned(ibin_in, istub_in_);

  // info for phi / bend cuts on inner
  FPGAWord ifinephiinner = innervmstub.finephi();
  assert(ifinephiinner == innervmstub.finephi());
  int innerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(innervmstub.stub()->layerdisk()))) +
                     phi_in_ * (1 << settings_->nfinephi(1, iSeed_)) + ifinephiinner.value();
  if (phi_in_ >= 2 * settings_->nvmte(1, iSeed_))                     
      innerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(innervmstub.stub()->layerdisk()))) +
                     (phi_in_ - 24) * (1 << settings_->nfinephi(1, iSeed_)) + ifinephiinner.value(); // 24 = 3 * settings_.nvmte(1, iSeed) 

  int idphi_in = innerfinephi - trpdata_.middlefinephi_;
  bool inrange_in = true;
  if (iSeed_ == 8 || iSeed_ == 9) {
    inrange_in = (idphi_in < (1 << (nbitsfinephiinnerdiff_ - 1))) && (idphi_in >= -(1 << (nbitsfinephiinnerdiff_ - 1)));
  }  
  int idphi_in_for_index = idphi_in & ((1 << nbitsfinephiinnerdiff_) - 1);

  // check if r/z of inner stubs is within projection range
  int rzbin_in = (innervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1));
  if ((trpdata_.start_in_ != ibin_in))
    rzbin_in += (1 << NFINERZBITS);

  bool rz_in_ok = !(rzbin_in < trpdata_.rzbinfirst_in_ ||
                    rzbin_in - trpdata_.rzbinfirst_in_ > trpdata_.rzdiffmax_in_);

  // phi / bend cuts on inner
  bool pass_pt_cut_inner = true;
  FPGAWord innerbend = innervmstub.bend();
  if (iSeed_ == 8 || iSeed_ == 9) {
    int ptinnerindex = (idphi_in_for_index << innerbend.nbits()) + innerbend.value();
    int ptmiddleinindex = (idphi_in_for_index << trpdata_.middlebend_.nbits()) + trpdata_.middlebend_.value();
    pass_pt_cut_inner = pttablemiddlein_->lookup(ptmiddleinindex) && pttableinner_->lookup(ptinnerindex);
  }

  if (!(rz_in_ok && inrange_in && pass_pt_cut_inner)) {
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "Inner stub rejected";
    }
    trpdata_.bad_inner_stub_[(inmem_ * trpdata_.max_nbins_ + trpdata_.start_in_ + next_in_) * trpdata_.max_nstubs_ + istub_in_] = true;
  } else {  // condition on inner stub satisfied

    candtriplet_ =
        std::tuple<const Stub*, const Stub*, const Stub*>(innervmstub.stub(), trpdata_.stub_, outervmstub.stub());
    goodtriplet_ = true;
  }

  // go to next projection (looping through all inner stubs for each outer stub)
  advanceIndices();
  if (idle_) return;
}


void TripletEngineUnit::advanceIndices() {
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
        if (nproj_out_ >= trpdata_.projbin_out_.size()) {  // if gone through all out proj bins, reset everything and stop engine unit
          istub_in_ = 0;
          istub_out_ = 0;
          nproj_in_ = 0;
          nproj_out_ = 0;
//           if (iSeed_ == 8) std::cout << "METRIC_TEU_EXACT_TRIPLETS: " << n_triplets_tested_ << std::endl;
          idle_ = true;
          return;
        }
        // get next out proj bin
        std::tie(next_out_, outmem_, nstub_out_, phi_out_) = trpdata_.projbin_out_[nproj_out_];
      }
    }
    // get next in proj bin
    std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[nproj_in_];
  }
}

void TripletEngineUnit::advanceInnerBin() {

  nproj_in_++;
  istub_in_ = 0;
  if (nproj_in_ >= trpdata_.projbin_in_.size()) {  // if gone through all in proj bins, move to next out stub
    istub_out_++;
    nproj_in_ = 0;
    if (istub_out_ >= nstub_out_) {  // if gone through all out stubs, move to next out proj bin
      nproj_out_++;
      istub_out_ = 0;
      if (nproj_out_ >= trpdata_.projbin_out_.size()) {  // if gone through all out proj bins, reset everything and stop engine unit
        istub_in_ = 0;
        istub_out_ = 0;
        nproj_in_ = 0;
        nproj_out_ = 0;
//         if (iSeed_ == 8) std::cout << "METRIC_TEU_EXACT_TRIPLETS: " << n_triplets_tested_ << std::endl;
        idle_ = true;
        return;
      }
      // get next out proj bin
      std::tie(next_out_, outmem_, nstub_out_, phi_out_) = trpdata_.projbin_out_[nproj_out_];
    }
  }
  // get next in proj bin
  std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[nproj_in_];
}

void TripletEngineUnit::advanceOuterStub() {
  istub_in_ = 0;
  nproj_in_ = 0;
  istub_out_++;
  if (istub_out_ >= nstub_out_) {
    istub_out_ = 0;
    nproj_out_++;
    if (nproj_out_ >= trpdata_.projbin_out_.size()) {
      nproj_out_ = 0;
      idle_ = true;
      return;
    }
    std::tie(next_out_, outmem_, nstub_out_, phi_out_) = trpdata_.projbin_out_[nproj_out_];
  }
  std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[0];
}


// index of the current inner stub in bad_inner_stub_
size_t TripletEngineUnit::badInnerIndex() const {
  return (inmem_ * trpdata_.max_nbins_ + trpdata_.start_in_ + next_in_) * trpdata_.max_nstubs_ + istub_in_;
}

// From the current inner position, move forward (within this outer stub's inner loop) to the first
// inner stub NOT marked bad. Returns false if there is none left.
// Hardware: a find-first-set-bit over the bad-stub mask, one step.
bool TripletEngineUnit::nextGoodInner() {
  while (nproj_in_ < trpdata_.projbin_in_.size()) {
    std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[nproj_in_];
    while (istub_in_ < nstub_in_) {
      if (!trpdata_.bad_inner_stub_[badInnerIndex()])
        return true;
      istub_in_++;
    }
    nproj_in_++;
    istub_in_ = 0;
  }
  return false;
}
