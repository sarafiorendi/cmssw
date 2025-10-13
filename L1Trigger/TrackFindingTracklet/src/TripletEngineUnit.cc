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

  // see d10 is: in 6, mid 1, out 2

  bool applyPairCut = false; 
  bool lut_ok = false; // FIXME
  int rzbinfirst_out_new = 0;
  int rzdiffmax_out_new = 999;
//   // retrieve here the LUT for the triplet in case we are in seed 8
//   if ((trpdata_.stub_->layerdisk() == 2 && innervmstub.stub()->layerdisk() == 1) || 
//       (trpdata_.stub_->layerdisk() == 4 && innervmstub.stub()->layerdisk() == 3) ||
//       (trpdata_.stub_->layerdisk() == 1 && innervmstub.stub()->layerdisk() == 6) 
//      ) {

//     if (innervmstub.stub()->layerdisk() == 6) 
//       std::cout << "inside retrieving LUT for seed " 
//                 << innervmstub.stub()->layerdisk() << " " 
//                 << trpdata_.stub_->layerdisk() << " " 
//                 << outervmstub.stub()->layerdisk() << " " 
//                 << std::endl;   

    int middle_indexz = (((1 << (trpdata_.stub_->z().nbits() - 1)) + trpdata_.stub_->z().value()) >> (trpdata_.stub_->z().nbits() - settings_->vmrlutzbits(trpdata_.stub_->layerdisk()) ));
    int middle_indexr = (((1 << (trpdata_.stub_->r().nbits() - 1)) + trpdata_.stub_->r().value()) >> (trpdata_.stub_->r().nbits() - settings_->vmrlutrbits(trpdata_.stub_->layerdisk()) ));
    int inner_indexz = (((1 << (innervmstub.stub()->z().nbits() - 1)) + innervmstub.stub()->z().value()) >> (innervmstub.stub()->z().nbits() - settings_->vmrlutzbits(innervmstub.stub()->layerdisk()) ));
    int inner_indexr = (((1 << (innervmstub.stub()->r().nbits() - 1)) + innervmstub.stub()->r().value()) >> (innervmstub.stub()->r().nbits() - settings_->vmrlutrbits(innervmstub.stub()->layerdisk()) ));
    int shift_middle = settings_->vmrlutrbits(innervmstub.stub()->layerdisk()) + settings_->vmrlutzbits(innervmstub.stub()->layerdisk());
//     std::cout << ((middle_indexz << (settings_->vmrlutrbits(trpdata_.stub_->layerdisk()) + shift_middle)) + 
//                  (middle_indexr << shift_middle) + 
//                  (inner_indexz << settings_->vmrlutrbits(innervmstub.stub()->layerdisk())) + 
//                  inner_indexr )
//               << std::endl;
    int lutval_pair = trpdata_.outerpairtable_->lookup(
      (middle_indexz << (settings_->vmrlutrbits(trpdata_.stub_->layerdisk()) + shift_middle)) + 
      (middle_indexr << shift_middle) + 
      (inner_indexz << settings_->vmrlutrbits(innervmstub.stub()->layerdisk())) + 
      inner_indexr
    );
//       if (innervmstub.stub()->layerdisk() == 6) 
//         std::cout << "seed 10, lut = "  << lutval_pair
//                   << std::endl;   
//     

    if (lutval_pair != -1) {
      lut_ok = true;
//       if (innervmstub.stub()->layerdisk() == 6) 
//         std::cout << "seed 10, lut is ok" 
//                   << std::endl;   
      // retrieve cut value
      unsigned int lutwidth = settings_->lutwidthtabextended(0, 8); // always 21
      FPGAWord lookupbits(lutval_pair, lutwidth, true, __LINE__, __FILE__);
      int NFINERZBITS = 3;
      int nbitsrzbin_out = 3;    // N_RZBITS = 3; //number of bit for the r/z bins. it is 2 for seed 11
      rzbinfirst_out_new = lookupbits.bits(0, NFINERZBITS);   // NFINERZBITS = 3;   //number of bit for r or z within a r/z bin
      int next_out = lookupbits.bits(NFINERZBITS, 1);
      int start_out = lookupbits.bits(NFINERZBITS + 1, nbitsrzbin_out);  // first rz bin projection
      rzdiffmax_out_new = lookupbits.bits(NFINERZBITS + 1 + nbitsrzbin_out, NFINERZBITS);
      
      if (start_out !=  trpdata_.start_out_){
        applyPairCut = false;
        if (start_out == trpdata_.start_out_ + 1){
          rzbinfirst_out_new += (1 << NFINERZBITS);
          applyPairCut = true; // this was TRUE for pair cut
        }
      }  
//       bool print_csv_lut = true;
//       if (print_csv_lut){
//         std::cout << trpdata_.stub_->rapprox() << "," << trpdata_.stub_->zapprox() << "," 
//                   << innervmstub.stub()->rapprox() << "," << innervmstub.stub()->zapprox() << "," 
// 	            << rzbinfirst_out_new << "," << rzdiffmax_out_new << "," << start_out
// 	            << trpdata_.rzbinfirst_out_ << "," << trpdata_.rzdiffmax_out_   
// 	            << std::endl ;
//       }
    }
//   }// end retrieve LUT triplet

//   applyPairCut = false;

  // check if r/z of outer/inner stubs is within projection range
  int rzbin_out = (outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1));
  int rzbin_in = (innervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1));

  if (trpdata_.start_out_ != ibin_out)  // if we are in the the "next" bin, add 8 to the bin number 
    rzbin_out += (1 << NFINERZBITS);

  bool apply_outer_cut = true;  
  if (apply_outer_cut && ( rzbin_out < trpdata_.rzbinfirst_out_ || rzbin_out - trpdata_.rzbinfirst_out_ > trpdata_.rzdiffmax_out_)) {
    if (settings_->debugTracklet()) {
      edm::LogVerbatim("Tracklet") << "Outer stub rejected because of wrong r/z bin";
    }
  } else {  // condition on outer stub satisfied
    if ((trpdata_.start_in_ != ibin_in))
      rzbin_in += (1 << NFINERZBITS);

    bool apply_inner_cut = true;  
    if (apply_inner_cut && (rzbin_in < trpdata_.rzbinfirst_in_ || rzbin_in - trpdata_.rzbinfirst_in_ > trpdata_.rzdiffmax_in_)) {
      if (settings_->debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "Inner stub rejected because of wrong r/z bin";
      }
    } else {  // condition on both inner and outer stubs satisfied
      if (applyPairCut && 
//            ( (trpdata_.stub_->layerdisk() == 2 && innervmstub.stub()->layerdisk() == 1) || 
//              (trpdata_.stub_->layerdisk() == 4 && innervmstub.stub()->layerdisk() == 3) ||
//              (trpdata_.stub_->layerdisk() == 1 && innervmstub.stub()->layerdisk() == 6)
//              ) && 
           lut_ok){
//         if (innervmstub.stub()->layerdisk() == 6) std::cout << "cutting for seed 10" << std::endl;   
        if ( (rzbin_out < rzbinfirst_out_new || rzbin_out - rzbinfirst_out_new > rzdiffmax_out_new)){
          if (settings_->debugTracklet()) {
            edm::LogVerbatim("Tracklet") << "Outer stub rejected because of wrong r/z bin from pair";
          }
        }
        else{
            candtriplet_ =
                std::tuple<const Stub*, const Stub*, const Stub*>(innervmstub.stub(), trpdata_.stub_, outervmstub.stub());
            goodtriplet_ = true;
            
            FPGAWord phicorr = innervmstub.stub()->phicorr();
            int innerfinephi = phicorr.bits(phicorr.nbits() - 8, 8);
            
            std::cout << innervmstub.finephi().value() << " "
                      << (innervmstub.stub()->phicorr().value()) << " "
                      << (innerfinephi)
//                       << innervmstub.stub()->phicorr().bits()
//       int innerfinephi = phicorr.bits(phicorr.nbits() - nbitsfinephi_, nbitsfinephi_);
                      << std::endl;
                                              
            L1StubTriplet myTriplet;
            myTriplet.setStubRapprox(0, innervmstub.stub()->rapprox());
            myTriplet.setStubRapprox(1, trpdata_.stub_->rapprox());
            myTriplet.setStubRapprox(2, outervmstub.stub()->rapprox());
        
            myTriplet.setStubRValue(0, innervmstub.stub()->r().value());
            myTriplet.setStubRValue(1, trpdata_.stub_->r().value());
            myTriplet.setStubRValue(2, outervmstub.stub()->r().value());
          
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
            myTriplet.setFirstBinOutPair(rzbinfirst_out_new);
            myTriplet.setDiffMaxOut(trpdata_.rzdiffmax_out_);
            myTriplet.setDiffMaxIn(trpdata_.rzdiffmax_in_);
            myTriplet.setDiffMaxOutPair(rzdiffmax_out_new);
            myTriplet.setRZEffOut(rzbin_out);
            myTriplet.setRZEffIn(rzbin_in);
            myTriplet.setLargeBinIn(ibin_in);
            myTriplet.setLargeBinOutPair(ibin_out);
            
            foundtriplets.push_back(myTriplet);
        }  
      } // end apply pair cut and lut ok
      else{
            candtriplet_ =
                std::tuple<const Stub*, const Stub*, const Stub*>(innervmstub.stub(), trpdata_.stub_, outervmstub.stub());
            goodtriplet_ = true;
                                              
            L1StubTriplet myTriplet;
            myTriplet.setStubRapprox(0, innervmstub.stub()->rapprox());
            myTriplet.setStubRapprox(1, trpdata_.stub_->rapprox());
            myTriplet.setStubRapprox(2, outervmstub.stub()->rapprox());
        
            myTriplet.setStubRValue(0, innervmstub.stub()->r().value());
            myTriplet.setStubRValue(1, trpdata_.stub_->r().value());
            myTriplet.setStubRValue(2, outervmstub.stub()->r().value());
          
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
            myTriplet.setFirstBinOutPair(rzbinfirst_out_new);
            myTriplet.setDiffMaxOut(trpdata_.rzdiffmax_out_);
            myTriplet.setDiffMaxIn(trpdata_.rzdiffmax_in_);
            myTriplet.setDiffMaxOutPair(rzdiffmax_out_new);
            myTriplet.setRZEffOut(rzbin_out);
            myTriplet.setRZEffIn(rzbin_in);
            myTriplet.setLargeBinIn(ibin_in);
            myTriplet.setLargeBinOutPair(ibin_out);
            
            foundtriplets.push_back(myTriplet);
      }  
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
