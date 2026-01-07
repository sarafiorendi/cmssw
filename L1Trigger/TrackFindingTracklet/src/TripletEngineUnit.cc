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
                                     const TrackletLUT* pttablemiddlenew,
                                     const TrackletLUT* pttableouternew,
                                     std::vector<VMStubsTEMemory*> innervmstubs,
                                     std::vector<VMStubsTEMemory*> outervmstubs)
    : settings_(settings), 
      pttablemiddlenew_(pttablemiddlenew),
      pttableouternew_(pttableouternew),
      candtriplets_(3) {
  idle_ = true;
  layerdisk1_ = layerdisk1;
  layerdisk2_ = layerdisk2;
  layerdisk3_ = layerdisk3;
  iSeed_ = iSeed;
  iAllStub_ = iAllStub;
  nbitsfinephi_ = nbitsfinephi;
  nbitsfinephiouterdiff_ = nbitsfinephiouterdiff;
  nbitsfinephiinnerdiff_ = nbitsfinephiinnerdiff;
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
  std::tie(next_out_, outmem_, nstub_out_, phi_out_) = trpdata_.projbin_out_[0];
  std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[0];
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

//   std::cout << "\t\t [TEU] ibin_out/istub_out_ " << ibin_out << "/" << istub_out_ 
//             << "           ibin_in/istub_in_ " << ibin_in << "/" << istub_in_ 
//             << std::endl;

  // see d10 is: in 6, mid 1, out 2
  
  // set up needed info for phi / bend cuts on outer (for now)
  FPGAWord ifinephiouter = outervmstub.finephi();
  assert(ifinephiouter == outervmstub.finephi());
  // retrieve lut value 
  // first calculate the phi value, as done in TrackletLUT
  int outerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(outervmstub.stub()->layerdisk()))) +
                     phi_out_ * (1 << settings_->nfinephi(1, iSeed_)) + ifinephiouter.value();
                     
  if (phi_out_ >= 2 * settings_->nvmte(1, iSeed_))                     
      outerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(outervmstub.stub()->layerdisk()))) +
                     (phi_out_ - 24) * (1 << settings_->nfinephi(1, iSeed_)) + ifinephiouter.value(); // 24 = 3 * settings_.nvmte(1, iSeed) 

  // ireg corresponds to out_phi_region in the TPD
  // should be in the range from 0 to 7 and represent the 8 memories per each of the 4 large phi region
  // how to get it from here? in TPD it was
  // unsigned int out_phi_region = 
  //      (outervmstubs_[outmem]->phibin() - 1) - (outervmstubs_[outmem]->getName()[11] - 'A') * 8;
  int idphi_out = outerfinephi - trpdata_.middlefinephi_;
  bool inrange = (idphi_out < (1 << (nbitsfinephiouterdiff_ - 1))) && (idphi_out >= -(1 << (nbitsfinephiouterdiff_ - 1)));
  int idphi_out_for_index = idphi_out & ((1 << nbitsfinephiouterdiff_) - 1);

//   idphi_out = idphi_out & ((1 << nbitsfinephiouterdiff_) - 1);
//   std::cout << " idphi_out " << idphi_out << "   inrange = " << inrange  << std::endl;


  FPGAWord ifinephiinner = innervmstub.finephi();
  assert(ifinephiinner == innervmstub.finephi());
  int innerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(innervmstub.stub()->layerdisk()))) +
                     phi_in_ * (1 << settings_->nfinephi(1, iSeed_)) + ifinephiinner.value();
  if (phi_in_ >= 2 * settings_->nvmte(1, iSeed_))                     
      innerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(innervmstub.stub()->layerdisk()))) +
                     (phi_in_ - 24) * (1 << settings_->nfinephi(1, iSeed_)) + ifinephiinner.value(); // 24 = 3 * settings_.nvmte(1, iSeed) 

  int idphi_in = innerfinephi - trpdata_.middlefinephi_;
  bool inrange_in = (idphi_in < (1 << (nbitsfinephiinnerdiff_ - 1))) && (idphi_in >= -(1 << (nbitsfinephiinnerdiff_ - 1)));

  bool opposite_sign_dphi = (idphi_out * idphi_in > 0) ? 0 : 1;
  if (abs(idphi_out) < 10 && abs(idphi_in) < 10) opposite_sign_dphi = true;
//   idphi_in = idphi_in & ((1 << nbitsfinephiinnerdiff_) - 1);
//   if (iSeed_ == 8){
// //     std::cout << "   iAllStub_ " << iAllStub_ 
// //               << "   phi_out_ = " << phi_out_  
// //               << "   ifinephiouter.value() " << ifinephiouter.value()
// //               << std::endl;
// //     std::cout << " idphi_out " << idphi_out << " (cut is " << ( 1 << (nbitsfinephiouterdiff_ - 1)) << ")     inrange = " << inrange  << std::endl;
// //     std::cout << "   innerfinephi " << innerfinephi 
// //               << "   idphi_in " << idphi_in 
// //               << "   inrange_in = " << inrange_in  
// //               << "   nbitsfinephiinnerdiff_ = " << nbitsfinephiinnerdiff_  
// //               << std::endl;
//     std::cout << "   innerphi/innerfinephi " << innervmstub.stub()->phiapprox(0., 0) 
//               << "    " << innerfinephi 
// //               << "   middlephi " << trpdata_.stub_->phiapprox(0., 0) 
// //               << "   delta =  " << innervmstub.stub()->phiapprox(0., 0)  - trpdata_.stub_->phiapprox(0., 0) 
//               << "\n" << std::endl;
// //     std::cout << "   outerphi " << outervmstub.stub()->phiapprox(0., 0) 
// //               << "   middlephi " << trpdata_.stub_->phiapprox(0., 0) 
// //               << "   delta =  " << outervmstub.stub()->phiapprox(0., 0)  - trpdata_.stub_->phiapprox(0., 0) 
// //               << "   \t\t inrange = " << inrange 
// //               << "\n" << std::endl;
//   }
// 


  bool applyPairCut = true; 
  bool lut_ok = false; // FIXME
  int rzbinfirst_out_new = 0;
  int rzdiffmax_out_new = 999;

  int middle_indexz = (((1 << (trpdata_.stub_->z().nbits() - 1)) + trpdata_.stub_->z().value()) >> (trpdata_.stub_->z().nbits() - settings_->vmrlutzbits(trpdata_.stub_->layerdisk()) ));
  int middle_indexr = (((1 << (trpdata_.stub_->r().nbits() - 1)) + trpdata_.stub_->r().value()) >> (trpdata_.stub_->r().nbits() - settings_->vmrlutrbits(trpdata_.stub_->layerdisk()) ));
  int inner_indexz = (((1 << (innervmstub.stub()->z().nbits() - 1)) + innervmstub.stub()->z().value()) >> (innervmstub.stub()->z().nbits() - settings_->vmrlutzbits(innervmstub.stub()->layerdisk()) ));
  int inner_indexr = (((1 << (innervmstub.stub()->r().nbits() - 1)) + innervmstub.stub()->r().value()) >> (innervmstub.stub()->r().nbits() - settings_->vmrlutrbits(innervmstub.stub()->layerdisk()) ));
  int shift_middle = settings_->vmrlutrbits(innervmstub.stub()->layerdisk()) + settings_->vmrlutzbits(innervmstub.stub()->layerdisk());

  int lutval_pair = trpdata_.outerpairtable_->lookup(
    (middle_indexz << (settings_->vmrlutrbits(trpdata_.stub_->layerdisk()) + shift_middle)) + 
    (middle_indexr << shift_middle) + 
    (inner_indexz << settings_->vmrlutrbits(innervmstub.stub()->layerdisk())) + 
    inner_indexr
  );

  if (lutval_pair != -1) {
      lut_ok = true;

      // retrieve cut value
      unsigned int lutwidth = settings_->lutwidthtabextended(0, 8); // always 21
      FPGAWord lookupbits(lutval_pair, lutwidth, true, __LINE__, __FILE__);
      int NFINERZBITS = 3;
      int nbitsrzbin_out = 3;    // N_RZBITS = 3; //number of bit for the r/z bins. it is 2 for seed 11
      rzbinfirst_out_new = lookupbits.bits(0, NFINERZBITS);   // NFINERZBITS = 3;   //number of bit for r or z within a r/z bin
//       int next_out = lookupbits.bits(NFINERZBITS, 1);
      int start_out = lookupbits.bits(NFINERZBITS + 1, nbitsrzbin_out);  // first rz bin projection
      rzdiffmax_out_new = lookupbits.bits(NFINERZBITS + 1 + nbitsrzbin_out, NFINERZBITS);
      
      if (start_out !=  trpdata_.start_out_){
        applyPairCut = false;
        if (start_out == trpdata_.start_out_ + 1){
          rzbinfirst_out_new += (1 << NFINERZBITS);
          applyPairCut = true; // this was TRUE for pair cut
        }
      }  
//       applyPairCut = true;
//       bool print_csv_lut = true;
//       if (print_csv_lut){
//         std::cout << trpdata_.stub_->rapprox() << "," << trpdata_.stub_->zapprox() << "," 
//                   << innervmstub.stub()->rapprox() << "," << innervmstub.stub()->zapprox() << "," 
// 	            << rzbinfirst_out_new << "," << rzdiffmax_out_new << "," << start_out
// 	            << trpdata_.rzbinfirst_out_ << "," << trpdata_.rzdiffmax_out_   
// 	            << std::endl ;
//       }
  }

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
    
      
      FPGAWord outerbend = outervmstub.bend();
      FPGAWord innerbend = innervmstub.bend();
      bool pass_pt_cut = false;
      if (iSeed_ == 8){
        int ptouterindex = (idphi_out_for_index << outerbend.nbits()) + outerbend.value();
        int ptmiddleindex = (idphi_out_for_index << trpdata_.middlebend_.nbits()) + trpdata_.middlebend_.value();
//         std::cout << "ptouterindex = " << ptouterindex << std::endl;
//         std::cout << std::bitset<16>(ptouterindex) << std::endl;
//         std::cout << std::bitset<16>(ptmiddleindex) << std::endl;
//         std::cout << std::endl;
//         pass_pt_cut = pttableouternew_->lookup(ptouterindex);
        pass_pt_cut = pttablemiddlenew_->lookup(ptmiddleindex) && pttableouternew_->lookup(ptouterindex);
      }
      if (pass_pt_cut){
//         std::cout << "passed pt cut! " << std::endl;
      }
      // now apply cut on pT
//       if (iSeed_ == 8 && ( !inrange) ) {
//       if (iSeed_ == 8 && ( !inrange_in) ) {
//       if (iSeed_ == 8 && ( !(inrange && inrange_in) )) {
      if (iSeed_ == 8 && ( !(inrange && inrange_in && opposite_sign_dphi) )) {
//       if (iSeed_ == 8 && ( !(inrange && inrange_in && opposite_sign_dphi && pass_pt_cut) )) {
//       opposite_sign_dphi = true; // fake to emulate no cut for processing est 13:21
//       if (iSeed_ == 8 && ( !(opposite_sign_dphi) )) {
//       if (iSeed_ == 8 && ( ! (inrange && pass_pt_cut)) ) {
//       if (!(inrange && pttablemiddlenew_->lookup(ptmiddleindex) && pttableouternew_->lookup(ptouterindex))) {
        if (settings_->debugTracklet()) {
          edm::LogVerbatim("Tracklet") << " Stub pair rejected because of stub pt cut bends : "
                                       << settings_->benddecode(
                                              trpdata_.middlebend_.value(), layerdisk1_, trpdata_.stub_->isPSmodule())
                                       << " "
                                       << settings_->benddecode(outerbend.value(), layerdisk2_, outervmstub.isPSmodule())
                                       << " "
                                       << settings_->benddecode(innerbend.value(), layerdisk3_, innervmstub.isPSmodule());
        }
      } else {

        if (applyPairCut && lut_ok){ 
          // cut on triplet 
          if ( (rzbin_out < rzbinfirst_out_new || rzbin_out - rzbinfirst_out_new > rzdiffmax_out_new)){
            if (settings_->debugTracklet()) {
              edm::LogVerbatim("Tracklet") << "Outer stub rejected because of wrong r/z bin from pair";
            }
          }
          else {
            // passes pT cuts and triplet cut
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

            myTriplet.setStubPhi(0, innervmstub.stub()->phiapprox(0., 0));
            myTriplet.setStubPhi(1, trpdata_.stub_->phiapprox(0., 0));
            myTriplet.setStubPhi(2, outervmstub.stub()->phiapprox(0., 0));
          
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

          } // does not pass triplet cut
        } else { // triplet cut couldn't be defined, so still keep the candidate      
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
  
            myTriplet.setStubPhi(0, innervmstub.stub()->phiapprox(0, 0));
            myTriplet.setStubPhi(1, trpdata_.stub_->phiapprox(0, 0));
            myTriplet.setStubPhi(2, outervmstub.stub()->phiapprox(0, 0));
            
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
        std::tie(next_out_, outmem_, nstub_out_, phi_out_) = trpdata_.projbin_out_[nproj_out_];
      }
    }
    // get next in proj bin
    std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[nproj_in_];
  }
}
