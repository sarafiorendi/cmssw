#include "L1Trigger/TrackFindingTracklet/interface/TripletEngineUnit.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1StubTripletBuilder.h"

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
    : settings_(settings), 
      pttablemiddlenew_(pttablemiddle),
      pttableouternew_(pttableouter),
      pttablemiddleinnew_(pttablemiddlein),
      pttableinnernew_(pttableinner),
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

  if (trpdata_.bad_inner_stub_[(inmem_ * trpdata_.max_nbins_ + trpdata_.start_in_ + next_in_) * trpdata_.max_nstubs_ + istub_in_])
  {
//     std::cout << "  ####### skipping the stub at " << (inmem_ * trpdata_.max_nbins_ + trpdata_.start_in_ + next_in_) * trpdata_.max_nstubs_ + istub_in_ << " #########" << std::endl;
//     std::cout << "  inmem_" <<  inmem_  << std::endl;
//     std::cout << "  trpdata_.max_nbins_" <<  trpdata_.max_nbins_  << std::endl;
//     std::cout << "  trpdata_.start_in_" <<  trpdata_.start_in_  << std::endl;
//     std::cout << "  next_in_" <<  next_in_  << std::endl;
//     std::cout << "  trpdata_.max_nstubs_" <<  trpdata_.max_nstubs_  << std::endl;
//     std::cout << "  istub_in_" <<  istub_in_  << std::endl;
    // advance all indices
    advanceIndices();
    if (idle_) return;
  }
  
//   std::cout << "delta phi regions = " << phi_out_ - phi_in_ << " : " << phi_out_ << " : " << phi_in_ << std::endl;
  // I should get the minimum and max phi values of the inner and outer mem
  // if their max diff is larger than 0.22, just skip them
  // is this correct??
  int my_max_diff = 1000; // was 32  // 24 means some good seeds are removed; // 40 does nothing
  if (iSeed_ == 8){

    int tmp_phi_in = phi_in_;
    if (phi_in_ >= 2 * settings_->nvmte(1, iSeed_))
      tmp_phi_in = phi_in_ - 24;

    int min_innerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(1))) +   // 1 was innervmstub.stub()->layerdisk()
                           tmp_phi_in * (1 << settings_->nfinephi(1, iSeed_))  ;
    int max_innerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(1))) +   // 1 was innervmstub.stub()->layerdisk()
                           tmp_phi_in * (1 << settings_->nfinephi(1, iSeed_)) +  (1 << settings_->nfinephi(1, iSeed_));

    int tmp_phi_out = (phi_out_ >= 2 * settings_->nvmte(1, iSeed_)) ? phi_out_ - 24 : phi_out_;

    int min_outerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(1))) +   
                           tmp_phi_out * (1 << settings_->nfinephi(1, iSeed_))  ;
    int max_outerfinephi = iAllStub_ * (1 << (nbitsfinephi_ - settings_->nbitsallstubs(1))) +   
                           tmp_phi_out * (1 << settings_->nfinephi(1, iSeed_)) +  (1 << settings_->nfinephi(1, iSeed_));
//     std::cout << min_outerfinephi - max_innerfinephi << " : " << max_outerfinephi - min_innerfinephi  << std::endl ; 
    if ((min_outerfinephi - max_innerfinephi > my_max_diff) || (max_outerfinephi - min_innerfinephi > my_max_diff) )  {   
//       std::cout << "it happens " << min_outerfinephi - max_innerfinephi  << "    " <<  max_outerfinephi - min_innerfinephi << std::endl;
      // should mark this whole inner mem as bad and skip to next one               
      advanceIndicesInnerMem();  
    }  
  }  
  
  
  // try to see if new cut is useful
  if (iSeed_ == 8){
//   std::cout << nproj_out_ << "   but size is " << trpdata_.start_in_vec_.size() << std::endl;
    if (trpdata_.start_in_vec_[nproj_out_] > trpdata_.start_in_){
      std::cout << "new LUT is tighter: now " << trpdata_.start_in_vec_[nproj_out_] << "  instead of " << trpdata_.start_in_ << std::endl;
       advanceIndicesInnerMem();
      if (idle_) return;
    }
  }
    
  // get inner and outer projected stub for certain next value
  int ibin_out = trpdata_.start_out_ + next_out_;
  int ibin_in = trpdata_.start_in_ + next_in_;
  const VMStubTE& outervmstub = outervmstubs_[outmem_]->getVMStubTEBinned(ibin_out, istub_out_);
  const VMStubTE& innervmstub = innervmstubs_[inmem_]->getVMStubTEBinned(ibin_in, istub_in_);

//   std::cout << "loading from mem: " << trpdata_.start_in_ << " " << trpdata_.rzbinfirst_in_ << std::endl;

//   if (iSeed_ == 8) std::cout << "\t\t [TEU/" << trpdata_.istub_middle << "/" << outmem_ << "/"  << inmem_ 
//                              << "] ibin_out/istub_out_ "        << ibin_out << "/" << istub_out_ 
//                              << "           ibin_in/istub_in_ " << ibin_in << "/" << istub_in_ 
//                              << std::endl;

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
  if (!inrange_in){
    trpdata_.bad_inner_stub_[(inmem_ * trpdata_.max_nbins_ + trpdata_.start_in_ + next_in_) * trpdata_.max_nstubs_ + istub_in_] = true;
//     std::cout << " setting bad stub " << istub_in_ << std::endl;
  }  

  bool opposite_sign_dphi = (idphi_out * idphi_in > 0) ? 0 : 1;
  if (abs(idphi_out) < 10 && abs(idphi_in) < 10) opposite_sign_dphi = true;
//   if (!opposite_sign_dphi) std::cout << "\t will fail dphi "<< std::endl;
//   if (abs(idphi_out) < 10 && abs(idphi_in) < 10) opposite_sign_dphi = true;
  int idphi_in_for_index = idphi_in & ((1 << nbitsfinephiinnerdiff_) - 1);
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
//       std::cout << "setting stub " << istub_in_ << " as bad / at " << (inmem_ * trpdata_.max_nbins_ + trpdata_.start_in_ + next_in_) * trpdata_.max_nstubs_ + istub_in_ << std::endl;
//       std::cout << "\t  inmem_" <<  inmem_  << std::endl;
//       std::cout << "\t  trpdata_.max_nbins_" <<  trpdata_.max_nbins_  << std::endl;
//       std::cout << "\t  trpdata_.start_in_" <<  trpdata_.start_in_  << std::endl;
//       std::cout << "\t  next_in_" <<  next_in_  << std::endl;
//       std::cout << "\t  trpdata_.max_nstubs_" <<  trpdata_.max_nstubs_  << std::endl;
// 
// //       if (trpdata_.bad_inner_stub_[istub_in_] == true)
// //         std::cout << "wait it was already bad, should not happen " << std::endl;
      trpdata_.bad_inner_stub_[(inmem_ * trpdata_.max_nbins_ + trpdata_.start_in_ + next_in_) * trpdata_.max_nstubs_ + istub_in_] = true;
    } else {  // condition on both inner and outer stubs satisfied
    
      FPGAWord outerbend = outervmstub.bend();
      FPGAWord innerbend = innervmstub.bend();
      bool pass_pt_cut = false;
      if (iSeed_ == 8){
        int ptouterindex = (idphi_out_for_index << outerbend.nbits()) + outerbend.value();
        int ptmiddleindex = (idphi_out_for_index << trpdata_.middlebend_.nbits()) + trpdata_.middlebend_.value();
        pass_pt_cut = pttablemiddlenew_->lookup(ptmiddleindex) && pttableouternew_->lookup(ptouterindex);
      }
      if (pass_pt_cut){
//         std::cout << "passed pt cut! " << std::endl;
      }
      bool pass_pt_cut_inner = false;
      if (pass_pt_cut_inner){}
      if (iSeed_ == 8){
        int ptinnerindex = (idphi_in_for_index << innerbend.nbits()) + innerbend.value();
        int ptmiddleinindex = (idphi_in_for_index << trpdata_.middlebend_.nbits()) + trpdata_.middlebend_.value();
        pass_pt_cut_inner = pttablemiddleinnew_->lookup(ptmiddleinindex) && pttableinnernew_->lookup(ptinnerindex);
        if (!pass_pt_cut_inner)
          trpdata_.bad_inner_stub_[(inmem_ * trpdata_.max_nbins_ + trpdata_.start_in_ + next_in_) * trpdata_.max_nstubs_ + istub_in_] = true;
      }
      // now apply cut on pT
      //set all to true to only apply phi region cut
//       pass_pt_cut = true;
//       pass_pt_cut_inner = true;
//       inrange = true;
//       inrange_in = true;
//       opposite_sign_dphi = true;
      if (iSeed_ == 8 && ( !(inrange && inrange_in && opposite_sign_dphi && pass_pt_cut && pass_pt_cut_inner) )) {

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

            L1StubTriplet myTriplet = makeL1StubTriplet(
              innervmstub.stub(), trpdata_.stub_, outervmstub.stub(), 
              iSector, iTC, count_trpunits,
              trpdata_.rzbinfirst_out_, trpdata_.rzbinfirst_in_, rzbinfirst_out_new,
              trpdata_.rzdiffmax_out_, trpdata_.rzdiffmax_in_, rzdiffmax_out_new,
              rzbin_out, rzbin_in, ibin_out, ibin_in, 
              innervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1),
              outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1)
              );
          
            foundtriplets.push_back(myTriplet);

          } // does not pass triplet cut
        } else { // triplet cut couldn't be defined, so still keep the candidate      
            candtriplet_ =
                std::tuple<const Stub*, const Stub*, const Stub*>(innervmstub.stub(), trpdata_.stub_, outervmstub.stub());
            goodtriplet_ = true;

            L1StubTriplet myTriplet = makeL1StubTriplet(
              innervmstub.stub(), trpdata_.stub_, outervmstub.stub(), 
              iSector, iTC, count_trpunits,
              trpdata_.rzbinfirst_out_, trpdata_.rzbinfirst_in_, rzbinfirst_out_new,
              trpdata_.rzdiffmax_out_, trpdata_.rzdiffmax_in_, rzdiffmax_out_new,
              rzbin_out, rzbin_in, ibin_out, ibin_in, 
              innervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1),
              outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1)
            );
            foundtriplets.push_back(myTriplet);
        }    
      }  
    }
  }
  
  advanceIndices();
  if (idle_) return;

//   std::cout << "\t\t\t TEU maybe end of step, istub_in/out = " << istub_in_  << "  " << istub_out_ << std::endl;

}


void TripletEngineUnit::advanceIndices() {
  istub_in_++;
  if (istub_in_ >= nstub_in_) {
    nproj_in_++;
    istub_in_ = 0;
    if (nproj_in_ >= trpdata_.projbin_in_.size()) {
      istub_out_++;
      nproj_in_ = 0;
      if (istub_out_ >= nstub_out_) {
        nproj_out_++;
        istub_out_ = 0;
        if (nproj_out_ >= trpdata_.projbin_out_.size()) {
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
    std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[nproj_in_];
  }
}


void TripletEngineUnit::advanceIndicesInnerMem() {
    
  nproj_in_++;
  istub_in_ = 0;
  if (nproj_in_ >= trpdata_.projbin_in_.size()) {
    istub_out_++;
    nproj_in_ = 0;
    if (istub_out_ >= nstub_out_) {
      nproj_out_++;
      istub_out_ = 0;
      if (nproj_out_ >= trpdata_.projbin_out_.size()) {
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
  std::tie(next_in_, inmem_, nstub_in_, phi_in_) = trpdata_.projbin_in_[nproj_in_];
}