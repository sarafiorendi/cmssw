#include "L1Trigger/TrackFindingTracklet/interface/TrackletProcessorDisplaced.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllInnerStubsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1StubTripletBuilder.h"

#include <utility>
#include <tuple>

using namespace std;
using namespace trklet;

// TrackletProcessorDisplaced
//
// This module takes in collections of stubs within a phi region and a
// displaced seed name and tries to create that displaced seed out of the stubs
//
// Update: Claire Savard, Nov. 2024

TrackletProcessorDisplaced::TrackletProcessorDisplaced(string name, Settings const& settings, Globals* globals)
    : TrackletCalculatorDisplaced(name, settings, globals),
      trpbuffer_(CircularBuffer<TrpEData>(3), 0, 0, 0, 0),
      innerTable_(settings),
      innerThirdTable_(settings),
      outerPairTable_(settings),
      innerPairTable_(settings),
      testOuterInnerTable_(settings),
      useOuterRegiontable_(settings),
      useInnerRegiontable_(settings),
      pttablemiddle_(settings),
      pttableouter_(settings),
      pttablemiddlein_(settings),
      pttableinner_(settings),
      pttablemiddle_region_out_(settings),
      pttablemiddle_region_in_(settings) {

  innerallstubs_.clear();
  middleallstubs_.clear();
  outerallstubs_.clear();
  innervmstubs_.clear();
  outervmstubs_.clear();

  iAllStub_ = -1;
  // set layer/disk types based on input seed name
  initLayerDisksandISeedDisp(layerdisk1_, layerdisk2_, layerdisk3_, iSeed_);

  // as done in the prompt version, find min and max radii
  double rmiddle = -1.0;
  double router  = -1.0;
  double rinner  = -1.0;

  if (iSeed_ == Seed::L2L3L4 || iSeed_ == Seed::L4L5L6) {
    rmiddle = settings_.rmean(layerdisk1_); // middle
    router  = settings_.rmean(layerdisk2_); // outer
    rinner  = settings_.rmean(layerdisk3_); // inner
  } //else {
//     if (iSeed_ == Seed::L2L3D1) {
//       rmax = settings_.rmaxdiskvm();
//       rmin = settings_.rmean(layerdisk1_);
//     } else if (iSeed_ == Seed::L2D1) {
//       rmax = settings_.rmaxdiskvm();
//       rmin = settings_.rmean(layerdisk1_);
//     } else {
//       rmax = settings_.rmaxdiskvm();
//       rmin = rmax * settings_.zmean(layerdisk2_ - N_LAYER - 1) / settings_.zmean(layerdisk2_ - N_LAYER);
//     }
//   }
  // seed 11 missing

  // rmin and rmax correspond to the radius of the two layers/disks of the seed 
  // dphimax (in radians) is the max delta phi between the two stubs of the pair to form that seed in order to satisfy the min pt cut (maxrinv)
  // double dphimax = asin(0.5 * settings_.maxrinv() * rmax) - asin(0.5 * settings_.maxrinv() * rmin);
  // I'll replace this with my calculation
  double rinv = settings_.rinvmaxDisplaced(); // assuming charge = 1
  std::vector<float> d0_vals = linspace(-10,10,101);
  
  double dphimax_out = compute_dphimax(d0_vals, rinv, rmiddle, router);
  double dphimax_in = compute_dphimax(d0_vals, rinv, rinner, rmiddle);
  double dphimax_pair = compute_dphimax(d0_vals, rinv, rinner, router);

  // number of fine phi bins in sector for the outer/inner layer
  nbitsfinephiouterdiff_ = compute_nfinephibins(layerdisk2_,iSeed_, dphimax_out);
  nbitsfinephiinnerdiff_ = compute_nfinephibins(layerdisk3_,iSeed_, dphimax_in);
//   nbitsfinephipairdiff_ = compute_nfinephibins(layerdisk3_,iSeed_, dphimax_pair); // not used

  nbitsfinephi_ = settings_.nbitsallstubs(layerdisk2_) + settings_.nbitsvmte(1, iSeed_) + settings_.nfinephi(1, iSeed_);
  
//   std::cout << "[TPD:Seed " << iSeed_ << "] max dphi middle-outer = " << dphimax_out 
//             << "   max dphi middle-inner = " << dphimax_in 
//             << " (r_middle, r_outer) " << rmiddle << "," << router << std::endl;
// 
//   std::cout << "[TPD] nfinephibins = " << nfinephibins_out  // should not depend on the cut
//             << "   dfinephi " << dfinephi_out               // should not depend on the cut 
//             << "   nbitsfinephi_ " << nbitsfinephi_         // does not depend on the cut
//             << "   nbitsfinephiouterdiff_ " << nbitsfinephiouterdiff_  // this should be the cut
//             << "   nbitsfinephiinnerdiff_ " << nbitsfinephiinnerdiff_  // this should be the cut
//             << std::endl;
    
//   if (iSeed_ == 10) std::cout << "ld123:" <<  layerdisk1_ << " " << layerdisk2_ << " " << layerdisk3_ << std::endl;
//   ld123:1 2 6

  // get projection tables
  unsigned int region = name.back() - 'A';
  innerTable_.initVMRTable(
      layerdisk1_, TrackletLUT::VMRTableType::inner, region, false);  //projection to next layer/disk
  innerThirdTable_.initVMRTable(
      layerdisk1_, TrackletLUT::VMRTableType::innerthird, region, false);  //projection to third layer/disk
  outerPairTable_.initVMRTableTriplet(
      layerdisk1_, layerdisk3_, TrackletLUT::VMRTableType::outerfrompair, region, false);  //projection to outer layer/disk from inner pair
  innerPairTable_.initVMRTableTriplet(
      layerdisk2_, layerdisk1_, TrackletLUT::VMRTableType::outerfrompair, region, false);  //projection to outer layer/disk from inner pair
  
  testOuterInnerTable_.initVMRTable(
      layerdisk2_, TrackletLUT::VMRTableType::innerthird, region, false);  //projection to third layer/disk      

//   if ((layerdisk1_== 1 && layerdisk2_== 2 && layerdisk3_== 6))
//     std::cout << "size of lut " << outerPairTable_.size() << std::endl;
  nbitszfinebintable_ = settings_.vmrlutzbits(layerdisk1_);
  nbitsrfinebintable_ = settings_.vmrlutrbits(layerdisk1_);

  for (unsigned int ilayer = 0; ilayer < N_LAYER; ilayer++) {
    vector<TrackletProjectionsMemory*> tmp(settings_.nallstubs(ilayer), nullptr);
    trackletprojlayers_.push_back(tmp);
  }

  for (unsigned int idisk = 0; idisk < N_DISK; idisk++) {
    vector<TrackletProjectionsMemory*> tmp(settings_.nallstubs(idisk + N_LAYER), nullptr);
    trackletprojdisks_.push_back(tmp);
  }

  // set TC index
  iTC_ = region;
  TCIndex_ = (iSeed_ << settings.nbitsseed()) + iTC_;

  maxStep_ = settings_.maxStep("TPD");
//   maxStep_ = 10;
}

void TrackletProcessorDisplaced::addOutputProjection(TrackletProjectionsMemory*& outputProj, MemoryBase* memory) {
  outputProj = dynamic_cast<TrackletProjectionsMemory*>(memory);
  assert(outputProj != nullptr);
}

void TrackletProcessorDisplaced::addOutput(MemoryBase* memory, string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }

  if (output == "trackpar") {
    auto* tmp = dynamic_cast<TrackletParametersMemory*>(memory);
    assert(tmp != nullptr);
    trackletpars_ = tmp;
    return;
  }

  if (output.substr(0, 7) == "projout") {
    //output is on the form 'projoutL2PHIC' or 'projoutD3PHIB'
    auto* tmp = dynamic_cast<TrackletProjectionsMemory*>(memory);
    assert(tmp != nullptr);

    constexpr unsigned layerdiskPosInprojout = 8;
    constexpr unsigned phiPosInprojout = 12;

    unsigned int layerdisk = output[layerdiskPosInprojout] - '1';  //layer or disk counting from 0
    unsigned int phiregion = output[phiPosInprojout] - 'A';        //phiregion counting from 0

    if (output[7] == 'L') {
      assert(layerdisk < N_LAYER);
      assert(phiregion < trackletprojlayers_[layerdisk].size());
      //check that phiregion not already initialized
      assert(trackletprojlayers_[layerdisk][phiregion] == nullptr);
      trackletprojlayers_[layerdisk][phiregion] = tmp;
      return;
    }

    if (output[7] == 'D') {
      assert(layerdisk < N_DISK);
      assert(phiregion < trackletprojdisks_[layerdisk].size());
      //check that phiregion not already initialized
      assert(trackletprojdisks_[layerdisk][phiregion] == nullptr);
      trackletprojdisks_[layerdisk][phiregion] = tmp;
      return;
    }
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find output : " << output;
}

void TrackletProcessorDisplaced::addInput(MemoryBase* memory, string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }

  if (input == "thirdallstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    innerallstubs_.push_back(tmp);
    return;
  }
  if (input == "firstallstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    middleallstubs_.push_back(tmp);

    iAllStub_ = tmp->getName()[8] - 'A';
    unsigned int iTP = iAllStub_;
//     std::cout << "[TPD" << getName() << "] middlestub from " << tmp->getName()
//             << " will initiate table for seed " << iSeed_ 
//             << "   iAllStub_ = " << iAllStub_ 
//             << "   (from " << tmp->getName()[8] 
// // //             << "and nfinephibins = " << nfinephibins_out  // should not depend on the cut
// // //             << "   dfinephi " << dfinephi_out               // should not depend on the cut 
// //             << "   nbitsfinephi_ " << nbitsfinephi_         // does not depend on the cut
//             << "   nbitsfinephiouterdiff_ " << nbitsfinephiouterdiff_  // this should be the cut
//             << std::endl;

    pttablemiddle_region_out_.initTPlut(true, iSeed_, layerdisk1_, layerdisk2_, nbitsfinephiouterdiff_, iTP);
    useOuterRegiontable_.initDisplacedOuterTPregionlut(
      iSeed_, layerdisk1_, layerdisk2_, iAllStub_, nbitsfinephiouterdiff_, nbitsfinephi_, pttablemiddle_region_out_, iTP);
      // iTP is only used in the name of the LUT table, if written out
      // iAllStub is used to define the outerfinephi
    pttablemiddle_region_in_.initDisplacedTPlutForInner(true, iSeed_, layerdisk1_, layerdisk3_, nbitsfinephiinnerdiff_, iTP);
    useInnerRegiontable_.initDisplacedOuterTPregionlut(
      iSeed_, layerdisk1_, layerdisk3_, iAllStub_, nbitsfinephiinnerdiff_, nbitsfinephi_, pttablemiddle_region_in_, iTP);

    return;
  }
  if (input == "secondallstubin") {
    auto* tmp = dynamic_cast<AllStubsMemory*>(memory);
    assert(tmp != nullptr);
    outerallstubs_.push_back(tmp);
    return;
  }

  if (input == "thirdvmstubin") {
    auto* tmp = dynamic_cast<VMStubsTEMemory*>(memory);
    assert(tmp != nullptr);
    innervmstubs_.push_back(tmp);

    unsigned int iTP = getName()[7] - 'A';
    pttablemiddlein_.initDisplacedTPlutForInner(true, iSeed_, layerdisk1_, layerdisk3_, nbitsfinephiinnerdiff_, iTP);
    pttableinner_.initDisplacedTPlutForInner(false, iSeed_, layerdisk1_, layerdisk3_, nbitsfinephiinnerdiff_, iTP);
    return;
  }
  
  if (input == "secondvmstubin") {
    auto* tmp = dynamic_cast<VMStubsTEMemory*>(memory);
    assert(tmp != nullptr);
    outervmstubs_.push_back(tmp);

    unsigned int iTP = getName()[7] - 'A';
//     std::cout << "[" << getName() << "] outerstub from " << tmp->getName()
//               << " will initiate table for seed " << iSeed_ 
//               << "   nbitsfinephiouterdiff_ " << nbitsfinephiouterdiff_  // this should be the cut
//               << std::endl;
    pttablemiddle_.initTPlut(true, iSeed_, layerdisk1_, layerdisk2_, nbitsfinephiouterdiff_, iTP);
    pttableouter_.initTPlut(false, iSeed_, layerdisk1_, layerdisk2_, nbitsfinephiouterdiff_, iTP);
    return;
  }

  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " Could not find input : " << input;
}

void TrackletProcessorDisplaced::execute(unsigned int iSector, double phimin, double phimax, std::vector<L1StubTriplet>& foundtriplets_, std::vector<L1StubTriplet>& acceptedtriplets_) {
  phimin_ = phimin;
  phimax_ = phimax;
  iSector_ = iSector;

  unsigned int countall = 0;
  unsigned int countsel = 0;
  int donecount = 0;
  
//   if (iSeed_ == 8){
//   std::cout << getName() 
//             << "   iAllStub_ = " << iAllStub_ 
//             << std::endl;
//   }
//   std::cout << "TPD: iSector = " << iSector 
//             << "\nTPD: iSeed_ = " << iSeed_ 
//             << "\nTPD: iAllStub_ = " << iAllStub_ 
//             << "\nTPD: iTC_ = " << iTC_ 
//             << "\nTPD: TCIndex_ = " << TCIndex_ 
//             << std::endl;

  // set the triplet engine units and buffer
  TripletEngineUnit trpunit(&settings_, 
                            layerdisk1_, 
                            layerdisk2_, 
                            layerdisk3_, 
                            iSeed_, 
                            iAllStub_,
                            nbitsfinephi_,
                            nbitsfinephiouterdiff_,
                            nbitsfinephiinnerdiff_,
                            &pttablemiddle_,
                            &pttableouter_,
                            &pttablemiddlein_,
                            &pttableinner_,
                            innervmstubs_, 
                            outervmstubs_);
  trpunits_.resize(settings_.trpunits(iSeed_), trpunit);
  trpbuffer_ = tuple<CircularBuffer<TrpEData>, unsigned int, unsigned int, unsigned int, unsigned int>(
      CircularBuffer<TrpEData>(3), 0, 0, 0, middleallstubs_.size());

  // reset the trpunits
  for (auto& trpunit : trpunits_) {
    trpunit.reset();
  }

  // reset the Circular buffer (first element of trpbuffer_)
  std::get<0>(trpbuffer_).reset();
  // reset the first int to 0
  std::get<1>(trpbuffer_) = 0;
  // copy fourth element into the third
  std::get<2>(trpbuffer_) = std::get<3>(trpbuffer_);

  TrpEData trpdata;
  TrpEData trpdata__;
  TrpEData trpdata___;
  bool goodtrpdata = false;
  bool goodtrpdata__ = false;
  bool goodtrpdata___ = false;

  bool trpbuffernearfull;
  for (unsigned int istep = 0; istep < maxStep_; istep++) {
    CircularBuffer<TrpEData>& trpdatabuffer = std::get<0>(trpbuffer_);
    trpbuffernearfull = trpdatabuffer.nearfull();
//     if (iSeed_ == 8) std::cout << "istep: " << istep << " , nearfull " << trpbuffernearfull << std::endl;

    //
    // First block here checks if there is a trpunit with data that should be used
    // to calculate the tracklet parameters
    //

    // set pointer to the last filled trpunit
    TripletEngineUnit* trpunitptr = nullptr;
    int count_trpunits_block1 = 0;
    int the_trpunit_being_read = 0;    
    for (auto& trpunit : trpunits_) {
      trpunit.setNearFull();
      if (!trpunit.empty()) {
        trpunitptr = &trpunit;
        the_trpunit_being_read = count_trpunits_block1;
//         if (iSeed_ == 8) std::cout << "first block, reading trpunit " << count_trpunits_block1 << std::endl;
      }
      count_trpunits_block1++;      
    }
  
    bool studyTrunc = false;
    bool do_block1 = true;
    if (studyTrunc)
      if (iSeed_ == 8 && istep >= 108) do_block1 = false;
    if (do_block1){
      if (trpunitptr != nullptr) {
        auto stubtriplet = trpunitptr->read();
  
        countall++;
  
        const Stub* innerFPGAStub = std::get<0>(stubtriplet);
        const Stub* middleFPGAStub = std::get<1>(stubtriplet);
        const Stub* outerFPGAStub = std::get<2>(stubtriplet);
  
        const L1TStub* innerStub = innerFPGAStub->l1tstub();
        const L1TStub* middleStub = middleFPGAStub->l1tstub();
        const L1TStub* outerStub = outerFPGAStub->l1tstub();
  
        if (settings_.debugTracklet()) {
          edm::LogVerbatim("Tracklet") << "TrackletProcessorDisplaced execute " << getName() << "[" << iSector_ << "]";
        }
  //       if (iSeed_ == 8) std::cout << "first block, checking one triplet"<< std::endl;
  
        L1StubTriplet myTripletNow = makeL1StubTriplet(innerFPGAStub, middleFPGAStub, outerFPGAStub, iSector, iTC_, the_trpunit_being_read);
//         myTripletNow.setStubRapprox(0, innerFPGAStub->rapprox());
//         myTripletNow.setStubRapprox(1, middleFPGAStub->rapprox());
//         myTripletNow.setStubRapprox(2, outerFPGAStub->rapprox());
//   
//         myTripletNow.setStubZapprox(0, innerFPGAStub->zapprox());
//         myTripletNow.setStubZapprox(1, middleFPGAStub->zapprox());
//         myTripletNow.setStubZapprox(2, outerFPGAStub->zapprox());
//   
//         myTripletNow.setStubBend(0, innerFPGAStub->bend().value());
//         myTripletNow.setStubBend(1, middleFPGAStub->bend().value());
//         myTripletNow.setStubBend(2, outerFPGAStub->bend().value());
//   
//         myTripletNow.setStubPhi(0, innerFPGAStub->phiapprox(0., 0));
//         myTripletNow.setStubPhi(1, middleFPGAStub->phiapprox(0., 0));
//         myTripletNow.setStubPhi(2, outerFPGAStub->phiapprox(0., 0));
//   
//         myTripletNow.setStubIndex(0, innerFPGAStub->stubindex().value());
//         myTripletNow.setStubIndex(1, middleFPGAStub->stubindex().value());
//         myTripletNow.setStubIndex(2, outerFPGAStub->stubindex().value());
//   
//         myTripletNow.setStubLayerdisk(0, innerFPGAStub->layerdisk());
//         myTripletNow.setStubLayerdisk(1, middleFPGAStub->layerdisk());
//         myTripletNow.setStubLayerdisk(2, outerFPGAStub->layerdisk());
//   
//         myTripletNow.setSector(iSector);
//         myTripletNow.setRegion(iTC_);
//         myTripletNow.setTPDUnit(the_trpunit_being_read);
        
        // check if the seed made from the 3 stubs is valid
        bool accept = false;
        if (iSeed_ == Seed::L2L3L4 || iSeed_ == Seed::L4L5L6)
          accept = LLLSeeding(innerFPGAStub, innerStub, middleFPGAStub, middleStub, outerFPGAStub, outerStub);
        else if (iSeed_ == Seed::L2L3D1)
          accept = LLDSeeding(innerFPGAStub, innerStub, middleFPGAStub, middleStub, outerFPGAStub, outerStub);
        else if (iSeed_ == Seed::D1D2L2)
          accept = DDLSeeding(innerFPGAStub, innerStub, middleFPGAStub, middleStub, outerFPGAStub, outerStub);
  
        if (accept){
          acceptedtriplets_.push_back(myTripletNow);
          countsel++;
  //         if (iSeed_ == 8) std::cout << "\t accepted "  << countsel << std::endl;
        }  
  
        if (trackletpars_->nTracklets() >= settings_.ntrackletmax()) {
  //         std::cout << "Will break on number of tracklets in " << getName() <<  "  max is " << settings_.ntrackletmax() << std::endl;
          edm::LogVerbatim("Tracklet") << "Will break on number of tracklets in " << getName();
          assert(0);
          break;
        }
  
        if (settings_.debugTracklet()) {
          edm::LogVerbatim("Tracklet") << "TrackletProcessor execute done";
        }
      }
    }

    //
    // The second block fills the trpunit if data in buffer and process TripletEngineUnit step
    //
    //

    bool notemptytrpbuffer = !trpdatabuffer.empty();
    int count_trpunits = 0;    
    bool do_block2 = true;
//     if (studyTrunc)
//       if (iSeed_ == 8 && istep >= 108) do_block2 = false;
    if (do_block2){
//         if (iSeed_ == 8 ) {
//           if (!notemptytrpbuffer)
//             std::cout << "buffer is empty" << std::endl;
//           else
//             std::cout << "buffer has data" << std::endl;
//         }
        for (auto& trpunit : trpunits_) {
          if (trpunit.idle() && notemptytrpbuffer) {  // only fill one idle unit every step
            trpunit.init(std::get<0>(trpbuffer_).read());
            notemptytrpbuffer = false;  //prevent initializing another triplet engine unit
          }
//           if (iSeed_ == 8) std::cout << "  " << iTC_ << "," << istep << "," << count_trpunits << std::endl;
          trpunit.step(foundtriplets_, iSector, iTC_, count_trpunits);
          count_trpunits++;      
        }
//         if (iSeed_ == 8 && istep == (108 - 1) && !trpdatabuffer.empty()) 
//             std::cout << "will truncate block2: left " << std::endl;
    }

    //
    // The third block here checks if we have input stubs to process
    //
    //

    if (goodtrpdata___){
      trpdatabuffer.store(trpdata___);
//       if (iSeed_==8) std::cout << "finally storing into buffer" << std::endl;
    }  
    goodtrpdata = false;

    unsigned int& istub = std::get<1>(trpbuffer_);
    unsigned int& midmem = std::get<2>(trpbuffer_);
    unsigned int midmemend = std::get<4>(trpbuffer_);

//     std::cout << "istep: " << istep << " , looking at middle stub " << istub << " out of " << middleallstubs_[midmem]->nStubs() << std::endl;
//     if (iSeed_ == 8)
//       std::cout << "\t midmem: " << midmem << std::endl;
    bool do_next = true;
//     if (studyTrunc)
//       if (istep >= 108) do_next = false;
    if (do_next && (!trpbuffernearfull) && midmem < midmemend && istub < middleallstubs_[midmem]->nStubs()) {

//       if (iSeed_ == 8) std::cout << iTC_ << "," << istep << "," << istub << "," << midmem << "," << middleallstubs_[midmem]->nStubs() << "," << midmemend << std::endl;
      
      const Stub* stub = middleallstubs_[midmem]->getStub(istub);

      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "In " << getName() << " have middle stub";
      }

      bool negdisk = (stub->disk().value() < 0);  // check if disk in negative z region

      // get r/z index of the middle stub
      int indexz = (((1 << (stub->z().nbits() - 1)) + stub->z().value()) >> (stub->z().nbits() - nbitszfinebintable_));
      int indexr = -1;
      if (layerdisk1_ >= LayerDisk::D1) {  // if projecting from a disk
        if (negdisk)
          indexz = (1 << nbitszfinebintable_) - indexz;
        indexr = stub->rvalue();
        if (stub->isPSmodule()) {
          indexr = stub->rvalue() >> (stub->r().nbits() + 1 - nbitsrfinebintable_);
        }
      } else {  // else if from a layer
        indexr = (((1 << (stub->r().nbits() - 1)) + stub->rvalue()) >> (stub->r().nbits() - nbitsrfinebintable_));
      }

      // create lookupbits that define projections from middle stub
      int lutval = -1;
      const auto& lutshift = innerTable_.nbits();
      lutval = innerTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
//       auto sara_lutval = innerTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
      int lutval2 = innerThirdTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
//       if (lutval == -1 && iSeed_ == 8){
//         std::cout << "outer lut -1: z/r = " << stub->zapprox() << "/" << stub->rapprox() << std::endl;
//         std::cout << "index z/r = " << indexz << "/" <<indexr  << " -> " << (indexz << nbitsrfinebintable_) + indexr << std::endl;
//         std::cout << "lutval = " << sara_lutval 
//                   << " " << std::bitset<16>(sara_lutval) << std::endl;
//       }
      if (lutval != -1 && lutval2 != -1)
        lutval += (lutval2 << lutshift);

      if (lutval != -1) {
        unsigned int lutwidth = settings_.lutwidthtabextended(0, iSeed_);
        FPGAWord lookupbits(lutval, lutwidth, true, __LINE__, __FILE__);

        // get r/z bins for projection into outer layer/disk
        int nbitsrzbin_out = N_RZBITS;  //number of bits for the r/z bins
        if (iSeed_ == Seed::D1D2L2)
          nbitsrzbin_out--;
        int rzbinfirst_out = lookupbits.bits(0, NFINERZBITS);  // first small rz-bin projection
        int next_out = lookupbits.bits(NFINERZBITS, 1);
        int start_out = lookupbits.bits(NFINERZBITS + 1, nbitsrzbin_out);  // first large rz-bin projection
        int rzdiffmax_out = lookupbits.bits(NFINERZBITS + 1 + nbitsrzbin_out, NFINERZBITS);
        if (iSeed_ == Seed::D1D2L2 && negdisk)  // if projecting into disk
          start_out += (1 << nbitsrzbin_out);
        int last_out = start_out + next_out;  // last rz bin projection

        // get r/z bins for projection into third (inner) layer/disk
        int nbitsrzbin_in = N_RZBITS;
        int rzbinfirst_in = lookupbits.bits(lutshift, NFINERZBITS);
        int next_in = lookupbits.bits(lutshift + NFINERZBITS, 1);
        int start_in = lookupbits.bits(lutshift + NFINERZBITS + 1, nbitsrzbin_in);  // first rz bin projection
        int rzdiffmax_in = lookupbits.bits(lutshift + NFINERZBITS + 1 + nbitsrzbin_in, NFINERZBITS);
        // LUT doesn't know about z sign.
        // So, first, mirror index of large z bin wrt center (as 0-3 bins are for negative z, 4 to 7 on the positive z).
        // Then subtract next_in so that we take that into account
        if (iSeed_ == Seed::D1D2L2 && negdisk) {  // if projecting from disk into layer in the negative z region
          start_in = settings_.NLONGVMBINS() - 1 - start_in - next_in;
          if (next_in)
            rzbinfirst_in = settings_.NLONGVMBINS() - (rzbinfirst_in + rzdiffmax_in - settings_.NLONGVMBINS());
          else
            rzbinfirst_in = settings_.NLONGVMBINS() - 1 - rzbinfirst_in - rzdiffmax_in;
          if (rzbinfirst_in < 0)
            rzbinfirst_in = 0;
        }
//         int last_in = start_in + next_in;  // last large rz-bin projection


        // trial inner from outer, first outer bin
        // negative z handling!!
        unsigned int n_rbins_out = (1 << settings_.vmrlutrbits(layerdisk2_));
        unsigned int n_zbins_out = (1 << settings_.vmrlutzbits(layerdisk2_));

//         if (iSeed_ == 8){
//             std::cout << "[TPD:middle]: r/z = " << stub->rapprox() << "/" << stub->zapprox() 
//                       << "  index z/r = " << indexz << "/" <<indexr  << " -> total index = " << (indexz << nbitsrfinebintable_) + indexr 
//                       << std::endl;
//         }
            // this is for the first possible rzbin of the outer stub
//             int zbin_out_index = start_out * (8) + rzbinfirst_out;
//             std::vector<int> min_max_rzindex_out = find_rzbin_from_zbin(zbin_out_index, n_zbins_out);
//             
//             int shift_innerPairLUT = settings_.vmrlutrbits(layerdisk2_) + settings_.vmrlutzbits(layerdisk2_);
//             int lutval_innerPair_minr = innerPairTable_.lookup(
//               (indexz << (settings_.vmrlutrbits(layerdisk2_) + shift_innerPairLUT)) + 
//               (indexr << shift_innerPairLUT) + 
//               min_max_rzindex_out[0]
//             );
//             int lutval_innerPair_maxr = innerPairTable_.lookup(
//               (indexz << (settings_.vmrlutrbits(layerdisk2_) + shift_innerPairLUT)) + 
//               (indexr << shift_innerPairLUT) + 
//               min_max_rzindex_out[1]
//             );
// 
//             if (lutval_innerPair_minr != -1){
//               std::vector<int> inner_rz_from_pair_minr = extract_rz_from_lut(lutval_innerPair_minr, lutwidth, 0, nbitsrzbin_in);
//               if (true){
//                 std::cout << "from PAIR: " << inner_rz_from_pair_minr[0] << "," << inner_rz_from_pair_minr[1] << "," << inner_rz_from_pair_minr[2] << ",\n" 
//                           << "from MIDD: " << start_in      << "," << rzbinfirst_in      << "," << rzdiffmax_in << std::endl ;
//               }
//             }
//             if (lutval_innerPair_maxr != -1){
//               std::vector<int> inner_rz_from_pair_maxr = extract_rz_from_lut(lutval_innerPair_maxr, lutwidth, 0, nbitsrzbin_in);
//               if (true){
//                 std::cout << "from PAIRMAX: " << inner_rz_from_pair_maxr[0] << "," << inner_rz_from_pair_maxr[1] << "," << inner_rz_from_pair_maxr[2] << std::endl ;
//               }
//             }
//             std::vector<int> updated_inner_ranges_from_pair = choose_rzbins_range(inner_rz_from_pair_minr, 
//                                                                                   inner_rz_from_pair_maxr, 
//                                                                                   start_in, 
//                                                                                   rzbinfirst_in,
//                                                                                   rzdiffmax_in
//                                                                                   );
//             
//             // find smallest possible range for this specific ibin_out
//             // save it into the trpdata as a vector (same size of the n_bin_out)
//             trpdata.start_in_vec_.push_back(updated_inner_ranges_from_pair[0]);
//             
            
            // this is inner range from outer + BS
//             int lutval_inner_from_out_minr = testOuterInnerTable_.lookup(min_max_rzindex_out[0]);
//             int lutval_inner_from_out_maxr = testOuterInnerTable_.lookup(min_max_rzindex_out[1]);
//             unsigned int lutwidth_out = settings_.lutwidthtabextended(0, iSeed_);
//             const auto& lutshift_out = testOuterInnerTable_.nbits();
//             std::vector<int> inner_rz_outer_min = extract_rz_from_lut(lutval_inner_from_out_minr, lutwidth_out, lutshift_out, nbitsrzbin_in);
//             std::vector<int> inner_rz_outer_max = extract_rz_from_lut(lutval_inner_from_out_maxr, lutwidth_out, lutshift_out, nbitsrzbin_in);
//             
//             std::vector<int> updated_inner_ranges = choose_rzbins_range(inner_rz_outer_min, 
//                                                                         inner_rz_outer_max, 
//                                                                         start_in, 
//                                                                         rzbinfirst_in,
//                                                                         rzdiffmax_in
//                                                                         );
//             // if I want to udpate the inner ranges using outer + BS                                                                                
// //             start_in = updated_inner_ranges[0];
// //             rzbinfirst_in = updated_inner_ranges[1];
// //             rzdiffmax_in = updated_inner_ranges[2];
// //             next_in = updated_inner_ranges[3];
//         }
        int last_in = start_in + next_in;  // last large rz-bin projection

        FPGAWord phicorr = stub->phicorr();
        // nbitsallstubs is the same for all layers except L1
        // the other two pars only depend on the seed, not on the layer 
        nbitsfinephi_ = settings_.nbitsallstubs(layerdisk2_) + settings_.nbitsvmte(1, iSeed_) + settings_.nfinephi(1, iSeed_);
        int middlefinephi = phicorr.bits(phicorr.nbits() - nbitsfinephi_, nbitsfinephi_);
        FPGAWord middlebend = stub->bend();

//         unsigned int useregindex = (middlefinephi); // valid for both the inner and outer LUTs as it depends only on the phi of the middle stub
        // use pt cut in phiregion
        unsigned int nbendbitsmiddle = 3;
        if (iSeed_ == Seed::L4L5L6) {
          nbendbitsmiddle = 4;
        }
        unsigned int useregindex = (middlefinephi << nbendbitsmiddle) + middlebend.value();

        int usereg_out = -1;
        usereg_out = useOuterRegiontable_.lookup(useregindex);
        int usereg_in = -1;
        usereg_in = useInnerRegiontable_.lookup(useregindex);

        // fill trpdata with projection info of middle stub
        trpdata.stub_ = stub;
        trpdata.istub_middle = istub; // only for understanding
        trpdata.rzbinfirst_out_ = rzbinfirst_out;
        trpdata.rzdiffmax_out_ = rzdiffmax_out;
        trpdata.rzbinfirst_in_ = rzbinfirst_in;
        trpdata.rzdiffmax_in_ = rzdiffmax_in;
        trpdata.start_out_ = start_out;
        trpdata.start_in_ = start_in;
        trpdata.outerpairtable_ = &outerPairTable_; 
        trpdata.middlefinephi_ = middlefinephi;
        trpdata.middlebend_ = middlebend;
        // fill projection bins info for single engine unit
        trpdata.projbin_out_.clear();
        trpdata.projbin_in_.clear();
        trpdata.start_in_vec_.clear();
        trpdata.bad_inner_stub_.clear();
        
//         std::cout << "[TPD] size of out memories to be read (outervmstubs_.size()) = " << outervmstubs_.size() << 
//                      "\n \t settings_.nvmte(1, iSeed_ = " << iSeed_ << ") = " << settings_.nvmte(1, iSeed_) << 
//                      std::endl;
       
        std::string mask = "";
        int count_nstubs_out_sara = 0;
        int count_nstubs_in_sara = 0;
        for (int ibin_out = start_out; ibin_out <= last_out; ibin_out++) {
//           if (iSeed_ == 8)
//             std::cout << "\t\t ibin_out " << ibin_out << std::endl;
          // look into all memories for the outer layer 
          // outervmstubs_.size() is the number of outer memories, as from the wiring file, e.g.
          // instance A has 9 outer memories, instance B -> 6, instance C -> 6, instance D -> 10...
          // each outer memory covers 1/settings_.nvmte(1, iSeed_) of one of the 4 large phi regions of the outer layer (1/8*1/4 of L4 for seed 8 for example)
          for (unsigned int outmem = 0; outmem < outervmstubs_.size(); outmem++) {
              // for each memory, check if its region is compatible 
            unsigned int out_phi_region = (outervmstubs_[outmem]->phibin() - 1) - (outervmstubs_[outmem]->getName()[11] - 'A') * 8;
//             double this_phimin, this_phimax;
//             outervmstubs_[outmem]->getPhiRange(this_phimin, this_phimax, iSeed_, 1);
//             if (iSeed_ == 8) std::cout << outervmstubs_[outmem]->getName() << " -> " 
//                                           " this_phimin = " << this_phimin << "," << 
//                                           " this_phimax = " << this_phimax << 
//                                           std::endl;
            
            if (iSeed_ == Seed::L2L3L4) {
              char mem_reg_out  = outervmstubs_[outmem]->getName()[11];
              char mem_reg_mid = middleallstubs_[midmem]->getName()[8];

              int diff_reg = mem_reg_out - mem_reg_mid;
              if (diff_reg == 1) {
                  out_phi_region += 8;
              } 
              else if (diff_reg == -1) {
                  out_phi_region += 16;
              } 
              if (std::abs(diff_reg) < 2){
                if (usereg_out != -1 && !(usereg_out & (1 << out_phi_region))) {
                  continue;
                }  
              }  
            } 
            int nstubs_out = outervmstubs_[outmem]->nVMStubsBinned(ibin_out);
            if (iSeed_ == 8) count_nstubs_out_sara += nstubs_out;
//             if (iSeed_ == 8) std::cout << iTC_ << "," << istub << "," << outmem << "," << nstubs_out << std::endl;
            if (nstubs_out > 0){
              mask = "1" + mask;
              
              
              // refine ibins of inner stub
              if (iSeed_ == 8) {
//                 std::cout << ibin_out << "," << start_out << "," << rzbinfirst_out << "," << start_in << "," << rzbinfirst_in  << std::endl;
                int zbin_out_index = start_out * (8) + rzbinfirst_out;
                std::vector<int> min_max_rzindex_out = find_rzbin_from_zbin(zbin_out_index, n_zbins_out);
                
//                 int shift_innerPairLUT = settings_.vmrlutrbits(layerdisk2_) + settings_.vmrlutzbits(layerdisk2_);
                int shift_innerPairLUT = settings_.vmrlutrbits(layerdisk2_) + settings_.vmrlutzbits(layerdisk2_);
                int lutval_innerPair_minr = innerPairTable_.lookup(
                   (min_max_rzindex_out[0] << shift_innerPairLUT) +  
                   (indexz << settings_.vmrlutrbits(layerdisk1_)) + 
                   indexr
                );
//                 std::cout << "outer index in new LUT = " << 
//                             ((min_max_rzindex_out[0] << shift_innerPairLUT) + (indexz << settings_.vmrlutrbits(layerdisk1_)) + indexr) << std::endl; 
                int lutval_innerPair_maxr = innerPairTable_.lookup(
                   (min_max_rzindex_out[1] << shift_innerPairLUT) +  
                   (indexz << settings_.vmrlutrbits(layerdisk1_)) + 
                   indexr
                );
    
                std::vector<int> inner_rz_from_pair_maxr{-1,-1,-1,-1};
                std::vector<int> inner_rz_from_pair_minr{-1,-1,-1,-1};
                
                if (lutval_innerPair_minr != -1){
//                   std::cout << "inner from PAIR_minr: " ; 
                  inner_rz_from_pair_minr = extract_rz_from_lut(lutval_innerPair_minr, lutwidth, 0, nbitsrzbin_in);
                  if (false){
                    std::cout << "from PAIR: " << inner_rz_from_pair_minr[0] << "," << inner_rz_from_pair_minr[1] << "," << inner_rz_from_pair_minr[2] << ",\n" 
                              << "from MIDD: " << start_in      << "," << rzbinfirst_in      << "," << rzdiffmax_in << std::endl ;
                  }
                }
                if (lutval_innerPair_maxr != -1){
//                   std::vector<int> 
//                   std::cout << "inner from PAIR_maxr: " ; 
                  inner_rz_from_pair_maxr = extract_rz_from_lut(lutval_innerPair_maxr, lutwidth, 0, nbitsrzbin_in);
                  if (false)
                    std::cout << "from PAIRMAX: " << inner_rz_from_pair_maxr[0] << "," << inner_rz_from_pair_maxr[1] << "," << inner_rz_from_pair_maxr[2] << std::endl ;
                }
                std::vector<int> updated_inner_ranges_from_pair = choose_rzbins_range(inner_rz_from_pair_minr, 
                                                                                      inner_rz_from_pair_maxr, 
                                                                                      start_in, 
                                                                                      rzbinfirst_in,
                                                                                      rzdiffmax_in,
                                                                                      next_in
                                                                                      );
//                   if (true)
//                     std::cout << "FINAL: " << updated_inner_ranges_from_pair[0] << "," << updated_inner_ranges_from_pair[1] << "," << updated_inner_ranges_from_pair[2] << std::endl ;
                
                // find smallest possible range for this specific ibin_out
                // save it into the trpdata as a vector (same size of the n_bin_out)
                trpdata.start_in_vec_.emplace_back(updated_inner_ranges_from_pair[0]);

                //check if i'm doing ok
//                 int lutval_inner_from_out_minr = testOuterInnerTable_.lookup(min_max_rzindex_out[0]);
//                 int lutval_inner_from_out_maxr = testOuterInnerTable_.lookup(min_max_rzindex_out[1]);
//                 unsigned int lutwidth_out = settings_.lutwidthtabextended(0, iSeed_);
//                 const auto& lutshift_out = innerPairTable_.nbits();
//                 std::vector<int> inner_rz_outer_min = extract_rz_from_lut(lutval_inner_from_out_minr, lutwidth_out, lutshift_out, nbitsrzbin_in);
//                 std::vector<int> inner_rz_outer_max = extract_rz_from_lut(lutval_inner_from_out_maxr, lutwidth_out, lutshift_out, nbitsrzbin_in);
//                 
//                 std::cout << "NOW FOR INNER FROM OUTER: "  << std::endl ;
//                 std::vector<int> updated_inner_ranges = choose_rzbins_range(inner_rz_outer_min, 
//                                                                             inner_rz_outer_max, 
//                                                                             start_in, 
//                                                                             rzbinfirst_in,
//                                                                             rzdiffmax_in,
//                                                                             next_in
//                                                                             );
                    

              }
              
              trpdata.projbin_out_.emplace_back(tuple<int, int, int, int>(ibin_out - start_out, outmem, nstubs_out, out_phi_region)); 
            } else {
              mask = "0" + mask;
            }
          } // end loop outmem
        }
                
        int max_nstubs = 0;        
        unsigned int max_nbins  = 0;        
        int count_in_mem = 0;
        for (int ibin_in = start_in; ibin_in <= last_in; ibin_in++) {
          for (unsigned int inmem = 0; inmem < innervmstubs_.size(); inmem++) {
            // for each memory, check if its region is compatible 
            unsigned int in_phi_region = (innervmstubs_[inmem]->phibin() - 1) - (innervmstubs_[inmem]->getName()[11] - 'A') * 8;
            if (iSeed_ == Seed::L2L3L4) {
              char mem_reg_in  = innervmstubs_[inmem]->getName()[11];
              char mem_reg_mid = middleallstubs_[midmem]->getName()[8];

              int diff_reg = mem_reg_in - mem_reg_mid;
              if (diff_reg == 1) {
                  in_phi_region += 8;
              } 
              else if (diff_reg == -1) {
                  in_phi_region += 16;
              } 
              if (std::abs(diff_reg) < 2){
                if (usereg_in != -1 && !(usereg_in & (1 << in_phi_region))) {
                  continue;
                }  
              }  
            }
            int nstubs_in = innervmstubs_[inmem]->nVMStubsBinned(ibin_in);
            if (iSeed_ == 8) count_nstubs_in_sara += nstubs_in;
            if (nstubs_in > 0){
              trpdata.projbin_in_.emplace_back(tuple<int, int, int, int>(ibin_in - start_in, inmem, nstubs_in, in_phi_region));
              
              count_in_mem += 1;
              if (nstubs_in > max_nstubs)
                max_nstubs = nstubs_in;
              if (innervmstubs_[inmem]->nBin() > max_nbins)
                max_nbins = innervmstubs_[inmem]->nBin(); 
            }    
          }
        }

        size_t total_inn_stubs = innervmstubs_.size() * max_nbins * max_nstubs;
        trpdata.bad_inner_stub_.assign(total_inn_stubs, false);
        trpdata.max_nbins_ = max_nbins;
        trpdata.max_nstubs_ = max_nstubs;

//         if (iSeed_ == 8) std::cout << iTC_ << "," << istub << "," << count_nstubs_out_sara << "," << count_nstubs_in_sara << std::endl;

        if (!trpdata.projbin_in_.empty() && !trpdata.projbin_out_.empty()) {
          goodtrpdata = true;
//           if (iSeed_ == 8) std::cout << "done for middlestub " << istub 
//                     << ": have n in/out mem: " << trpdata.projbin_in_.size() 
//                     << "/" << trpdata.projbin_out_.size()<< std::endl;
        }
      } // end lutval != -1

      istub++;
      if (istub >= middleallstubs_[midmem]->nStubs()) {
        istub = 0;
        midmem++;
      }

//       if (iSeed_ == 8 && istep == (108 - 1)) 
//         std::cout << "will truncate block3" << std::endl;
    
    } else if (do_next && (!trpbuffernearfull) && midmem < midmemend && istub == 0){
      midmem++;
//       if (iSeed_ == 8) {
//         std::cout << iTC_ << "," << istep << ",increasing midmem to " << midmem;
//         if (midmem == midmemend) std::cout << "next is doing nothing" << std::endl;
//         else std::cout << std::endl;
//       }  
    }  
//     else {
//       if (iSeed_ == 8) std::cout << "block3 done" << std::endl;
//     }  

//     if (iSeed_ == 8 && goodtrpdata) std::cout << "block3, saving into buffer " << std::endl;
    goodtrpdata___ = goodtrpdata;
    trpdata___ = trpdata;

//    original version
//     goodtrpdata___ = goodtrpdata__;
//     goodtrpdata__ = goodtrpdata;
// 
//     trpdata___ = trpdata__;
//     trpdata__ = trpdata;

    //
    // stop looping over istep if done
    //

    bool done = true;

    if (midmem < midmemend || (!trpdatabuffer.empty())) {
      done = false;
    }

    for (auto& trpunit : trpunits_) {
      if (!(trpunit.idle() && trpunit.empty()))
        done = false;
    }

    if (done) {
      donecount++;
    }

    //FIXME This should be done cleaner... Not too hard, but need to check fully the TEBuffer and TEUnit buffer.
    if (donecount > 4) {
      break;
    }
  } // end istep

//   std::cout << iSeed_ << "," << iTC_ << "," << countall << "," << countsel  << std::endl;

  if (settings_.writeMonitorData("TPD")) {
    globals_->ofstream("trackletprocessordisplaced.txt")
        << getName() << " " << countall << " " << countsel << std::endl;
  }
}


double TrackletProcessorDisplaced::calc_phi_tmp(double r, double rho, double d0){
    
    double phi1 = -r/2/rho + d0/r + d0*d0/2/r/rho -2*d0*r/4/rho/rho + 1/6*pow(-r/2/rho + d0/r,3);
    return phi1;
}
double TrackletProcessorDisplaced::calc_deltaPhi(double rinv, double d0, double r1, double r2){
  
  if (r1 <= 0 || r2 <= 0 || r1 >= r2)
    return std::numeric_limits<double>::quiet_NaN();
    
  double rho = 1.0 / rinv;
  double phi1 = calc_phi_tmp(r1, rho, d0);
  double phi2 = calc_phi_tmp(r2, rho, d0);
  double delta_phi = phi2 - phi1;
    
  return delta_phi;
}

std::vector<int> TrackletProcessorDisplaced::find_rzbin_from_zbin(int zbin_out_index, unsigned int zbins_out){
    float zbin_width = 2*settings_.zlength() / (settings_.NLONGVMBINS() * settings_.NLONGVMBINS());
    float zvalue_out = -settings_.zlength() + zbin_width * zbin_out_index;
    float zbin_width_second_lut = 2*settings_.zlength()/zbins_out;
    int indexz_out = (zvalue_out + settings_.zlength()) / zbin_width_second_lut;
    int indexrz_out_minr = indexz_out*16 + 0;
    int indexrz_out_maxr = indexz_out*16 + 15;
    std::vector<int> min_max_rzindex = {indexrz_out_minr, indexrz_out_maxr};

//     std::cout << "  [TPD:outer] indexz:" 
//               << indexz_out << "  indexrz:" 
//               << indexrz_out_minr << " or  " 
//               << indexrz_out_maxr << "   zvalue=" 
//               << zvalue_out << "  " 
// //               << zbin_width << "  "    
// // //               << start_out << "  " 
// // //               << rzbinfirst_out << "  " 
// // //                       << zbin_width_second_lut << " " 
// // //                       << lutval_inner_from_out << " " 
// // //                       << std::bitset<16>(lutval_inner_from_out) << " " 
//               << std::endl;

    return min_max_rzindex;
}

std::vector<int> TrackletProcessorDisplaced::extract_rz_from_lut( int lutval, 
                                                                  unsigned int lutwidth_out, 
                                                                  const unsigned int lutshift_out,
                                                                  int nbitsrzbin_in
                                                                  ){
    FPGAWord lookupbits(lutval, lutwidth_out, true, __LINE__, __FILE__);
    int rzbinfirst_in = lookupbits.bits(lutshift_out, NFINERZBITS);                      // NFINERZBITS = 3, N_RZBITS = 3
    int start_in      = lookupbits.bits(lutshift_out + NFINERZBITS + 1, nbitsrzbin_in);  // first rz bin projection
    int rzdiffmax_in  = lookupbits.bits(lutshift_out + NFINERZBITS + 1 + nbitsrzbin_in, NFINERZBITS);
    int next_in       = lookupbits.bits(lutshift_out + NFINERZBITS, 1);

    std::vector<int> return_values = {start_in, rzbinfirst_in, rzdiffmax_in, next_in};
    float zbin_width = 2*settings_.zlength() / (settings_.NLONGVMBINS() * settings_.NLONGVMBINS());
//     std::cout << "-> z value: " << -settings_.zlength() + zbin_width *( start_in * (8) + rzbinfirst_in) << std::endl;;
    return return_values;                                                            
}

std::vector<int> TrackletProcessorDisplaced::choose_rzbins_range( std::vector<int> inner_rz_outer_min, 
                                                                  std::vector<int> inner_rz_outer_max, 
                                                                  int start_in, 
                                                                  int rzbinfirst_in,
                                                                  int rzdiffmax_in,
                                                                  int next_in
                                                                ){

//    int zbin_in_index_from_middle = start_in * (8) + rzbinfirst_in;
//    float zvalue_in_from_middle = -settings_.zlength() + zbin_width * zbin_in_index_from_middle;
//    std::cout << "from middle \t" << start_in << " " << rzbinfirst_in << " " << rzdiffmax_in << " " << next_in << "  / bin index/value" << zbin_in_index_from_middle << " " 
//              << zvalue_in_from_middle << " "  << std::endl;
//    std::cout << "from outer min r\t"   << start_in_from_out_minr << " "   << rzbinfirst_in_from_out_minr << " " << rzdiffmax_in_from_out_minr << " " << std::endl;
//    std::cout << "from outer max r\t"   << start_in_from_out_maxr << " "   << rzbinfirst_in_from_out_maxr << " " << rzdiffmax_in_from_out_maxr << " " << std::endl;

    int final_start_in, final_rzbinfirst_in, final_inner_start, final_next_in;
    int full_inner_index_from_outer;
    int full_inner_index_from_middle = start_in * (8) + rzbinfirst_in;
    int full_inner_index_from_outer_minr = inner_rz_outer_min[0] * (8) + inner_rz_outer_min[1];
    int full_inner_index_from_outer_maxr = inner_rz_outer_max[0] * (8) + inner_rz_outer_max[1];
    
    // here take the wider range
    if (full_inner_index_from_outer_maxr > full_inner_index_from_outer_minr && full_inner_index_from_outer_minr >=0){
      full_inner_index_from_outer = full_inner_index_from_outer_minr;
      final_start_in = inner_rz_outer_min[0];
      final_rzbinfirst_in = inner_rz_outer_min[1];
    } else if (full_inner_index_from_outer_maxr >= 0){ 
      full_inner_index_from_outer = full_inner_index_from_outer_maxr;
      final_start_in = inner_rz_outer_max[0];
      final_rzbinfirst_in = inner_rz_outer_max[1];
    } else {      
      return {start_in, rzbinfirst_in, rzdiffmax_in, next_in};
    }
    // here take the smaller range
    if (full_inner_index_from_middle > full_inner_index_from_outer){
      final_inner_start = full_inner_index_from_middle;
      final_start_in = start_in;
      final_rzbinfirst_in = rzbinfirst_in;
    } else { 
      final_inner_start = full_inner_index_from_outer;
    }
    
    int end_inner_index_from_middle = full_inner_index_from_middle + rzdiffmax_in;
    int end_inner_index_from_outer_minr = full_inner_index_from_outer_minr + inner_rz_outer_min[2];
    int end_inner_index_from_outer_maxr = full_inner_index_from_outer_maxr + inner_rz_outer_max[2];
    int end_inner_index_from_outer = end_inner_index_from_outer_minr > end_inner_index_from_outer_maxr ? 
                                      end_inner_index_from_outer_minr : end_inner_index_from_outer_maxr;
    
    int final_inner_end =  (end_inner_index_from_middle < end_inner_index_from_outer) ? 
                            end_inner_index_from_middle : end_inner_index_from_outer;
    if (final_inner_end <  final_inner_start)
      std::cout << "WARNING!!!!" << " final_inner_end = " << final_inner_end << " while final_inner_start = " << final_inner_start  << std::endl;            
    int final_rzdiffmax_in = final_inner_end - final_inner_start;           
    int final_next = final_rzdiffmax_in > 7 ? 1 : 0;                
    if (final_rzdiffmax_in > 7) 
      final_rzdiffmax_in = 7;
    
//   if (final_rzdiffmax_in < rzdiffmax_in){
//     std::cout << "reducing range: " << rzdiffmax_in << "  now  " <<  final_rzdiffmax_in << std::endl;            
//     std::cout << "reducing final_start_in: " << start_in << "  now  " <<  final_start_in << std::endl;            
//     std::cout << "reducing rzbinfirst_in: " << rzbinfirst_in << "  now  " <<  final_rzbinfirst_in << std::endl;            
//   }
//             
//   std::cout << "total \t" 
//             << full_inner_index_from_middle << " - " 
//             << end_inner_index_from_middle << "    " 
//             << full_inner_index_from_outer_minr << " -  " 
//             << end_inner_index_from_outer_minr << "  " 
//             << full_inner_index_from_outer_maxr << "    " 
//             << end_inner_index_from_outer_maxr << " ---> " 
//             << final_inner_start   << " - " 
//             << final_inner_end   
//             << std::endl;
//     if (final_start_in < 0 || final_rzbinfirst_in < 0 || final_rzdiffmax_in < 0 || final_next < 0)
//       std::cout << "in function \t" 
//             << full_inner_index_from_middle << " - " 
//             << end_inner_index_from_middle << "    " 
//             << full_inner_index_from_outer_minr << " - " 
//             << end_inner_index_from_outer_minr << "  " 
//             << full_inner_index_from_outer_maxr << " - " 
//             << end_inner_index_from_outer_maxr << " ---> " 
//                 << final_start_in << " - " 
//                 << final_rzbinfirst_in << " - " 
//                 << final_rzdiffmax_in << " -  " 
//                 << final_next << "   --->  ---> " 
//                 << final_start_in*8 + final_rzbinfirst_in << " - " 
//                 << final_start_in*8 + final_rzbinfirst_in + final_rzdiffmax_in<< " - " 
//                 << std::endl;
     
    std::vector<int> return_values = {final_start_in, final_rzbinfirst_in, final_rzdiffmax_in, final_next};
    return return_values;           
}


std::vector<float> TrackletProcessorDisplaced::linspace(float A, float B, int N) {
    std::vector<float> v(N);
    if (N == 1) return {A};
    float step = (B - A) / (N - 1);
    float val = A;
    std::generate(v.begin(), v.end(), [&]{ float tmp = val; val += step; return tmp; });
    return v;
}

double TrackletProcessorDisplaced::compute_dphimax (std::vector<float> d0_vals, double rinv, double r1, double r2) {
  std::vector<double> dphi_vals;
  dphi_vals.reserve(2 * d0_vals.size());

  for (double d0 : d0_vals) {
    dphi_vals.push_back(calc_deltaPhi( rinv, d0, r1, r2));
    dphi_vals.push_back(calc_deltaPhi(-rinv, d0, r1, r2));
  }

  auto max_it = std::max_element(
      dphi_vals.begin(),
      dphi_vals.end(),
      [](double a, double b) { return std::abs(a) < std::abs(b); });

  return std::abs(*max_it);
};

int TrackletProcessorDisplaced::compute_nfinephibins(int layerdisk, int iSeed_, double dphimax){
   
  int nfinephibins =
      settings_.nallstubs(layerdisk) * settings_.nvmte(1, iSeed_) * (1 << settings_.nfinephi(1, iSeed_));
  double dfinephi = settings_.dphisectorHG() / nfinephibins;
  // next line: passing layerdisk2, but nbitsallstubs is the same for all layers except L1, so it does not really matter for displaced seeds
  // it also means that nbitsfinephi_ can be used for outer and inner
  int nbins_out = 2.0 * (dphimax / dfinephi + 1.0);
  // find the number of bits needed to represent nbins
  return (log(nbins_out) / log(2.0) + 1);
  
}