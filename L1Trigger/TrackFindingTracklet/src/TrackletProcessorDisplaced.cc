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

  // set layer/disk types based on input seed name
  initLayerDisksandISeedDisp(layerdisk1_, layerdisk2_, layerdisk3_, iSeed_);

  // find the maximum deltaPhi (assuming max d0 of XX) expected for stubs in the seed
  double rmiddle = -1.0;
  double router  = -1.0;
  double rinner  = -1.0;

  if (iSeed_ == Seed::L2L3L4 || iSeed_ == Seed::L4L5L6) {
    rmiddle = settings_.rmean(layerdisk1_); // middle
    router  = settings_.rmean(layerdisk2_); // outer
    rinner  = settings_.rmean(layerdisk3_); // inner
  } //else {
  // add missing seeds 10 and 11
  // } 
  // rmin and rmax correspond to the radius of the two layers/disks of the seed 
  // dphimax (in radians) is the max delta phi between the two stubs of the pair to form that seed in order to satisfy the min pt cut (maxrinv)
  const double rinv = settings_.rinvmaxDisplaced(); // assuming charge = 1
  // FIXME: max/min d0 to be configurable
  std::vector<float> d0_vals = linspace(-10,10,101);
  double dphimax_out = compute_dphimax(d0_vals, rinv, rmiddle, router);
  double dphimax_in = compute_dphimax(d0_vals, rinv, rinner, rmiddle);

  // number of fine phi bins in sector for the outer/inner layer
  nbitsfinephiouterdiff_ = compute_nfinephibins(layerdisk2_,iSeed_, dphimax_out);
  nbitsfinephiinnerdiff_ = compute_nfinephibins(layerdisk3_,iSeed_, dphimax_in);

  nbitsfinephi_ = settings_.nbitsallstubs(layerdisk2_) + settings_.nbitsvmte(1, iSeed_) + settings_.nfinephi(1, iSeed_);

  // get projection tables
  unsigned int region = name.back() - 'A';
  innerTable_.initVMRTable(
      layerdisk1_, TrackletLUT::VMRTableType::inner, region, false);  //projection to next layer/disk
  innerThirdTable_.initVMRTable(
      layerdisk1_, TrackletLUT::VMRTableType::innerthird, region, false);  //projection to third layer/disk

  nbitszfinebintable_ = settings_.vmrlutzbits(layerdisk1_);
  nbitsrfinebintable_ = settings_.vmrlutrbits(layerdisk1_);

  // set TC index
  iTC_ = region;
  TCIndex_ = (iSeed_ << settings.nbitsseed()) + iTC_;

//   maxStep_ = settings_.maxStep("TPD");
  maxStep_ = 108;
  if (maxStep_ < 200) 
    std::cout << "WARNING: Truncation ON! maxStep_ = " << maxStep_ << std::endl;
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

    pttablemiddle_region_out_.initTPlut(true, iSeed_, layerdisk1_, layerdisk2_, nbitsfinephiouterdiff_, iTP, false);
    useOuterRegiontable_.initTPregionlut(
      iSeed_, layerdisk1_, layerdisk2_, iAllStub_, nbitsfinephiouterdiff_, nbitsfinephi_, pttablemiddle_region_out_, iTP);
      // iTP is only used in the name of the LUT table, if written out
      // iAllStub is used to define the outerfinephi
    pttablemiddle_region_in_.initTPlut(true, iSeed_, layerdisk1_, layerdisk3_, nbitsfinephiinnerdiff_, iTP, true);
    useInnerRegiontable_.initTPregionlut(
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
    pttablemiddlein_.initTPlut(true, iSeed_, layerdisk1_, layerdisk3_, nbitsfinephiinnerdiff_, iTP, true);
    pttableinner_.initTPlut(false, iSeed_, layerdisk1_, layerdisk3_, nbitsfinephiinnerdiff_, iTP, true);
    return;
  }
  if (input == "secondvmstubin") {
    auto* tmp = dynamic_cast<VMStubsTEMemory*>(memory);
    assert(tmp != nullptr);
    outervmstubs_.push_back(tmp);

    unsigned int iTP = getName()[7] - 'A';
    pttablemiddle_.initTPlut(true, iSeed_, layerdisk1_, layerdisk2_, nbitsfinephiouterdiff_, iTP, false);
    pttableouter_.initTPlut(false, iSeed_, layerdisk1_, layerdisk2_, nbitsfinephiouterdiff_, iTP, false);
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

  // reset the tebuffer
  std::get<0>(trpbuffer_).reset();
  std::get<1>(trpbuffer_) = 0;
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
      }
      count_trpunits_block1++;      
    }

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

      L1StubTriplet myTripletNow = makeL1StubTriplet(innerFPGAStub, middleFPGAStub, outerFPGAStub, iSector, iTC_, the_trpunit_being_read);

      // check if the seed made from the 3 stubs is valid
      bool accept = false;
      if (iSeed_ == Seed::L2L3L4 || iSeed_ == Seed::L4L5L6)
        accept = LLLSeeding(innerFPGAStub, innerStub, middleFPGAStub, middleStub, outerFPGAStub, outerStub);
      else if (iSeed_ == Seed::L2L3D1)
        accept = LLDSeeding(innerFPGAStub, innerStub, middleFPGAStub, middleStub, outerFPGAStub, outerStub);
      else if (iSeed_ == Seed::D1D2L2)
        accept = DDLSeeding(innerFPGAStub, innerStub, middleFPGAStub, middleStub, outerFPGAStub, outerStub);

      if (accept) {
        acceptedtriplets_.push_back(myTripletNow);
        countsel++;
      }
      
      if (trackletpars_->nTracklets() >= settings_.ntrackletmax()) {
        edm::LogVerbatim("Tracklet") << "Will break on number of tracklets in " << getName();
        assert(0);
        break;
      }

      if (settings_.debugTracklet()) {
        edm::LogVerbatim("Tracklet") << "TrackletProcessor execute done";
      }
    }

    //
    // The second block fills the trpunit if data in buffer and process TripletEngineUnit step
    //
    //

    bool notemptytrpbuffer = !trpdatabuffer.empty();
    int count_trpunits = 0;    
    for (auto& trpunit : trpunits_) {
      if (trpunit.idle() && notemptytrpbuffer) {  // only fill one idle unit every step
        trpunit.init(std::get<0>(trpbuffer_).read());
        notemptytrpbuffer = false;  //prevent initializing another triplet engine unit
      }
      trpunit.step(foundtriplets_, iSector, iTC_, count_trpunits);
      count_trpunits++;            
    }

    //
    // The third block here checks if we have input stubs to process
    //
    //

    if (goodtrpdata___)
      trpdatabuffer.store(trpdata___);
    goodtrpdata = false;

    unsigned int& istub = std::get<1>(trpbuffer_);
    unsigned int& midmem = std::get<2>(trpbuffer_);
    unsigned int midmemend = std::get<4>(trpbuffer_);

    if ((!trpbuffernearfull) && midmem < midmemend && istub < middleallstubs_[midmem]->nStubs()) {

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
      int lutval2 = innerThirdTable_.lookup((indexz << nbitsrfinebintable_) + indexr);
      // TO BE CHECKED      
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
        int last_in = start_in + next_in;  // last large rz-bin projection

        // find the phi region where outer / inner stubs are expected, based on middle stub phi
        // and on the LUTs
        FPGAWord phicorr = stub->phicorr();
        // nbitsallstubs is the same for all layers except L1
        // the other two pars only depend on the seed, not on the layer 
        int middlefinephi = phicorr.bits(phicorr.nbits() - nbitsfinephi_, nbitsfinephi_);
        FPGAWord middlebend = stub->bend();

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
        trpdata.rzbinfirst_out_ = rzbinfirst_out;
        trpdata.rzdiffmax_out_ = rzdiffmax_out;
        trpdata.rzbinfirst_in_ = rzbinfirst_in;
        trpdata.rzdiffmax_in_ = rzdiffmax_in;
        trpdata.start_out_ = start_out;
        trpdata.start_in_ = start_in;
        trpdata.middlefinephi_ = middlefinephi;
        trpdata.middlebend_ = middlebend;

        // fill projection bins info for single engine unit
        trpdata.projbin_out_.clear();
        trpdata.projbin_in_.clear();
        trpdata.bad_inner_stub_.clear();

        for (int ibin_out = start_out; ibin_out <= last_out; ibin_out++) {
          // outervmstubs_.size() is the number of outer memories as from the wiring file, e.g.
          // instance A has 9 outer memories, instance B -> 6, instance C -> 6, instance D -> 10...
          // each outer memory covers 1/settings_.nvmte(1, iSeed_) of one of the 4 large phi regions of the outer layer (1/8*1/4 of L4 for seed 8 for example)
          for (unsigned int outmem = 0; outmem < outervmstubs_.size(); outmem++) {
            // for each memory, check if its region is compatible 
            unsigned int out_phi_region = (outervmstubs_[outmem]->phibin() - 1) - (outervmstubs_[outmem]->getName()[11] - 'A') * 8;
            if (iSeed_ == Seed::L2L3L4) {
              char mem_reg_out  = outervmstubs_[outmem]->getName()[11];
              char mem_reg_mid = middleallstubs_[midmem]->getName()[8];

              int diff_reg = mem_reg_out - mem_reg_mid;
              if (diff_reg == 1) 
                  out_phi_region += 8; // FIXME: not use hardcoded number
              else if (diff_reg == -1) 
                  out_phi_region += 16;
              if (std::abs(diff_reg) < 2){
                if (usereg_out != -1 && !(usereg_out & (1 << out_phi_region))) {
                  continue;
                }  
              }  
            } 
            // check if the memory is not empty
            unsigned int nstubs_out = outervmstubs_[outmem]->nVMStubsBinned(ibin_out);
            if (nstubs_out > 0)
              trpdata.projbin_out_.emplace_back(tuple<int, int, int, int>(ibin_out - start_out, outmem, nstubs_out, out_phi_region));
          }
        }


        unsigned int max_nstubs = 0;        
        unsigned int max_nbins  = 0;        
        for (int ibin_in = start_in; ibin_in <= last_in; ibin_in++) {
          for (unsigned int inmem = 0; inmem < innervmstubs_.size(); inmem++) {
            // for each memory, check if its region is compatible 
            unsigned int in_phi_region = (innervmstubs_[inmem]->phibin() - 1) - (innervmstubs_[inmem]->getName()[11] - 'A') * 8;
            if (iSeed_ == Seed::L2L3L4) {
              char mem_reg_in  = innervmstubs_[inmem]->getName()[11];
              char mem_reg_mid = middleallstubs_[midmem]->getName()[8];

              int diff_reg = mem_reg_in - mem_reg_mid;
              if (diff_reg == 1) 
                  in_phi_region += 8;
              else if (diff_reg == -1) 
                  in_phi_region += 16;
              if (std::abs(diff_reg) < 2){
                if (usereg_in != -1 && !(usereg_in & (1 << in_phi_region))) {
                  continue;
                }  
              }  
            }
            unsigned int nstubs_in = innervmstubs_[inmem]->nVMStubsBinned(ibin_in);
            if (nstubs_in > 0) {
              trpdata.projbin_in_.emplace_back(tuple<int, int, int, int>(ibin_in - start_in, inmem, nstubs_in, in_phi_region));
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

        if (!trpdata.projbin_in_.empty() && !trpdata.projbin_out_.empty()) {
          goodtrpdata = true;
        }
      }

      istub++;
      if (istub >= middleallstubs_[midmem]->nStubs()) {
        istub = 0;
        midmem++;
      }

    } else if ((!trpbuffernearfull) && midmem < midmemend && istub == 0)
      midmem++;

    goodtrpdata___ = goodtrpdata__;
    goodtrpdata__ = goodtrpdata;

    trpdata___ = trpdata__;
    trpdata__ = trpdata;

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
  }

  if (settings_.writeMonitorData("TPD")) {
    globals_->ofstream("trackletprocessordisplaced.txt")
        << getName() << " " << countall << " " << countsel << std::endl;
  }
}

double TrackletProcessorDisplaced::compute_dphimax (std::vector<float> d0_vals, const double rinv, double r1, double r2) {
  std::vector<double> dphi_vals;
  dphi_vals.reserve(2 * d0_vals.size());

  for (double d0 : d0_vals) {
    dphi_vals.push_back(compute_deltaPhi( rinv, d0, r1, r2));
    dphi_vals.push_back(compute_deltaPhi(-rinv, d0, r1, r2));
  }

  if (dphi_vals.empty())
    return 9999.0;
  
  auto max_it = std::max_element(
      dphi_vals.begin(),
      dphi_vals.end(),
      [](double a, double b) { return std::abs(a) < std::abs(b); });

  return std::abs(*max_it);
};

double TrackletProcessorDisplaced::compute_phi(double r, double rho, double d0){
    return (-r/2./rho + d0/r + d0*d0/2./r/rho -2*d0*r/4./rho/rho + 1.0/6.0*pow(-r/2/rho + d0/r,3));
}

double TrackletProcessorDisplaced::compute_deltaPhi(const double rinv, double d0, double r1, double r2){
  
  if (r1 <= 0 || r2 <= 0 || r1 >= r2 || std::abs(rinv) < 1e-12)
    return std::numeric_limits<double>::quiet_NaN();
    
  double rho = 1.0 / rinv;
  double phi1 = compute_phi(r1, rho, d0);
  double phi2 = compute_phi(r2, rho, d0);
  double delta_phi = phi2 - phi1;
  return delta_phi;
}

int TrackletProcessorDisplaced::compute_nfinephibins(int layerdisk, int iSeed_, double dphimax){
   
  int nfinephibins =
      settings_.nallstubs(layerdisk) * settings_.nvmte(1, iSeed_) * (1 << settings_.nfinephi(1, iSeed_));
  double dfinephi = settings_.dphisectorHG() / nfinephibins;
  // next line: passing layerdisk2, but nbitsallstubs is the same for all layers except L1, so it does not really matter for displaced seeds
  // it also means that nbitsfinephi_ can be used for outer and inner
  int nbins_out = 2.0 * (dphimax / dfinephi + 1.0);
  // find the number of bits needed to represent nbins
//   int bits = 0;
//   while ((1 << bits) < nbins_out) bits++;
  int bits = log(nbins_out) / log(2.0) + 1;
  
  return bits;
  
}
std::vector<float> TrackletProcessorDisplaced::linspace(float A, float B, int N) {
    std::vector<float> v(N);
    if (N == 1) return {A};
    float step = (B - A) / (N - 1);
    float val = A;
    std::generate(v.begin(), v.end(), [&]{ float tmp = val; val += step; return tmp; });
    return v;
}
