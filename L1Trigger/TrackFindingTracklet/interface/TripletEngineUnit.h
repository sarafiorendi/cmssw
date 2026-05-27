#ifndef L1Trigger_TrackFindingTracklet_interface_TripletEngineUnit_h
#define L1Trigger_TrackFindingTracklet_interface_TripletEngineUnit_h

#include "L1Trigger/TrackFindingTracklet/interface/VMStubsTEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/CircularBuffer.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletLUT.h"
#include "DataFormats/L1TrackTrigger/interface/L1StubTriplet.h"

#include <cassert>
#include <vector>

namespace trklet {

  class Settings;
  class Stub;
  class L1TStub;

  struct TrpEData {
    const Stub* stub_;
    int start_out_;
    int start_in_;
    int rzbinfirst_out_;
    int rzdiffmax_out_;
    int rzbinfirst_in_;
    int rzdiffmax_in_;
    int middlefinephi_;
    FPGAWord middlebend_;
    unsigned int max_nbins_;
    unsigned int max_nstubs_;
    std::vector<uint8_t> bad_inner_stub_;
    std::vector<uint8_t> rzbinfirst_in_pair_;
    std::vector<uint8_t> rzstart_in_pair_;
    std::vector<uint8_t> rzdiffmax_in_pair_;
    std::vector<uint8_t> rznext_in_pair_;
    std::vector<std::tuple<int, int, int, int> > projbin_out_;  // next z/r bin; outer stub mem; nstub
    std::vector<std::tuple<int, int, int, int> > projbin_in_;   // next z/r bin; inner stub mem; nstub
  };

  class TripletEngineUnit {
  public:
    TripletEngineUnit(const Settings* const settings,
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
                      const TrackletLUT* pairphiwindow,
                      std::vector<VMStubsTEMemory*> innervmstubs,
                      std::vector<VMStubsTEMemory*> outervmstubs,
                      const std::vector<unsigned int>& innervmrow
                      );

    ~TripletEngineUnit() = default;

    void init(const TrpEData& trpdata);

    bool getGoodTriplet() { return goodtriplet__; }

    bool empty() const { return candtriplets_.empty(); }

    const std::tuple<const Stub*, const Stub*, const Stub*>& read() { return candtriplets_.read(); }

    const std::tuple<const Stub*, const Stub*, const Stub*>& peek() const { return candtriplets_.peek(); }

    bool idle() const { return idle_; }

    void setNearFull() { nearfull_ = candtriplets_.nearfull(); }

    void reset();

    void advanceIndices();
    void advanceInnerBin();
    void advanceOuterStub();

    bool innerEntryCompatible(unsigned int, int, int, int, int) const;
    size_t badInnerIndex() const ;
    bool nextGoodInner() ;

    void step(std::vector<L1StubTriplet>&, unsigned int, int, int);    

    // for debugging
    unsigned int nTested() const { return n_triplets_tested_; }
    unsigned int nSkippedBadInner() const { return n_skipped_bad_; }
    unsigned int nRejOuter() const { return n_rej_outer_; }
    unsigned int nRejOuterFromBend() const { return n_rej_outer_pair_bend_; }
    unsigned int nRejPair() const { return n_rej_pair_; }
    unsigned int nRejInner() const { return n_rej_inner_; }
    unsigned int nRejMemInner() const { return n_rej_in_mem_pair_; }
//     unsigned int nSkippedEntries() const { return n_skipped_entries_; }

    const Stub* innerStub() const { return trpdata_.stub_; }

  private:
    std::vector<VMStubsTEMemory*> innervmstubs_;
    std::vector<VMStubsTEMemory*> outervmstubs_;
    TrpEData trpdata_;
    const Settings* settings_;
    unsigned int layerdisk1_;
    unsigned int layerdisk2_;
    unsigned int layerdisk3_;
    unsigned int iSeed_;
    unsigned int iAllStub_;
    unsigned int nbitsfinephi_;
    unsigned int nbitsfinephiouterdiff_;
    unsigned int nbitsfinephiinnerdiff_;
    bool nearfull_;  //initialized at start of each processing step

    //unsigned int memory slot
    unsigned int nmem_out_ = 0;
    unsigned int nmem_in_ = 0;
    unsigned int istub_out_ = 0;
    unsigned int istub_in_ = 0;
    unsigned int next_out_ = 0;
    unsigned int next_in_ = 0;
    unsigned int nstub_out_ = 0;
    unsigned int nstub_in_ = 0;
    unsigned int outmem_ = 0;
    unsigned int inmem_ = 0;
    unsigned int nproj_out_ = 0;
    unsigned int nproj_in_ = 0;
    unsigned int phi_out_ = 0;
    unsigned int phi_in_ = 0;
    
    bool idle_ = true;

    int n_triplets_tested_; 
    int n_skipped_bad_;
    int n_rej_outer_;
    int n_rej_outer_pair_bend_;
    int n_rej_pair_;
    int n_rej_inner_;
    int n_rej_in_mem_pair_;

    const TrackletLUT* pttablemiddle_;
    const TrackletLUT* pttableouter_;
    const TrackletLUT* pttablemiddlein_;
    const TrackletLUT* pttableinner_;
    const TrackletLUT* pairphiwindow_;

    std::vector<int> innervmshift_;  // per inner VM memory: log2(middle fine-phi bins per VM bin)
    
    std::tuple<const Stub*, const Stub*, const Stub*> candtriplet_, candtriplet__;
    bool goodtriplet_ = false;
    bool goodtriplet__ = false;

    //save the candidate matches
    CircularBuffer<std::tuple<const Stub*, const Stub*, const Stub*> > candtriplets_;
  };  // TripletEngineUnit

};  // namespace trklet
#endif
