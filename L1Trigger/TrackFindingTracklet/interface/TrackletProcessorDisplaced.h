// TrackletProcessorDisplaced: This class performs the tasks of the TrackletEngineDisplaced+TripletEngine+TrackletCalculatorDisplaced.
#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletProcessorDisplaced_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletProcessorDisplaced_h

#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorDisplaced.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletLUT.h"
#include "L1Trigger/TrackFindingTracklet/interface/CircularBuffer.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletParametersMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TripletEngineUnit.h"

#include <vector>
#include <tuple>
#include <map>

namespace trklet {

  class Settings;
  class Globals;
  class MemoryBase;
  class AllStubsMemory;
  class AllInnerStubsMemory;
  class VMStubsTEMemory;
  class StubPairsMemory;

  class TrackletProcessorDisplaced : public TrackletCalculatorDisplaced {
  public:
    TrackletProcessorDisplaced(std::string name, Settings const& settings, Globals* globals);

    ~TrackletProcessorDisplaced() override = default;

    void addOutputProjection(TrackletProjectionsMemory*& outputProj, MemoryBase* memory);

    void addOutput(MemoryBase* memory, std::string output) override;

    void addInput(MemoryBase* memory, std::string input) override;

    void execute(unsigned int iSector, double phimin, double phimax, std::vector<L1StubTriplet>&, std::vector<L1StubTriplet>&);

    double calc_phi_tmp(double r, double rho, double d0);
    double calc_deltaPhi(double rho, double d0, double r1, double r2);
    double compute_dphimax(std::vector<float> d0_vals, double rinv, double r1, double r2);
    int compute_nfinephibins(int layerdisk, int iSeed_, double dphimax);

    std::vector<float> linspace(float A, float B, int N);
    std::vector<int> find_rzbin_from_zbin(int zbin_out_index, unsigned int zbins_out);

    std::vector<int> choose_rzbins_range(std::vector<int>, 
                                         std::vector<int>, 
                                         int start_in, 
                                         int rzbinfirst_in,
                                         int rzdiffmax_in,
                                         int next_in
                                         );
    std::vector<int> extract_rz_from_lut(int lutval_inner_from_out_minr, 
                                        unsigned int lutwidth_out, 
                                        const unsigned int lutshift_out,
                                        int nbitsrzbin_in);
    

  private:
    int iTC_;
    int iAllStub_;
    unsigned int maxStep_;

    std::tuple<CircularBuffer<TrpEData>, unsigned int, unsigned int, unsigned int, unsigned int> trpbuffer_;
    std::vector<TripletEngineUnit> trpunits_;

    unsigned int layerdisk1_;
    unsigned int layerdisk2_;
    unsigned int layerdisk3_;

    int firstphibits_;
    int secondphibits_;
    int thirdphibits_;

    int nbitsfinephi_; 
    int nbitsfinephiouterdiff_; 
    int nbitsfinephiinnerdiff_; 
    int nbitsfinephipairdiff_; 

    int nbitszfinebintable_;
    int nbitsrfinebintable_;

    TrackletLUT innerTable_;       //projection to next layer/disk
    TrackletLUT innerThirdTable_;  //projection to third disk/layer
    TrackletLUT outerPairTable_;   //projection to outer disk/layer
    TrackletLUT innerPairTable_;   //projection to outer disk/layer
    TrackletLUT testOuterInnerTable_;  //projection from outer to inner

    TrackletLUT useOuterRegiontable_;   // phi LUT
    TrackletLUT useInnerRegiontable_;   // phi LUT

    TrackletLUT pttablemiddle_;
    TrackletLUT pttableouter_;
    TrackletLUT pttablemiddlein_;
    TrackletLUT pttableinner_;

    TrackletLUT pttablemiddle_region_out_;
    TrackletLUT pttablemiddle_region_in_;

    std::vector<VMStubsTEMemory*> innervmstubs_;
    std::vector<VMStubsTEMemory*> outervmstubs_;
  };

};  // namespace trklet
#endif
