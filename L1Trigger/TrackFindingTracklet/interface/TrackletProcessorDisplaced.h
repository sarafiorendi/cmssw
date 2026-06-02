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

    std::vector<float> linspace(float A, float B, int N);

    double compute_phi(double r, double rho, double d0);

    double compute_deltaPhi(double rho, double d0, double r1, double r2);

    double compute_dphimax(std::vector<float> d0_vals, double rinv, double r1, double r2);

    int compute_nfinephibins(int layerdisk, int iSeed_, double dphimax);

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

    int nbitszfinebintable_;
    int nbitsrfinebintable_;

    TrackletLUT innerTable_;       //projection to next layer/disk
    TrackletLUT innerThirdTable_;  //projection to third disk/layer

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

    std::vector<AllStubsMemory*> innerallstubs_;
    std::vector<AllStubsMemory*> middleallstubs_;
    std::vector<AllStubsMemory*> outerallstubs_;
  };

};  // namespace trklet
#endif
