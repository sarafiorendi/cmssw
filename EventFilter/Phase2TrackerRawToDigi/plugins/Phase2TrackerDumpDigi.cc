#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
//#include "DataFormats/FEDRawData/src/fed_header.h"
//#include "DataFormats/FEDRawData/src/fed_trailer.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
// #include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerDigiToRaw.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"


#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"

#include <sstream>
// #include <iomanip>
// #include <ext/algorithm>

using namespace std;
using namespace Phase2Tracker;


class Phase2TrackerDumpDigi : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  /// constructor
  Phase2TrackerDumpDigi(const edm::ParameterSet& pset);
  /// default constructor
  ~Phase2TrackerDumpDigi() override = default;
  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {};
  void analyze(const edm::Event&, const edm::EventSetup&);

private:
  const edm::ESGetToken<Phase2TrackerCabling, Phase2TrackerCablingRcd> ph2CablingESToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::EDGetTokenT<Phase2TrackerCluster1DCollectionNew> token_;
  const Phase2TrackerCabling* cabling_ = nullptr;
  const TrackerTopology* tTopo_ = nullptr;
  const TrackerGeometry* tGeom_ = nullptr;
  std::map<int, std::pair<int, int>> stackMap_;
};

Phase2TrackerDumpDigi::Phase2TrackerDumpDigi(const edm::ParameterSet& pset)
    : ph2CablingESToken_(esConsumes<Phase2TrackerCabling, Phase2TrackerCablingRcd, edm::Transition::BeginRun>()),
      geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
      topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
      token_(
          consumes<Phase2TrackerCluster1DCollectionNew>(pset.getParameter<edm::InputTag>("ProductLabel"))) {
}

void Phase2TrackerDumpDigi::beginRun(edm::Run const& run, edm::EventSetup const& es) {
    cabling_ = &es.getData(ph2CablingESToken_);
    tGeom_ = &es.getData(geomToken_);
    tTopo_ = &es.getData(topoToken_);

    // build map of upper and lower for each module
//     for (auto iu = tGeom_->detUnits().begin(); iu != tGeom_->detUnits().end(); ++iu) {
//       unsigned int detId_raw = (*iu)->geographicalId().rawId();
//       DetId detId = DetId(detId_raw);
//       if (detId.det() == DetId::Detector::Tracker) {
//         if (tTopo_->isLower(detId) != 0) {
//           stackMap_[tTopo_->stack(detId)].first = detId;
//         }
//         if (tTopo_->isUpper(detId) != 0) {
//           stackMap_[tTopo_->stack(detId)].second = detId;
//         }
//       }
//     }  // end loop on detunits
  }

void Phase2TrackerDumpDigi::analyze(const edm::Event& event, const edm::EventSetup& es) {
//     std::unique_ptr<FEDRawDataCollection> buffers(new FEDRawDataCollection);
  edm::Handle<Phase2TrackerCluster1DCollectionNew> clusters_handle;
  event.getByToken(token_, clusters_handle);
  for (const auto& DSVItr : *clusters_handle) {
    // Getting the id of detector unit
    uint32_t rawid(DSVItr.detId());
    DetId detId(rawid);
    const GeomDetUnit* geomDetUnit(tGeom_->idToDetUnit(detId));
    if (!geomDetUnit)
      continue;
    for (const auto& clusterItr : DSVItr) { 
    
      // temporary restriction to 2S modules
      TrackerGeometry::ModuleType mType = tGeom_->getDetectorType(detId);
      if (mType != TrackerGeometry::ModuleType::Ph2SS)
        continue;
      
      std::cout << "detId: " << detId.rawId() << std::endl;
      // I could find this number in the DetId sensors list e.g. https://cms-tklayout.web.cern.ch/cms-tklayout/layouts/repository-git-dev/OT616_200_IT404/index.html
      MeasurementPoint mpCluster(clusterItr.center(), clusterItr.column() + 0.5);
      Local3DPoint localPosCluster = geomDetUnit->topology().localPosition(mpCluster);
      Global3DPoint globalPosCluster = geomDetUnit->surface().toGlobal(localPosCluster);
      std::cout << "\t cluster r position: " << globalPosCluster.perp()  << std::endl; 
      std::cout << "\t cluster global z position: " << globalPosCluster.z()  << std::endl; 
    }  
  }      

}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Phase2TrackerDumpDigi);
