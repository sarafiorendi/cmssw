// system includes
#include <utility>
#include <vector>

// user includes
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDBuffer.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDChannel.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDRawChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDZSChannelUnpacker.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/utils.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "CondFormats/DataRecord/interface/Phase2TrackerCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/Phase2TrackerCabling.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#define LOGPRINT edm::LogPrint("Phase2TrackerFEDTestAnalyzer")

/**
   @class Phase2TrackerFEDTestAnalyzer 
   @brief Analyzes contents of FED_test_ collection
*/

class Phase2TrackerFEDTestAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  typedef std::pair<uint16_t, uint16_t> Fed;
  typedef std::vector<Fed> Feds;
  typedef std::vector<uint16_t> Channels;
  typedef std::map<uint16_t, Channels> ChannelsMap;

  Phase2TrackerFEDTestAnalyzer(const edm::ParameterSet&);
  ~Phase2TrackerFEDTestAnalyzer();

  void beginJob();
  void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {};
  void endJob();

private:
  std::map<int, std::pair<int, int>> stackMap_;
  const edm::ESGetToken<Phase2TrackerCabling, Phase2TrackerCablingRcd> ph2CablingESToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::EDGetTokenT<FEDRawDataCollection> token_;
  const TrackerGeometry* tkGeom_ = nullptr;
  const TrackerTopology* tTopo_ = nullptr;
  const Phase2TrackerCabling* cabling_ = nullptr;
};

using namespace Phase2Tracker;
using namespace std;

// -----------------------------------------------------------------------------
//
Phase2TrackerFEDTestAnalyzer::Phase2TrackerFEDTestAnalyzer(const edm::ParameterSet& pset)
  : ph2CablingESToken_(esConsumes<Phase2TrackerCabling, Phase2TrackerCablingRcd, edm::Transition::BeginRun>()),
    geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
    topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
    token_(consumes<FEDRawDataCollection>(pset.getParameter<edm::InputTag>("ProductLabel"))) {
  LogDebug("Phase2TrackerFEDTestAnalyzer") << "[Phase2TrackerFEDTestAnalyzer::" << __func__ << "]"
                                           << "Constructing object...";
}

// -----------------------------------------------------------------------------
//
Phase2TrackerFEDTestAnalyzer::~Phase2TrackerFEDTestAnalyzer() {
  LogDebug("Phase2TrackerFEDTestAnalyzer") << "[Phase2TrackerFEDTestAnalyzer::" << __func__ << "]"
                                           << " Destructing object...";
}

// -----------------------------------------------------------------------------
//
void Phase2TrackerFEDTestAnalyzer::beginJob() {
  LogDebug("Phase2TrackerFEDTestAnalyzer") << "[Phase2TrackerFEDTestAnalyzer::" << __func__ << "]";
}
// -----------------------------------------------------------------------------
//

void Phase2TrackerFEDTestAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& es) {
  // fetch cabling from event setup
  cabling_ = &es.getData(ph2CablingESToken_);

  // FIXME: build map of stacks to compensate for missing trackertopology methods
  tkGeom_ = &es.getData(geomToken_);
  tTopo_ = &es.getData(topoToken_);

  for (auto iu = tkGeom_->detUnits().begin(); iu != tkGeom_->detUnits().end(); ++iu) {
    unsigned int detId_raw = (*iu)->geographicalId().rawId();
    DetId detId = DetId(detId_raw);
    if (detId.det() == DetId::Detector::Tracker) {
      // build map of upper and lower for each module
      if (tTopo_->isLower(detId) != 0) {
        stackMap_[tTopo_->stack(detId)].first = detId;
      }
      if (tTopo_->isUpper(detId) != 0) {
        stackMap_[tTopo_->stack(detId)].second = detId;
      }
    }
  }  // end loop on detunits
}

// -----------------------------------------------------------------------------
//
void Phase2TrackerFEDTestAnalyzer::endJob() {
  LogDebug("Phase2TrackerFEDTestAnalyzer") << "[Phase2TrackerFEDTestAnalyzer::" << __func__ << "]";
}

// -----------------------------------------------------------------------------
//
void Phase2TrackerFEDTestAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  // Retrieve FEDRawData collection
  edm::Handle<FEDRawDataCollection> buffers;
  event.getByToken(token_, buffers);

  // Analyze strip tracker FED buffers in data
  std::vector<int> feds = cabling_->listFeds();
  for (int fedIndex : feds) {
    const FEDRawData& fed = buffers->FEDData(fedIndex);
    if (fed.size() == 0)
      continue;
    if (fed.size() != 0 && fedIndex >= Phase2Tracker::FED_ID_MIN && fedIndex <= Phase2Tracker::FED_ID_MAX) {
      // construct buffer
      Phase2Tracker::Phase2TrackerFEDBuffer buffer(fed.data(), fed.size());
      // Skip FED if buffer is not a valid tracker FEDBuffer
      if (buffer.isValid() == 0) {
        LogTrace("Phase2TrackerDigiProducer") << "[Phase2Tracker::Phase2TrackerDigiProducer::" << __func__ << "]: \n";
        LogTrace("Phase2TrackerDigiProducer") << "Skipping invalid buffer for FED nr " << fedIndex << endl;
        continue;
      }

      LOGPRINT << " -------------------------------------------- ";
      LOGPRINT << " buffer debug ------------------------------- ";
      LOGPRINT << " -------------------------------------------- ";
      LOGPRINT << " buffer size : " << buffer.bufferSize();
      LOGPRINT << " fed id      : " << fedIndex;
      LOGPRINT << " -------------------------------------------- ";
      LOGPRINT << " tracker header debug ------------------------";
      LOGPRINT << " -------------------------------------------- ";

      Phase2TrackerFEDHeader tr_header = buffer.trackerHeader();
      LOGPRINT << " Version  : " << hex << setw(2) << (int)tr_header.getDataFormatVersion();
      LOGPRINT << " Mode     : " << hex << setw(2) << (int)tr_header.getDebugMode();
      LOGPRINT << " Type     : " << hex << setw(2) << (int)tr_header.getEventType();
      LOGPRINT << " Readout  : " << hex << setw(2) << (int)tr_header.getReadoutMode();
      LOGPRINT << " Status   : " << hex << setw(16) << (int)tr_header.getGlibStatusCode();
      LOGPRINT << " FE stat  : ";
      for (int i = MAX_FE_PER_FED - 1; i >= 0; i--) {
//       for (int i = 15; i >= 0; i--) {
        if ((tr_header.frontendStatus())[i]) {
          LOGPRINT << "1";
        } else {
          LOGPRINT << "0";
        }
      }
      LOGPRINT << endl;
      LOGPRINT << " Nr CBC   : " << hex << setw(16) << (int)tr_header.getNumberOfCBC() << endl;
      LOGPRINT << " FE/Chip status : ";
      std::vector<Phase2TrackerFEDFEDebug> all_fe_debug = tr_header.CBCStatus();
      std::vector<Phase2TrackerFEDFEDebug>::iterator FE_it;
      for (FE_it = all_fe_debug.begin(); FE_it < all_fe_debug.end(); FE_it++) {
        if (FE_it->IsOn()) {
          LOGPRINT << " FE L1ID: " << endl;
          LOGPRINT << "    " << hex << setw(4) << FE_it->getFEL1ID()[0] << dec << endl;
          LOGPRINT << "    " << hex << setw(4) << FE_it->getFEL1ID()[1] << dec << endl;
          for (int i = 0; i < 16; i++) {
            LOGPRINT << " Chip Error" << hex << setw(1) << FE_it->getChipError(i) << dec << endl;
            LOGPRINT << " Chip L1ID " << hex << setw(4) << FE_it->getChipL1ID(i) << dec << endl;
            LOGPRINT << " Chip PA   " << hex << setw(4) << FE_it->getChipPipelineAddress(i) << dec << endl;
          }
        }
      }
      LOGPRINT << endl;
      LOGPRINT << " -------------------------------------------- " << endl;
      LOGPRINT << " Payload  ----------------------------------- " << endl;
      LOGPRINT << " -------------------------------------------- " << endl;

      // loop channels
      int ichan = 0;
      for (int ife = 0; ife < 16; ife++) {
        for (int icbc = 0; icbc < 16; icbc++) {
          const Phase2TrackerFEDChannel& channel = buffer.channel(ichan);
          if (channel.length() > 0) {
            LOGPRINT << dec << " reading channel : " << icbc << " on FE " << ife;
            LOGPRINT << dec << " with length  : " << (int)channel.length();
            Phase2TrackerFEDRawChannelUnpacker unpacker = Phase2TrackerFEDRawChannelUnpacker(channel);
            while (unpacker.hasData()) {
              LOGPRINT << (unpacker.stripOn() ? "1" : "_");
              unpacker++;
            }
            LOGPRINT << "\n";
          }
          ichan++;   
        }
      }  // end loop on channels
    }
  }
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Phase2TrackerFEDTestAnalyzer);
