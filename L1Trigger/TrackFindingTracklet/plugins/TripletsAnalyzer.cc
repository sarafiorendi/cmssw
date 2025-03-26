// EDAnalyzer producing a small ntuple containing properties of Phase2TrackerCluster1D, 
// to be used to debug the clusters-to-raw and raw-to-cluster steps

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/L1TrackTrigger/interface/L1StubTriplet.h"

#include "TTree.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

class TripletsAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  TripletsAnalyzer(const edm::ParameterSet& pset);
  ~TripletsAnalyzer() override;
  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {};
  void analyze(const edm::Event&, const edm::EventSetup&);

private:
  virtual void beginJob();
  virtual void endJob();
  const edm::EDGetTokenT<std::vector<L1StubTriplet>> inputTripletsToken_;

  edm::Service<TFileService> fs_;
  TTree* outTree_;
  ofstream logfile_;

  long unsigned int eventNumber_;
  
  float inner_r_;
  float middle_r_;
  float outer_r_;

  float inner_z_;
  float middle_z_;
  float outer_z_;

  float inner_bend_;
  float middle_bend_;
  float outer_bend_;

  int inner_rzbin_;
  int middle_rzbin_;
  int outer_rzbin_;

  int inner_index_;
  int middle_index_;
  int outer_index_;

  int inner_layerdisk_;
  int middle_layerdisk_;
  int outer_layerdisk_;

  int sector_;
  int region_;
  int tpdunit_;
  
  int firstbin_out_;
  int firstbin_in_;
  int diffmax_out_;
  int diffmax_in_;
  int rzeff_out_;
  int rzeff_in_;

  int inner_rval_;
  int middle_rval_;
  int outer_rval_;
  
};

TripletsAnalyzer::TripletsAnalyzer(const edm::ParameterSet& pset)
    : 
      inputTripletsToken_(consumes<std::vector<L1StubTriplet>>(pset.getParameter<edm::InputTag>("InputTriplets"))) {
  // Initialize the log file
//   logfile_.open("TripletsAnalyzer_output.txt");
//   if (!logfile_.is_open()) {
//     throw cms::Exception("OutputFileError") << "Failed to open log file for writing.";
//   }
}

TripletsAnalyzer::~TripletsAnalyzer() {
  // Close the log file
//   logfile_.close();
}

void TripletsAnalyzer::beginJob() {
  outTree_ = fs_->make<TTree>("tripletsTree", "tripletsTree");

  outTree_->Branch("eventNumber", &eventNumber_, "eventNumber/L");

  outTree_->Branch("inner_r", &inner_r_, "inner_r/F");
  outTree_->Branch("middle_r", &middle_r_, "middle_r/F");
  outTree_->Branch("outer_r", &outer_r_, "outer_stub_r/F");

  outTree_->Branch("inner_z", &inner_z_, "inner_z/F");
  outTree_->Branch("middle_z", &middle_z_, "middle_z/F");
  outTree_->Branch("outer_z", &outer_z_, "outer_stub_z/F");

  outTree_->Branch("inner_bend", &inner_bend_, "inner_bend/F");
  outTree_->Branch("middle_bend", &middle_bend_, "middle_bend/F");
  outTree_->Branch("outer_bend", &outer_bend_, "outer_stub_bend/F");

  outTree_->Branch("inner_rzbin", &inner_rzbin_, "inner_rzbin/I");
  outTree_->Branch("middle_rzbin", &middle_rzbin_, "middle_rzbin/I");
  outTree_->Branch("outer_rzbin", &outer_rzbin_, "outer_stub_rzbin/I");

  outTree_->Branch("inner_index", &inner_index_, "inner_index/i");
  outTree_->Branch("middle_index", &middle_index_, "middle_index/i");
  outTree_->Branch("outer_index", &outer_index_, "outer_stub_index/i");

  outTree_->Branch("inner_layerdisk", &inner_layerdisk_, "inner_layerdisk/i");
  outTree_->Branch("middle_layerdisk", &middle_layerdisk_, "middle_layerdisk/i");
  outTree_->Branch("outer_layerdisk", &outer_layerdisk_, "outer_stub_layerdisk/i");

  outTree_->Branch("sector", &sector_, "sector/i");
  outTree_->Branch("region", &region_, "region/i");
  outTree_->Branch("tpdunit", &tpdunit_, "tpdunit/i");

  outTree_->Branch("firstbin_out", &firstbin_out_, "firstbin_out/I");
  outTree_->Branch("firstbin_in", &firstbin_in_, "firstbin_in/I");
  outTree_->Branch("diffmax_out", &diffmax_out_, "diffmax_out/I");
  outTree_->Branch("diffmax_in", &diffmax_in_, "diffmax_in/I");
  outTree_->Branch("rzeff_out", &rzeff_out_, "rzeff_out/I");
  outTree_->Branch("rzeff_in", &rzeff_in_, "rzeff_in/I");

  outTree_->Branch("inner_rval", &inner_rval_, "inner_rval/I");
  outTree_->Branch("middle_rval", &middle_rval_, "middle_rval/I");
  outTree_->Branch("outer_rval", &outer_rval_, "outer_stub_rval/I");

}
void TripletsAnalyzer::endJob() {
//   outTree_->GetDirectory()->cd();
//   outTree_->Write();
}

void TripletsAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& es) {
}

void TripletsAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& es) {
  edm::Handle<std::vector<L1StubTriplet>> triplets_handle;
  event.getByToken(inputTripletsToken_, triplets_handle);

  std::stringstream output;
  output << "size of triplets: " << triplets_handle.product()->size() << std::endl;
  
  eventNumber_ = event.id().event();

  int count_triplets = 0;
  
  for (const auto& iTriplet : *triplets_handle) {
    inner_r_ = iTriplet.getStubRapprox(0);
    middle_r_ = iTriplet.getStubRapprox(1);
    outer_r_ = iTriplet.getStubRapprox(2);

    inner_z_ = iTriplet.getStubZapprox(0);
    middle_z_ = iTriplet.getStubZapprox(1);
    outer_z_ = iTriplet.getStubZapprox(2);

    inner_bend_ = iTriplet.getStubBend(0);
    middle_bend_ = iTriplet.getStubBend(1);
    outer_bend_ = iTriplet.getStubBend(2);

    inner_rzbin_ = iTriplet.getStubRZbin(0);
    middle_rzbin_ = iTriplet.getStubRZbin(1);
    outer_rzbin_ = iTriplet.getStubRZbin(2);

    inner_index_ = iTriplet.getStubIndex(0);
    middle_index_ = iTriplet.getStubIndex(1);
    outer_index_ = iTriplet.getStubIndex(2);

    inner_layerdisk_ = iTriplet.getStubLayerdisk(0);
    middle_layerdisk_ = iTriplet.getStubLayerdisk(1);
    outer_layerdisk_ = iTriplet.getStubLayerdisk(2);

    sector_ = iTriplet.getSector();
    region_ = iTriplet.getRegion();
    tpdunit_ = iTriplet.getTPDUnit();

    firstbin_in_ = iTriplet.getFirstBinIn();
    diffmax_in_  = iTriplet.getDiffMaxIn();
    rzeff_in_    = iTriplet.getRZEffIn();
    firstbin_out_= iTriplet.getFirstBinOut();
    diffmax_out_ = iTriplet.getDiffMaxOut();
    rzeff_out_   = iTriplet.getRZEffOut();

    inner_rval_ = iTriplet.getStubRValue(0);
    middle_rval_ = iTriplet.getStubRValue(1);
    outer_rval_ = iTriplet.getStubRValue(2);

    outTree_->Fill();  // Fill the tree with current cluster data
  }
  // Output to terminal and log file
//   logfile_ << output.str();
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TripletsAnalyzer);
