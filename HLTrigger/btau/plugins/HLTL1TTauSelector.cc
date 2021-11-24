#include "HLTL1TTauSelector.h"

// Framework
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

using namespace std;
using namespace edm;
using namespace l1t;

// constructors
HLTL1TTauSelector::HLTL1TTauSelector(const edm::ParameterSet& iConfig)
    : theSource(iConfig.getParameter<InputTag>("InputObjects")),
      theL1MinPt(iConfig.getParameter<double>("L1MinPt")),
      theL1MaxEta(iConfig.getParameter<double>("L1MaxEta")),
      centralBxOnly_(iConfig.getParameter<bool>("CentralBxOnly")) {
  tauCollToken_ = consumes<TauBxCollection>(theSource);

  produces<TauBxCollection>();
}

// destructor
HLTL1TTauSelector::~HLTL1TTauSelector() = default;

void HLTL1TTauSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("InputObjects", edm::InputTag("hltGmtStage2Digis"));
  desc.add<double>("L1MinPt", -1.);
  desc.add<double>("L1MaxEta", 5.0);
  desc.add<bool>("CentralBxOnly", true);
  descriptions.add("hltL1TTauSelector", desc);
}

void HLTL1TTauSelector::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  const std::string metname = "Tau|HLTL1TTauSelector";

  unique_ptr<TauBxCollection> output(new TauBxCollection());

  // Tau particles
  edm::Handle<TauBxCollection> tauColl;
  iEvent.getByToken(tauCollToken_, tauColl);
  LogTrace(metname) << "Number of taus " << tauColl->size() << endl;

  for (int ibx = tauColl->getFirstBX(); ibx <= tauColl->getLastBX(); ++ibx) {
    if (centralBxOnly_ && (ibx != 0))
      continue;
    for (auto it = tauColl->begin(ibx); it != tauColl->end(ibx); it++) {

      float pt = it->pt();
      float eta = it->eta();

      if (pt < theL1MinPt || fabs(eta) > theL1MaxEta)
        continue;

      LogTrace(metname) << "L1 Tau Found";
      LogTrace(metname) << "Pt = " << pt << " GeV/c";
      LogTrace(metname) << "eta = " << eta;

      output->push_back(ibx, *it);
    }
  }

  iEvent.put(std::move(output));
}