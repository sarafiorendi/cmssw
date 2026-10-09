#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TTree.h>
#include <cmath>
#include <vector>

class SimTrackNtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SimTrackNtupleMaker(const edm::ParameterSet& iConfig);
  ~SimTrackNtupleMaker() override = default;

private:
  void beginJob() override;
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override {}

  edm::InputTag SimTrackInputTag_;
  edm::EDGetTokenT<std::vector<SimTrack>> SimTrackToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> SimVertexToken_;
  edm::InputTag GenParticleInputTag_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticleToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> getTokenBField_;

  bool available_ = false;
  TTree* eventTree_ = nullptr;

  std::vector<float> m_st_pt;
  std::vector<float> m_st_eta;
  std::vector<float> m_st_phi;
  std::vector<int> m_st_charge;
  std::vector<int> m_st_pdgid;
  std::vector<int> m_st_mompdgid;
  std::vector<int> m_st_status;         // status of the associated GenParticle
  std::vector<int> m_st_eventid;        // 0 = signal event, nonzero = pileup interaction index / valid for TrackingParticle, not available for simTrack
  std::vector<int> m_st_bunchcrossing;  // in-time (0) vs out-of-time pileup bunch crossing / valid for TrackingParticle, not available for simTrack?
  std::vector<float> m_st_vx;           // production vertex position
  std::vector<float> m_st_vy;
  std::vector<float> m_st_vz;
  std::vector<float> m_st_lxy;          // not used
  std::vector<float> m_st_lz;           // not used
};

SimTrackNtupleMaker::SimTrackNtupleMaker(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");

  SimTrackInputTag_ = iConfig.getParameter<edm::InputTag>("SimTrackInputTag");
  SimTrackToken_ = consumes<std::vector<SimTrack>>(SimTrackInputTag_);
  SimVertexToken_ = consumes<std::vector<SimVertex>>(SimTrackInputTag_);
  GenParticleInputTag_ = iConfig.getParameter<edm::InputTag>("GenParticleInputTag");
  GenParticleToken_ = consumes<std::vector<reco::GenParticle>>(GenParticleInputTag_);
  getTokenBField_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
}

void SimTrackNtupleMaker::beginJob() {
  edm::Service<TFileService> fs;
  available_ = fs.isAvailable();
  if (!available_)
    return;  // No ROOT file open.

  eventTree_ = fs->make<TTree>("eventTree", "All SimTracks, no cuts, with parent pdgId");
  eventTree_->Branch("st_pt", &m_st_pt);
  eventTree_->Branch("st_eta", &m_st_eta);
  eventTree_->Branch("st_phi", &m_st_phi);
  eventTree_->Branch("st_charge", &m_st_charge);
  eventTree_->Branch("st_pdgid", &m_st_pdgid);
  eventTree_->Branch("st_parent_pdgid", &m_st_mompdgid);
  eventTree_->Branch("st_status", &m_st_status);
  eventTree_->Branch("st_eventid", &m_st_eventid);
  eventTree_->Branch("st_bunchcrossing", &m_st_bunchcrossing);
  eventTree_->Branch("st_vx", &m_st_vx);
  eventTree_->Branch("st_vy", &m_st_vy);
  eventTree_->Branch("st_vz", &m_st_vz);
  eventTree_->Branch("st_lxy", &m_st_lxy);
  eventTree_->Branch("st_lz", &m_st_lz);
}

void SimTrackNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (!available_)
    return;  

  m_st_pt.clear();
  m_st_eta.clear();
  m_st_phi.clear();
  m_st_charge.clear();
  m_st_pdgid.clear();
  m_st_mompdgid.clear();
  m_st_status.clear();
  m_st_eventid.clear();
  m_st_bunchcrossing.clear();
  m_st_vx.clear();
  m_st_vy.clear();
  m_st_vz.clear();
  m_st_lxy.clear();
  m_st_lz.clear();

  edm::Handle<std::vector<SimTrack>> SimTrackHandle;
  iEvent.getByToken(SimTrackToken_, SimTrackHandle);

  edm::Handle<std::vector<SimVertex>> SimVertexHandle;
  iEvent.getByToken(SimVertexToken_, SimVertexHandle);

  edm::Handle<std::vector<reco::GenParticle>> GenParticleHandle;
  iEvent.getByToken(GenParticleToken_, GenParticleHandle);
  
  for (const auto& tp : *SimTrackHandle) {
    const float tmp_st_pt = tp.momentum().pt();
    const float tmp_st_charge = tp.charge();
    const float tmp_st_eta = tp.momentum().eta();
    const float tmp_st_phi = tp.momentum().phi();
    const int gpIndex = tp.genpartIndex() - 1;
    if (gpIndex < 0 || gpIndex > int(GenParticleHandle->size()))
      continue;
    auto gp = GenParticleHandle->at(gpIndex);
    
    int mom_pdgId = -1;
    int gp_mom_index = -1;
    if (gp.numberOfMothers()>0)
      gp_mom_index = gp.motherRef(0).key();
    if (gp_mom_index > 0 && gp_mom_index < int(GenParticleHandle->size())){
      auto gp_mom = GenParticleHandle->at(gp_mom_index);
      mom_pdgId = gp_mom.pdgId();
    }

    const int tmp_st_pdgid = gp.pdgId();
    const int tmp_st_status = gp.status();
    
    bool proceed = false;
    if ((abs(tmp_st_pdgid) >  1000000 && 
         abs(tmp_st_pdgid) <  1000040) ||
        (abs(tmp_st_pdgid) >  2000000 && 
         abs(tmp_st_pdgid) <  2000016) || 
        (abs(tmp_st_pdgid) >  1570 && 
         abs(tmp_st_pdgid) <  1578) || 
        (abs(tmp_st_pdgid) >  10 && 
         abs(tmp_st_pdgid) <  16) 
       ) 
        proceed = true;

    if ((abs(mom_pdgId) >  1000000 && 
         abs(mom_pdgId) <  1000040) ||
        (abs(mom_pdgId) >  2000000 && 
         abs(mom_pdgId) <  2000016) || 
        (abs(mom_pdgId) >  1570 && 
         abs(mom_pdgId) <  1578) ||
        (abs(mom_pdgId) >  10 && 
         abs(mom_pdgId) <  16) ||
        (abs(mom_pdgId) >  24 && 
         abs(mom_pdgId) <  36) 
       ) 
        proceed = true;
    if (!proceed) continue;    

    auto vertexIndex = tp.vertIndex();
    if (vertexIndex > -1 && vertexIndex < int(SimVertexHandle->size())) {
        auto sim_vtx = SimVertexHandle->at(vertexIndex);
        m_st_vx.push_back(sim_vtx.position().x());
        m_st_vy.push_back(sim_vtx.position().y());
        m_st_vz.push_back(sim_vtx.position().z());
    }

    m_st_pt.push_back(tmp_st_pt);
    m_st_eta.push_back(tmp_st_eta);
    m_st_phi.push_back(tmp_st_phi);
    m_st_charge.push_back(static_cast<int>(tmp_st_charge));
    m_st_pdgid.push_back(tmp_st_pdgid);
    m_st_mompdgid.push_back(mom_pdgId);
    m_st_status.push_back(tmp_st_status);
    
//     m_st_eventid.push_back(tp.eventId().event());
//     m_st_bunchcrossing.push_back(tp.eventId().bunchCrossing());

  }  // end loop over sim tracks
  eventTree_->Fill();
}

DEFINE_FWK_MODULE(SimTrackNtupleMaker);
