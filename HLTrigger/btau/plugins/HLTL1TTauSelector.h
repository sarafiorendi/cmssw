#ifndef HLTTrigger_HLTL1TTauSelector_HLTL1TTauSelector_H
#define HLTTrigger_HLTL1TTauSelector_HLTL1TTauSelector_H

//-------------------------------------------------
//
/**  \class HLTL1TTauSelector
 * 
 *   HLTL1TTauSelector:
 *   Simple selector to output a subset of L1 muon collection 
 *   
 *   based on HLTL1TTauSelector
 *
 *
 */
//
//--------------------------------------------------

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"

// Data Formats
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}  // namespace edm

class HLTL1TTauSelector : public edm::global::EDProducer<> {
public:
  /// Constructor
  explicit HLTL1TTauSelector(const edm::ParameterSet&);

  /// Destructor
  ~HLTL1TTauSelector() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  edm::InputTag theSource;

  edm::EDGetTokenT<l1t::TauBxCollection> tauCollToken_;

  const double theL1MinPt;
  const double theL1MaxEta;

  /// use central bx only taus
  bool centralBxOnly_;
};

#endif
