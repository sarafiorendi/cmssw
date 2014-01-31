import FWCore.ParameterSet.Config as cms

HLTMuonPFIsoFilter = cms.EDFilter( "HLTMuonPFIsoFilter",
                            saveTags = cms.bool( False ),
                            PreviousCandTag = cms.InputTag( "hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q" ),
                            MinN = cms.int32( 1 ),
                            MaxIso = cms.double( 0.2 ),
                            CandTag = cms.InputTag( "hltL3MuonCandidates" ),
                            onlyCharged = cms.bool(False),
                            applyRhoCorrection = cms.bool(True),
                            EffectiveArea = cms.double(0.22766),
                            RhoTag = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuonsPF"),
                            DepTag = cms.VInputTag( cms.InputTag('hltMuPFIsoValueCharged03'), cms.InputTag('hltMuPFIsoValueGamma03'), cms.InputTag('hltMuPFIsoValueNeutral03') )
                            )
                                   
