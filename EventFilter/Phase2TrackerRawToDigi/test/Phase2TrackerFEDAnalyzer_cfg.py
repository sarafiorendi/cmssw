import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("RawToDigi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
# use this to read testbeam .dat files
# process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring( 'file:digi2raw.root')
)


# use this to use hand-made testbeam cabling
process.load('DummyCablingTxt_cfi')

process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')

process.Phase2TrackerRawAna = cms.EDAnalyzer(
    'Phase2TrackerFEDTestAnalyzer',
    ProductLabel = cms.InputTag("Phase2TrackerDigiToRawProducer")
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# Temporary change until we switch to D110 geometry.
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun4_realistic_v3', '')


# process.p = cms.Path(process.Phase2TrackerRawAna)
process.e = cms.EndPath(process.Phase2TrackerRawAna)
