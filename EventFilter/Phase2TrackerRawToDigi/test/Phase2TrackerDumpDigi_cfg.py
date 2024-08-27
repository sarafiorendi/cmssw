import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("DumpDigi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Input source
# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring("/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt1p5To8/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/3ce31040-55a5-4469-8ee2-16d050bb6ade.root")
#  )

process.load('P2TrackerCabling_cfi')
process.Phase2TrackerDumpDigi = cms.EDAnalyzer(
    'Phase2TrackerDumpDigi',
    ProductLabel = cms.InputTag("siPhase2Clusters")
)
### were created by 
### https://github.com/cms-sw/cmssw/blob/bbbd522740e9dec3d103ceadd44fc0361310e22d/RecoLocalTracker/SiPhase2Clusterizer/python/phase2TrackerClusterizer_cfi.py#L4

### test our digis after digi-raw-digi ###
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:raw2digi.root")    
)
process.Phase2TrackerDumpDigi.ProductLabel = cms.InputTag("Phase2TrackerDigiProducer", "Sparsified", "RawToDigi")
### end test our digis after digi-raw-digi ###


## Geometry D98
## from https://github.com/cms-sw/cmssw/blob/master/Configuration/Geometry/python/dict2026Geometry.py
##    ("O9","T32","C17","M10","F8","I16") : "D98",
## https://github.com/cms-sw/cmssw/blob/master/Configuration/Geometry/python/dict2026Geometry.py#L221C14-L221C86
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

# Temporary change until we switch to D110 geometry.
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun4_realistic_v3', '')
process.p = cms.EndPath(process.Phase2TrackerDumpDigi)