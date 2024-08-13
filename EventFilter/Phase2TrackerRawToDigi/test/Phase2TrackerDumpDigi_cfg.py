import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("DigiToRaw")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )

# Input source
process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring("/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/0b2b0b0b-f312-48a8-9d46-ccbadc69bbfd.root")
    fileNames = cms.untracked.vstring("/store/relval/CMSSW_14_0_0_pre2/RelValDisplacedSingleMuFlatPt1p5To8/GEN-SIM-DIGI-RAW/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/3ce31040-55a5-4469-8ee2-16d050bb6ade.root")
    
 )

process.load('DummyCablingTxt_cfi')

process.Phase2TrackerDumpDigi = cms.EDAnalyzer(
    'Phase2TrackerDumpDigi',
    ProductLabel = cms.InputTag("siPhase2Clusters")
)

process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# Temporary change until we switch to D110 geometry.
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun4_realistic_v3', '')

# process.out = cms.OutputModule(
#     "PoolOutputModule",
#     fileName = cms.untracked.string('digi2raw.root'),
#     outputCommands = cms.untracked.vstring(
#       # 'drop *',
#       'keep *_Phase2TrackerDigiToRawProducer_*_*'
#       )
#     )

process.p = cms.EndPath(process.Phase2TrackerDumpDigi)

# process.e = cms.EndPath(process.out)

# Automatic addition of the customisation function from #SLHCUpgradeSimulations.Configuration.combinedCustoms
#from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023tilted4021
#call to customisation function cust_2023tilted4021 imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
#process = cust_2023tilted4021(process)

