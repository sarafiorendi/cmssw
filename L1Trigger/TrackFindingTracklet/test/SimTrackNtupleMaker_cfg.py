############################################################
# Minimal cfg to run SimTrackNtupleMaker.
#
# This reads an existing GEN-SIM-DIGI-RAW (or later) file that already
# contains the "mix"/"MergedTrackTruth" TrackingParticle collection
# (i.e. any file usable as input to L1TrackNtupleMaker_cfg.py already
# has this) and just dumps every TrackingParticle, no cuts, with its
# parent's pdgId, into TPNtuple.root:eventTree.
#
# It does NOT run any L1 tracking, stub, or DTC sequences -- this is
# deliberately just the TrackingParticle collection, nothing else.
############################################################

import FWCore.ParameterSet.Config as cms

process = cms.Process("TPNtuple")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '133X_mcRun4_realistic_v1', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10))

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        # point this at your displacedSUSY GEN-SIM-DIGI-RAW (or later) file(s)
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/29321e3e-1f9c-4b24-b5d8-5f2d04afcb18.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/313a7230-0ca3-43a0-90c7-d9b19e602b49.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/36f286ec-3efd-4d81-8a93-c94fb0f58772.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/385f32c7-f907-4395-903f-4a0d1d54aefb.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/45e2c8d0-0602-4983-9c29-4b4870e9b663.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/4a2d4946-3139-4fac-a158-d8faceee69a2.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/52aea98f-af10-4484-977c-8a997a5a52a7.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/60d96272-1dd9-4c3b-93ee-5e55954c49e6.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/6e9856df-76e8-4c8e-9ac4-fa2772e9430d.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/7c3101fe-8265-497f-9d2a-f0289ac13d66.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/86e055a5-0dc4-46e2-a7b0-07db90c81b0d.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/9063604d-4ff5-49fd-993f-49388cf44863.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/a1d1db3c-bf59-4dc7-8342-7692a589af34.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/a9204407-39a0-43e5-bd5a-c3d600eab045.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/ade8ad66-3f91-430d-9d34-42c37ddff106.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/be339224-92ca-4797-871b-3e7769193908.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/cb26fad2-19b9-43d6-b6db-01ea7028145e.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/d0c02317-e9d6-4741-8bec-cb14867f61c2.root',
    '/store/relval/CMSSW_20_0_0_pre1/RelValDisplacedSUSY_14TeV/GEN-SIM/150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v3/2590000/d55258e9-72fa-4375-83c9-a4c27f9dde25.root',
    ),
)

process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("TPNtuple_SimTrack_HGGH.root"), closeFileFast=cms.untracked.bool(True)
)

process.tpNtuple = cms.EDAnalyzer(
    "SimTrackNtupleMaker",
    SimTrackInputTag=cms.InputTag("g4SimHits", ""),
    GenParticleInputTag=cms.InputTag("genParticles", ""),
#     TrackingParticleInputTag=cms.InputTag("mix", "MergedTrackTruth"),
)

process.p = cms.Path(process.tpNtuple)
process.schedule = cms.Schedule(process.p)