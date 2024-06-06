import FWCore.ParameterSet.Config as cms
import math

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
Phase2OTValidateTrackingParticles = DQMEDAnalyzer('Phase2OTValidateTrackingParticles',
    TopFolderName = cms.string('SiTrackerPhase2V'),
    trackingParticleToken = cms.InputTag("mix","MergedTrackTruth"), #tracking particles
    MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"), #truth stub associator
    MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), #truth track associator
    MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterInclusive"), #truth cluster associator
    #MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"), #truth cluster associator
    L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
    L1Tk_minNStub = cms.int32(4),       # L1 tracks with >= 4 stubs
    L1Tk_maxChi2dof = cms.double(25.0),# L1 tracks with Chi2 <= X
    TP_minNStub = cms.int32(4),      # require TP to have >= X number of stubs associated with it
    TP_minNLayersStub = cms.int32(4),   # require TP to have >= X number of layers hit with stubs
    TP_minPt = cms.double(2.0),      # only save TPs with pt > X GeV
    TP_maxEta = cms.double(2.4),     # only save TPs with |eta| < X
    TP_maxVtxZ = cms.double(15.0),     # only save TPs with |z0| < X cm

# deltaZ vs count
    TH1delta_Z = cms.PSet(
        Nbinsx = cms.int32(49),
        xmin = cms.double(-0.45),
        xmax = cms.double(0.45)
        ),

# deltaR vs count
    TH1delta_R = cms.PSet(
        Nbinsx = cms.int32(100),
        xmin = cms.double(0.0),
        xmax = cms.double(0.5)
        ),

# tp_phi vs count
    TH1tp_phi = cms.PSet(
        Nbinsx = cms.int32(100),
        xmin = cms.double(-4),
        xmax = cms.double(4)
        ),

# stub count for fake rate
    TH1Stubs = cms.PSet(
        Nbinsx = cms.int32(19),
        xmin = cms.double(-10.0),
        xmax = cms.double(10.0)
        ),

# TTStub Barrel Layers
    TH1TTStub_Layers = cms.PSet(
        Nbinsx = cms.int32(7),
        xmin = cms.double(0.5),
        xmax = cms.double(7.5)
        ),

# tracking particles vs eta
    TH1TrackParts_Eta = cms.PSet(
        Nbinsx = cms.int32(45),
        xmax = cms.double(3),
        xmin = cms.double(-3)
        ),

# tracking particles vs phi
    TH1TrackParts_Phi = cms.PSet(
        Nbinsx = cms.int32(60),
        xmax = cms.double(math.pi),
        xmin = cms.double(-math.pi)
        ),

# tracking particles vs pT
    TH1TrackParts_Pt = cms.PSet(
        Nbinsx = cms.int32(45),
        xmax = cms.double(100),
        xmin = cms.double(0)
        ),

# tracking particles vs pT_relative
    TH1Res_ptRel = cms.PSet(
        Nbinsx = cms.int32(200),
        xmax = cms.double(0.5),
        xmin = cms.double(-0.5)
        ),

# tracking particles vs pT (for efficiency)
    TH1Effic_pt = cms.PSet(
        Nbinsx = cms.int32(50),
        xmax = cms.double(100),
        xmin = cms.double(0)
        ),

# tracking particles vs pT (for efficiency)
    TH1Effic_pt_zoom = cms.PSet(
        Nbinsx = cms.int32(50),
        xmax = cms.double(10),
        xmin = cms.double(0)
        ),

# tracking particles vs eta (for efficiency)
    TH1Effic_eta = cms.PSet(
        Nbinsx = cms.int32(50),
        xmax = cms.double(2.5),
        xmin = cms.double(-2.5)
        ),

# tracking particles vs d0 (for efficiency)
    TH1Effic_d0 = cms.PSet(
        Nbinsx = cms.int32(50),
        xmax = cms.double(2),
        xmin = cms.double(-2)
        ),

# tracking particles vs VtxR/vxy (for efficiency)
    TH1Effic_VtxR = cms.PSet(
        Nbinsx = cms.int32(50),
        xmax = cms.double(5),
        xmin = cms.double(-5)
        ),

# tracking particles vs z0 (for efficiency)
    TH1Effic_VtxZ = cms.PSet(
        Nbinsx = cms.int32(50),
        xmax = cms.double(30),
        xmin = cms.double(-30)
        ),

# Stub associated to eta
    TH1TP_eta = cms.PSet(
        Nbinsx = cms.int32(50),
        xmax = cms.double(2.5),
        xmin = cms.double(-2.5)
        ),

# Stub associated tp pT
    TH1TP_pt = cms.PSet(
        Nbinsx = cms.int32(39),
        xmax = cms.double(40.),
        xmin = cms.double(0.)
        ),

# Stub associated tp dxy
    TH1TP_dxy = cms.PSet(
        Nbinsx = cms.int32(99),
        xmax = cms.double(100.),
        xmin = cms.double(0.)
        ),

# position of stub in z
    TH1Stub_z = cms.PSet(
        Nbinsx = cms.int32(299),
        xmax = cms.double(300.0),
        xmin = cms.double(-300.0)
        ),

# counts vs stub rawBend
    TH1Stub_rawBend = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(8.0),
        xmin = cms.double(-8.0)
        ),

# counts vs stub bend offset
    TH1Stub_bendOffset = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(6.0),
        xmin = cms.double(-6.0)
        ),

# stub inner cluster position
    TH1Stub_inClusPos = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(120.0),
        xmin = cms.double(0.0)
        ),

# counts vs stub bendFE
    TH1Stub_bendFE = cms.PSet(
        Nbinsx = cms.int32(29),
        xmax = cms.double(15.0),
        xmin = cms.double(-15.0)
        ),

# counts vs track bend
    TH1Track_Bend = cms.PSet(
        Nbinsx = cms.int32(29),
        xmax = cms.double(15.0),
        xmin = cms.double(-15.0)
        ),

# counts vs tp bend
    TH1Stub_tpBend = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(8.0),
        xmin = cms.double(-8.0)
        ),

# events vs num of stubs
    TH1StubInEvent = cms.PSet(
        Nbinsx = cms.int32(16),
        xmax = cms.double(16.0),
        xmin = cms.double(0.0)
        ),

# delta bend
    TH1Bend_Res = cms.PSet(
        Nbinsx = cms.int32(399),
        xmax = cms.double(5.0),
        xmin = cms.double(-5.0)
        ),

# delta z
    TH1Z_Res = cms.PSet(
        Nbinsx = cms.int32(399),
        xmax = cms.double(5.5),
        xmin = cms.double(-5.5)
        ),

# delta z endcap only
    TH1Z_Res_Endcap = cms.PSet(
        Nbinsx = cms.int32(121),
        xmax = cms.double(1.0),
        xmin = cms.double(-1.0)
        ),

# delta phi
    TH1Phi_Res = cms.PSet(
        Nbinsx = cms.int32(79),
        xmax = cms.double(0.1),
        xmin = cms.double(-0.1)
        ),

# 2D histo x coordsB vs x coordsC
    TH2coordsBvscoordsC = cms.PSet(
        Nbinsx = cms.int32(99),
        xmin = cms.double(-80.0),
        xmax = cms.double(80.0),
        Nbinsy = cms.int32(99),
        ymin = cms.double(-80.0),
        ymax = cms.double(80.0)
        ),

# 2D histo tiltAngle vs Z0
    TH2tiltAngleVsZ0 = cms.PSet(
        Nbinsx = cms.int32(49),
        xmin = cms.double(-120.0),
        xmax = cms.double(120.0),
        Nbinsy = cms.int32(27),
        ymin = cms.double(-1.8),
        ymax = cms.double(1.8)
        ),

# 2D histo deltaR vs deltaZ
    TH2DeltaRVsDeltaZ = cms.PSet(
        Nbinsx = cms.int32(99),
        xmin = cms.double(-130.0),
        xmax = cms.double(130.0),
        Nbinsy = cms.int32(64),
        ymin = cms.double(0.0),
        ymax = cms.double(65.0)
        ),

# 2D histo Z0 vs deltaZ
    TH2Z0VsDeltaZ = cms.PSet(
        Nbinsx = cms.int32(99),
        xmin = cms.double(-0.5),
        xmax = cms.double(0.5),
        Nbinsy = cms.int32(59),
        ymin = cms.double(-120.0),
        ymax = cms.double(120.0)
        ),

# 2D histo R0 vs deltaR
    TH2R0VsDeltaR = cms.PSet(
        Nbinsx = cms.int32(49),
        xmin = cms.double(0.0),
        xmax = cms.double(0.5),
        Nbinsy = cms.int32(29),
        ymin = cms.double(0.0),
        ymax = cms.double(60.0)
        ),

# 2D histo track_phi vs stub_phi
    TH2TpPhiVsStubPhi = cms.PSet(
        Nbinsx = cms.int32(99),
        xmin = cms.double(-3.5),
        xmax = cms.double(3.5),
        Nbinsy = cms.int32(99),
        ymin = cms.double(-3.5),
        ymax = cms.double(3.5)
        ),

# 2D histo track_z vs stub_z
    TH2TpZVsStubZ = cms.PSet(
        Nbinsx = cms.int32(59),
        xmin = cms.double(-120.0),
        xmax = cms.double(120.0),
        Nbinsy = cms.int32(59),
        ymin = cms.double(-120.0),
        ymax = cms.double(120.0)
        ),

# 2D histo modMaxR vs modMinR
    TH2MaxRvsMinR = cms.PSet(
        Nbinsx = cms.int32(59),
        xmin = cms.double(0.0),
        xmax = cms.double(120.0),
        Nbinsy = cms.int32(59),
        ymin = cms.double(0.0),
        ymax = cms.double(120.0)
        ),

# 2D histo track bend vs stub bend
    TH2TrackVsStub = cms.PSet(
        Nbinsx = cms.int32(29), # Binning for track bend
        xmin = cms.double(-7.0),
        xmax = cms.double(7.0),
        Nbinsy = cms.int32(29), # Binning for stub bend
        ymin = cms.double(-7.0),
        ymax = cms.double(7.0)
        ),

# genuine stubs in barrel layers
    TH1Barrel_Layers = cms.PSet(
        Nbinsx = cms.int32(6),
        xmax = cms.double(7.0),
        xmin = cms.double(1.0)
        ),

# genuine stubs in endcap layers
    TH1Endcap_Layers = cms.PSet(
        Nbinsx = cms.int32(6),
        xmax = cms.double(6.5),
        xmin = cms.double(0.5)
        ),

# tracking particles vs relative pT (for resolution plots)
    TH1Res_pt = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(0.2),
        xmin = cms.double(-0.2)
        ),

# tracking particles vs eta (for resolution)
    TH1Res_eta = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(0.01),
        xmin = cms.double(-0.01)
        ),

# tracking particles vs phi (for resolution)
    TH1Res_phi = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(0.01),
        xmin = cms.double(-0.01)
        ),

# tracking particles vs phi (for resolution)
    TH1StubRes_phi = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(0.01),
        xmin = cms.double(-0.01)
        ),

# tracking particles vs z0 (for resolution)
    TH1Res_VtxZ = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(1.0),
        xmin = cms.double(-1.0)
        ),

# tracking particles vs d0 (for resolution)
    TH1Res_d0 = cms.PSet(
        Nbinsx = cms.int32(100),
        xmax = cms.double(0.05),
        xmin = cms.double(-0.05)
        ),
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(OuterTrackerMonitorTrackingParticles, trackingParticleToken = "mixData:MergedTrackTruth")