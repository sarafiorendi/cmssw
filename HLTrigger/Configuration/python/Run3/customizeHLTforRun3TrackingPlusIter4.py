import copy
import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.common import *
from HLTrigger.Configuration.customizeHLTforPatatrack import *
from Configuration.ProcessModifiers.pixelNtupletFit_cff import pixelNtupletFit
from Configuration.ProcessModifiers.gpu_cff import gpu

def customizeHLTforRun3Tracking(process):
    
    process.extend(pixelNtupletFit)
    process.extend(gpu)

    process = customizeHLTforPatatrackTriplets(process)    
    if hasattr(process,'hltPixelTracksCUDA'):
        process.hltPixelTracksCUDA.includeJumpingForwardDoublets = cms.bool(True)
        process.hltPixelTracksCUDA.idealConditions               = cms.bool(False)
        process.hltPixelTracksCUDA.fillStatistics                = cms.bool(True)
        process.hltPixelTracksCUDA.useSimpleTripletCleaner       = cms.bool(False)
    if hasattr(process,'hltPixelTracksSoA'):
        process.hltPixelTracksSoA.cpu.includeJumpingForwardDoublets = cms.bool(True)
        process.hltPixelTracksSoA.cpu.idealConditions               = cms.bool(False)
        process.hltPixelTracksSoA.cpu.fillStatistics                = cms.bool(True)
        process.hltPixelTracksSoA.cpu.useSimpleTripletCleaner       = cms.bool(False)

    if hasattr(process,'hltPixelTracks'):
        process.hltPixelTracks.minNumberOfHits = cms.int32(0)
        process.hltPixelTracks.minQuality = cms.string('loose')

    if hasattr(process,'HLTIter0PSetTrajectoryFilterIT'):
        process.HLTIter0PSetTrajectoryFilterIT.minHitsMinPt        = cms.int32(3)
        process.HLTIter0PSetTrajectoryFilterIT.minimumNumberOfHits = cms.int32(3)

    if hasattr(process,'hltSiStripRawToClustersFacility'):
        process.hltSiStripRawToClustersFacility.onDemand = cms.bool( True )

    if hasattr(process,'hltIter0PFLowPixelSeedsFromPixelTracks'):
        process.hltIter0PFLowPixelSeedsFromPixelTracks.includeFourthHit = cms.bool(True)

    if hasattr(process,'hltIter0PFlowTrackCutClassifier'):
        process.hltIter0PFlowTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
            src = cms.InputTag("hltIter0PFlowCtfWithMaterialTracks"),
            beamspot = cms.InputTag("hltOnlineBeamSpot"),
            vertices = cms.InputTag("hltTrimmedPixelVertices"),
            qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
            mva = cms.PSet(
                minPixelHits = cms.vint32(0, 0, 0),
                maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
                dr_par = cms.PSet(
                    d0err = cms.vdouble(0.003, 0.003, 0.003),
                    dr_par2 = cms.vdouble(3.40282346639e+38, 0.6, 0.6),
                    dr_par1 = cms.vdouble(3.40282346639e+38, 0.8, 0.8),
                    dr_exp = cms.vint32(4, 4, 4),
                    d0err_par = cms.vdouble(0.001, 0.001, 0.001)
                ),
                maxLostLayers = cms.vint32(1, 1, 1),
                min3DLayers = cms.vint32(0, 0, 0),
                dz_par = cms.PSet(
                    dz_par1 = cms.vdouble(3.40282346639e+38, 0.75, 0.75),
                    dz_par2 = cms.vdouble(3.40282346639e+38, 0.5, 0.5),
                    dz_exp = cms.vint32(4, 4, 4)
                ),
                minNVtxTrk = cms.int32(3),
                maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
                minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
                maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
                maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
                maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
                minLayers = cms.vint32(3, 3, 3)
            ),
            ignoreVertices = cms.bool(False)
        )
    
    if hasattr(process,'hltMergedTracks'):
        process.hltMergedTracks = process.hltIter0PFlowTrackSelectionHighPurity.clone()

    process.HLTIterativeTrackingIteration0Task = cms.Sequence(
        process.hltIter0PFLowPixelSeedsFromPixelTracks +
        process.hltIter0PFlowCkfTrackCandidates +
        process.hltIter0PFlowCtfWithMaterialTracks +
        process.hltIter0PFlowTrackCutClassifier +
        process.hltMergedTracks
    )
    if hasattr(process,'HLTIterativeTrackingIteration0'):
        process.HLTIterativeTrackingIteration0 = cms.Sequence( process.HLTIterativeTrackingIteration0Task )
    
    if hasattr(process,'HLTIterativeTrackingIter02'):
        process.HLTIterativeTrackingIter02 = cms.Sequence( process.HLTIterativeTrackingIteration0 )
    
    if hasattr(process,'MC_ReducedIterativeTracking_v12'):
        process.MC_ReducedIterativeTracking_v12 = cms.Path( 
            process.HLTBeginSequence +
            process.hltPreMCReducedIterativeTracking +
            process.HLTDoLocalPixelSequence +
            process.HLTRecopixelvertexingSequence +
            process.HLTDoLocalStripSequence +
            process.HLTIterativeTrackingIter02 +
            process.HLTEndSequence
        )

    return process


### Default standard tracking + displaced iter-4 for taus using hltTrimmedPixelVertices
def customizeHLTforRun3TrackingPlusIter4ForTau(process):

    process.extend(pixelNtupletFit)
    process.extend(gpu)

    process = customizeHLTforPatatrackTriplets(process)    
    if hasattr(process,'hltPixelTracksCUDA'):
        process.hltPixelTracksCUDA.includeJumpingForwardDoublets = cms.bool(True)
        process.hltPixelTracksCUDA.idealConditions               = cms.bool(False)
        process.hltPixelTracksCUDA.fillStatistics                = cms.bool(True)
        process.hltPixelTracksCUDA.useSimpleTripletCleaner       = cms.bool(False)
    if hasattr(process,'hltPixelTracksSoA'):
        process.hltPixelTracksSoA.cpu.includeJumpingForwardDoublets = cms.bool(True)
        process.hltPixelTracksSoA.cpu.idealConditions               = cms.bool(False)
        process.hltPixelTracksSoA.cpu.fillStatistics                = cms.bool(True)
        process.hltPixelTracksSoA.cpu.useSimpleTripletCleaner       = cms.bool(False)

    if hasattr(process,'hltPixelTracks'):
        process.hltPixelTracks.minNumberOfHits = cms.int32(0)
        process.hltPixelTracks.minQuality = cms.string('loose')

    if hasattr(process,'HLTIter0PSetTrajectoryFilterIT'):
        process.HLTIter0PSetTrajectoryFilterIT.minHitsMinPt        = cms.int32(3)
        process.HLTIter0PSetTrajectoryFilterIT.minimumNumberOfHits = cms.int32(3)

    if hasattr(process,'hltSiStripRawToClustersFacility'):
        process.hltSiStripRawToClustersFacility.onDemand = cms.bool( True )

    if hasattr(process,'hltIter0PFLowPixelSeedsFromPixelTracks'):
        process.hltIter0PFLowPixelSeedsFromPixelTracks.includeFourthHit = cms.bool(True)

    if hasattr(process,'hltIter0PFlowTrackCutClassifier'):
        process.hltIter0PFlowTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
            src = cms.InputTag("hltIter0PFlowCtfWithMaterialTracks"),
            beamspot = cms.InputTag("hltOnlineBeamSpot"),
            vertices = cms.InputTag("hltTrimmedPixelVertices"),
            qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
            mva = cms.PSet(
                minPixelHits = cms.vint32(0, 0, 0),
                maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
                dr_par = cms.PSet(
                    d0err = cms.vdouble(0.003, 0.003, 0.003),
                    dr_par2 = cms.vdouble(3.40282346639e+38, 0.6, 0.6),
                    dr_par1 = cms.vdouble(3.40282346639e+38, 0.8, 0.8),
                    dr_exp = cms.vint32(4, 4, 4),
                    d0err_par = cms.vdouble(0.001, 0.001, 0.001)
                ),
                maxLostLayers = cms.vint32(1, 1, 1),
                min3DLayers = cms.vint32(0, 0, 0),
                dz_par = cms.PSet(
                    dz_par1 = cms.vdouble(3.40282346639e+38, 0.75, 0.75),
                    dz_par2 = cms.vdouble(3.40282346639e+38, 0.5, 0.5),
                    dz_exp = cms.vint32(4, 4, 4)
                ),
                minNVtxTrk = cms.int32(3),
                maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
                minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
                maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
                maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
                maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
                minLayers = cms.vint32(3, 3, 3)
            ),
            ignoreVertices = cms.bool(False)
        )
    
    if hasattr(process,'hltMergedTracks'):
        process.hltMergedTracks = process.hltIter0PFlowTrackSelectionHighPurity.clone()

    process.HLTIterativeTrackingIteration0Task = cms.Sequence(
        process.hltIter0PFLowPixelSeedsFromPixelTracks +
        process.hltIter0PFlowCkfTrackCandidates +
        process.hltIter0PFlowCtfWithMaterialTracks +
        process.hltIter0PFlowTrackCutClassifier +
        process.hltMergedTracks
    )
    if hasattr(process,'HLTIterativeTrackingIteration0'):
        process.HLTIterativeTrackingIteration0 = cms.Sequence( process.HLTIterativeTrackingIteration0Task )
    
    if hasattr(process,'HLTIterativeTrackingIter02'):
        process.HLTIterativeTrackingIter02 = cms.Sequence( process.HLTIterativeTrackingIteration0 )
    
    if hasattr(process,'MC_ReducedIterativeTracking_v12'):
        process.MC_ReducedIterativeTracking_v12 = cms.Path( 
            process.HLTBeginSequence +
            process.hltPreMCReducedIterativeTracking +
            process.HLTDoLocalPixelSequence +
            process.HLTRecopixelvertexingSequence +
            process.HLTDoLocalStripSequence +
            process.HLTIterativeTrackingIter02 +
            process.HLTEndSequence
        )
    ### Iter-4
    process.hltDisplacedhltIter4ClustersRefRemovalForTau = cms.EDProducer( 
        "TrackClusterRemover",
        trajectories = cms.InputTag( "hltMergedTracks" ),
        trackClassifier = cms.InputTag( '','QualityMasks' ),
        pixelClusters = cms.InputTag( "hltSiPixelClusters" ),
        stripClusters = cms.InputTag( "hltSiStripRawToClustersFacility" ),
        oldClusterRemovalInfo = cms.InputTag( "" ),
        TrackQuality = cms.string( "highPurity" ),
        maxChi2 = cms.double( 16.0 ),
        minNumberOfLayersWithMeasBeforeFiltering = cms.int32( 0 ),
        overrideTrkQuals = cms.InputTag( "" )
    )
    process.hltDisplacedhltIter4MaskedMeasurementTrackerEventForTau = cms.EDProducer( 
        "MaskedMeasurementTrackerEventProducer",
        src = cms.InputTag( "hltSiStripClusters" ),
        OnDemand = cms.bool( False ),
        clustersToSkip = cms.InputTag( "hltDisplacedhltIter4ClustersRefRemovalForTau" )
    )
    process.hltDisplacedhltIter4PixelLessLayerTripletsForTau = cms.EDProducer( 
        "SeedingLayersEDProducer",
        layerList = cms.vstring( 'TIB1+TIB2+MTIB3',        # 250 340 430
                                 'TIB1+TIB2+MTID1_pos',    # 250 340 277
                                 'TIB1+TIB2+MTID1_neg',
                                 'TID1_pos+TID2_pos+TID3_pos',   # 277 367 447
                                 'TID1_neg+TID2_neg+TID3_neg',
                                 'TID1_pos+TID2_pos+MTID3_pos',   
                                 'TID1_neg+TID2_neg+MTID3_neg' ),
        BPix = cms.PSet(  ),
        FPix = cms.PSet(  ),
        TIB = cms.PSet( 
            TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
            skipClusters = cms.InputTag( "hltDisplacedhltIter4ClustersRefRemovalForTau" ),
            clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) )
        ),
        TID = cms.PSet( 
            TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
            minRing = cms.int32( 1 ),
            skipClusters = cms.InputTag( "hltDisplacedhltIter4ClustersRefRemovalForTau" ),
            useRingSlector = cms.bool( True ),
            clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
            maxRing = cms.int32( 2 )
        ),
        TOB = cms.PSet(  ),
        TEC = cms.PSet( 
            TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
            minRing = cms.int32( 1 ),
            skipClusters = cms.InputTag( "hltDisplacedhltIter4ClustersRefRemovalForTau" ),
            useRingSlector = cms.bool( True ),
            clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
            maxRing = cms.int32( 2 )
        ),
        MTIB = cms.PSet( 
            TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
            skipClusters = cms.InputTag( "hltDisplacedhltIter4ClustersRefRemovalForTau" ),
            clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) )
        ),
        MTID = cms.PSet( 
            TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
            minRing = cms.int32( 3 ),
            skipClusters = cms.InputTag( "hltDisplacedhltIter4ClustersRefRemovalForTau" ),
            useRingSlector = cms.bool( True ),
            clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
            maxRing = cms.int32( 3 )
        ),
        MTOB = cms.PSet(  ),
        MTEC = cms.PSet( 
            TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
            minRing = cms.int32( 3 ),
            skipClusters = cms.InputTag( "hltDisplacedhltIter4ClustersRefRemovalForTau" ),
            useRingSlector = cms.bool( True ),
            clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
            maxRing = cms.int32( 3 )
        )
    )
    process.hltDisplacedhltIter4PFlowPixelLessTrackingRegionsForTau = cms.EDProducer( 
        "CandidateSeededTrackingRegionsEDProducer",
        RegionPSet = cms.PSet( 
            vertexCollection = cms.InputTag( "hltTrimmedPixelVertices" ),
            zErrorVetex = cms.double( 12.0 ),
            beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
            zErrorBeamSpot = cms.double( 15.0 ),
            maxNVertices = cms.int32( 10 ),
            maxNRegions = cms.int32( 100 ),
            nSigmaZVertex = cms.double( 3.0 ),
            nSigmaZBeamSpot = cms.double( 3.0 ),
            ptMin = cms.double( 0.8 ),
            mode = cms.string( "VerticesFixed" ),
            input = cms.InputTag( "hltL2TausForPixelIsolationL1TauSeeded" ),
            searchOpt = cms.bool( True ),
            whereToUseMeasurementTracker = cms.string( "ForSiStrips" ),
            originRadius = cms.double( 1.0 ),
            measurementTrackerName = cms.InputTag( "hltDisplacedhltIter4MaskedMeasurementTrackerEventForTau" ),
            precise = cms.bool( True ),
            deltaEta = cms.double( 0.5 ),
            deltaPhi = cms.double( 0.5 )
        )
    )
    process.hltDisplacedhltIter4PFlowPixelLessClusterCheckForTau = cms.EDProducer( 
        "ClusterCheckerEDProducer",
        doClusterCheck = cms.bool( False ),
        MaxNumberOfCosmicClusters = cms.uint32( 800000 ),
        ClusterCollectionLabel = cms.InputTag( "hltSiStripClusters" ),
        MaxNumberOfPixelClusters = cms.uint32( 40000 ),
        PixelClusterCollectionLabel = cms.InputTag( "hltSiPixelClusters" ),
        cut = cms.string( "" ),
        silentClusterCheck = cms.untracked.bool( False )
    )
    process.hltDisplacedhltIter4PFlowPixelLessHitDoubletsForTau = cms.EDProducer( 
        "HitPairEDProducer",
        seedingLayers = cms.InputTag( "hltDisplacedhltIter4PixelLessLayerTripletsForTau" ),
        trackingRegions = cms.InputTag( "hltDisplacedhltIter4PFlowPixelLessTrackingRegionsForTau" ),
        trackingRegionsSeedingLayers = cms.InputTag( "" ),
        clusterCheck = cms.InputTag( "hltDisplacedhltIter4PFlowPixelLessClusterCheckForTau" ),
        produceSeedingHitSets = cms.bool( False ),
        produceIntermediateHitDoublets = cms.bool( True ),
        maxElement = cms.uint32( 0 ),
        maxElementTotal = cms.uint32( 50000000 ),
        layerPairs = cms.vuint32( 0 )
    )
    process.hltDisplacedhltIter4PFlowPixelLessHitTripletsForTau = cms.EDProducer( 
        "MultiHitFromChi2EDProducer",
        doublets = cms.InputTag( "hltDisplacedhltIter4PFlowPixelLessHitDoubletsForTau" ),
        maxElement = cms.uint32( 100000 ),
        useFixedPreFiltering = cms.bool( False ),
        phiPreFiltering = cms.double( 0.3 ),
        extraHitRPhitolerance = cms.double( 0.0 ),
        extraHitRZtolerance = cms.double( 0.0 ),
        extraZKDBox = cms.double( 0.2 ),
        extraRKDBox = cms.double( 0.2 ),
        extraPhiKDBox = cms.double( 0.005 ),
        fnSigmaRZ = cms.double( 2.0 ),
        refitHits = cms.bool( True ),
        ClusterShapeHitFilterName = cms.string( "ClusterShapeHitFilter" ),
        TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        maxChi2 = cms.double( 5.0 ),
        chi2VsPtCut = cms.bool( True ),
        pt_interv = cms.vdouble( 0.4, 0.7, 1.0, 2.0 ),
        chi2_cuts = cms.vdouble( 3.0, 4.0, 5.0, 5.0 ),
        detIdsToDebug = cms.vint32( 0, 0, 0 )
    )
    process.hltDisplacedhltIter4PFlowPixelLessSeedsForTau = cms.EDProducer( 
        "SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
        seedingHitSets = cms.InputTag( "hltDisplacedhltIter4PFlowPixelLessHitTripletsForTau" ),
        propagator = cms.string( "PropagatorWithMaterialParabolicMf" ),
        SeedMomentumForBOFF = cms.double( 5.0 ),
        OriginTransverseErrorMultiplier = cms.double( 1.0 ),
        MinOneOverPtError = cms.double( 1.0 ),
        TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        magneticField = cms.string( "ParabolicMf" ),
        forceKinematicWithRegionDirection = cms.bool( False ),
        SeedComparitorPSet = cms.PSet( 
            FilterStripHits = cms.bool( False ),
            FilterPixelHits = cms.bool( False ),
            ClusterShapeHitFilterName = cms.string( "ClusterShapeHitFilter" ),
            FilterAtHelixStage = cms.bool( True ),
            ComponentName = cms.string( "PixelClusterShapeSeedComparitor" )
        )
    )
    process.hltDisplacedhltIter4PFlowCkfTrackCandidatesForTau = cms.EDProducer( 
        "CkfTrackCandidateMaker",
        RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
        TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
        cleanTrajectoryAfterInOut = cms.bool( False ),
        reverseTrajectories = cms.bool( False ),
        useHitsSplitting = cms.bool( False ),
        doSeedingRegionRebuilding = cms.bool( False ),
        maxNSeeds = cms.uint32( 100000 ),
        maxSeedsBeforeCleaning = cms.uint32( 1000 ),
        src = cms.InputTag( "hltDisplacedhltIter4PFlowPixelLessSeedsForTau" ),
        SimpleMagneticField = cms.string( "ParabolicMf" ),
        NavigationSchool = cms.string( "SimpleNavigationSchool" ),
        TrajectoryBuilder = cms.string( "" ),
        TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter4PSetTrajectoryBuilderIT" ) ),
        TransientInitialStateEstimatorParameters = cms.PSet( 
            propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
            numberMeasurementsForFit = cms.int32( 4 ),
            propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
        ),
        MeasurementTrackerEvent = cms.InputTag( "hltDisplacedhltIter4MaskedMeasurementTrackerEventForTau" )
    )
    process.hltDisplacedhltIter4PFlowCtfWithMaterialTracksForTau = cms.EDProducer( 
        "TrackProducer",
        useSimpleMF = cms.bool( True ),
        SimpleMagneticField = cms.string( "ParabolicMf" ),
        src = cms.InputTag( "hltDisplacedhltIter4PFlowCkfTrackCandidatesForTau" ),
        clusterRemovalInfo = cms.InputTag( "" ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        Fitter = cms.string( "hltESPFittingSmootherIT" ),
        useHitsSplitting = cms.bool( False ),
        alias = cms.untracked.string( "ctfWithMaterialTracks" ),
        TrajectoryInEvent = cms.bool( False ),
        TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        AlgorithmName = cms.string( "hltIterX" ),
        Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
        GeometricInnerState = cms.bool( True ),
        NavigationSchool = cms.string( "" ),
        MeasurementTracker = cms.string( "" ),
        MeasurementTrackerEvent = cms.InputTag( "hltDisplacedhltIter4MaskedMeasurementTrackerEventForTau" )
    )
    process.hltDisplacedhltIter4PFlowTrackSelectionHighPurityForTau = cms.EDProducer( 
        "AnalyticalTrackSelector",
        src = cms.InputTag( "hltDisplacedhltIter4PFlowCtfWithMaterialTracksForTau" ),
        keepAllTracks = cms.bool( False ),
        beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
        useVertices = cms.bool( True ),
        useVtxError = cms.bool( False ),
        vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
        vtxNumber = cms.int32( -1 ),
        vertexCut = cms.string( "tracksSize>=3" ),
        copyTrajectories = cms.untracked.bool( False ),
        copyExtras = cms.untracked.bool( True ),
        qualityBit = cms.string( "highPurity" ),
        chi2n_par = cms.double( 0.25 ),
        chi2n_no1Dmod_par = cms.double( 9999.0 ),
        res_par = cms.vdouble( 0.003, 0.001 ),
        d0_par1 = cms.vdouble( 1.0, 4.0 ),
        dz_par1 = cms.vdouble( 1.0, 4.0 ),
        d0_par2 = cms.vdouble( 1.0, 4.0 ),
        dz_par2 = cms.vdouble( 1.0, 4.0 ),
        applyAdaptedPVCuts = cms.bool( True ),
        max_d0 = cms.double( 100.0 ),
        max_z0 = cms.double( 100.0 ),
        nSigmaZ = cms.double( 3.0 ),
        minNumberLayers = cms.uint32( 5 ),
        minNumber3DLayers = cms.uint32( 0 ),
        maxNumberLostLayers = cms.uint32( 0 ),
        minHitsToBypassChecks = cms.uint32( 20 ),
        applyAbsCutsIfNoPV = cms.bool( False ),
        max_d0NoPV = cms.double( 100.0 ),
        max_z0NoPV = cms.double( 100.0 ),
        max_relpterr = cms.double( 9999.0 ),
        min_nhits = cms.uint32( 0 ),
        max_minMissHitOutOrIn = cms.int32( 99 ),
        max_lostHitFraction = cms.double( 1.0 ),
        max_eta = cms.double( 9999.0 ),
        min_eta = cms.double( -9999.0 )
    )
    process.hltIter4MergedWithIter0ForTau = cms.EDProducer( 
        "TrackListMerger",
        ShareFrac = cms.double( 0.19 ),
        FoundHitBonus = cms.double( 5.0 ),
        LostHitPenalty = cms.double( 20.0 ),
        MinPT = cms.double( 0.05 ),
        Epsilon = cms.double( -0.001 ),
        MaxNormalizedChisq = cms.double( 1000.0 ),
        MinFound = cms.int32( 3 ),
        TrackProducers = cms.VInputTag( 'hltMergedTracks','hltDisplacedhltIter4PFlowTrackSelectionHighPurityForTau' ),
        hasSelector = cms.vint32( 0, 0 ),
        indivShareFrac = cms.vdouble( 1.0, 1.0 ),
        selectedTrackQuals = cms.VInputTag( 'hltMergedTracks','hltDisplacedhltIter4PFlowTrackSelectionHighPurityForTau' ),
        setsToMerge = cms.VPSet( 
            cms.PSet(  pQual = cms.bool( False ),
                       tLists = cms.vint32( 0, 1 )
                   )
        ),
        trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
        allowFirstHitShare = cms.bool( True ),
        newQuality = cms.string( "confirmed" ),
        copyExtras = cms.untracked.bool( True ),
        writeOnlyTrkQuals = cms.bool( False ),
        copyMVA = cms.bool( False )
    )
    ### End of modules for iter-4 for taus
    
    ### Use new collection of merged tracks including iter-4 for taus to feed hltPFMuonMerging (temporary patch)
    if hasattr(process,'hltPFMuonMerging'):
        process.hltPFMuonMerging.TrackProducers     = cms.VInputTag( 'hltIterL3MuonTracks','hltIter4MergedWithIter0ForTau' )
        process.hltPFMuonMerging.selectedTrackQuals = cms.VInputTag( 'hltIterL3MuonTracks','hltIter4MergedWithIter0ForTau' )
    
    ### Define sequence for iter-4 for taus
    process.HLTIterativeTrackingIteration4ForTau = cms.Sequence( process.hltDisplacedhltIter4ClustersRefRemovalForTau + process.hltDisplacedhltIter4MaskedMeasurementTrackerEventForTau + process.hltDisplacedhltIter4PixelLessLayerTripletsForTau + process.hltDisplacedhltIter4PFlowPixelLessTrackingRegionsForTau + process.hltDisplacedhltIter4PFlowPixelLessClusterCheckForTau + process.hltDisplacedhltIter4PFlowPixelLessHitDoubletsForTau + process.hltDisplacedhltIter4PFlowPixelLessHitTripletsForTau + process.hltDisplacedhltIter4PFlowPixelLessSeedsForTau + process.hltDisplacedhltIter4PFlowCkfTrackCandidatesForTau + process.hltDisplacedhltIter4PFlowCtfWithMaterialTracksForTau + process.hltDisplacedhltIter4PFlowTrackSelectionHighPurityForTau + process.hltIter4MergedWithIter0ForTau )
    
    ### Add sequence for iter-4 for taus to HLTTrackReconstructionForPF
    if hasattr(process,'HLTTrackReconstructionForPF'):
        process.HLTTrackReconstructionForPF = cms.Sequence( 
            process.HLTDoLocalPixelSequence + 
            process.HLTRecopixelvertexingSequence + 
            process.HLTDoLocalStripSequence + 
            process.HLTIterativeTrackingIter02 + 
            ### Add displaced iteration
            process.HLTIterativeTrackingIteration4ForTau + 
            ###
            process.hltPFMuonMerging + 
            process.hltMuonLinks + 
            process.hltMuons 
        )
    
    ### End of iter-4 for taus customization

    return process

