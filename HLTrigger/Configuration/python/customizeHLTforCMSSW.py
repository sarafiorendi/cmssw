import FWCore.ParameterSet.Config as cms

#
# reusable functions
def producers_by_type(process, *types):
    return (module for module in process._Process__producers.values() if module._TypedParameterizable__type in types)
def filters_by_type(process, *types):
    return (filter for filter in process._Process__filters.values() if filter._TypedParameterizable__type in types)
def analyzers_by_type(process, *types):
    return (analyzer for analyzer in process._Process__analyzers.values() if analyzer._TypedParameterizable__type in types)

def esproducers_by_type(process, *types):
    return (module for module in process._Process__esproducers.values() if module._TypedParameterizable__type in types)

#
# one action function per PR - put the PR number into the name of the function

# example:
# def customiseFor12718(process):
#     for pset in process._Process__psets.values():
#         if hasattr(pset,'ComponentType'):
#             if (pset.ComponentType == 'CkfBaseTrajectoryFilter'):
#                 if not hasattr(pset,'minGoodStripCharge'):
#                     pset.minGoodStripCharge = cms.PSet(refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone'))
#     return process

#
# CMSSW version specific customizations
def customizeHLTforCMSSW(process, menuType="GRun"):

    import os
    cmsswVersion = os.environ['CMSSW_VERSION']

    if cmsswVersion >= "CMSSW_8_0":
#       process = customiseFor12718(process)
        pass

#   stage-2 changes if needed
    if ( (menuType == "Fake") or (menuType == "FULL") ):
        return process

#   set stage-2 L1T menu (as long as o2o is missing)
    if (menuType == "HIon"):
        from HLTrigger.Configuration.CustomConfigs import L1THIon
        process = L1THIon(process)
    elif (menuType == "PRef"):
        from HLTrigger.Configuration.CustomConfigs import L1TPRef
        process = L1TPRef(process)

#   remove modules/sequences/paths with legacy/stage-1 dependence

    badTypes = (
                'HLTEgammaL1MatchFilterRegional',
                'HLTEcalRecHitInAllL1RegionsProducer',
                'EgammaHLTCaloTowerProducer',
                'HLTMuonL1Filter',
                'HLTL1MuonSelector',
                'L2MuonSeedGenerator',
                'CaloTowerCreatorForTauHLT',
                'HLTPFJetL1MatchProducer',
                'IsolatedPixelTrackCandidateProducer',
#
                'HLTLevel1Activity',
#
                )
    print " "
    print "# Bad module types: ",badTypes

    badModules = [ ]
    for badType in badTypes:
        print " "
        print '## Bad module type: ',badType
        for module in producers_by_type(process,badType):
            label = module._Labelable__label
            badModules += [label]
            print '### producer label: ',label
        for module in filters_by_type(process,badType):
            label = module._Labelable__label
            badModules += [label]
            print '### filter   label: ',label
        for module in analyzers_by_type(process,badType):
            label = module._Labelable__label
            badModules += [label]
            print '### analyzer label: ',label

    print " "
    for label in badModules:
        setattr(process,label,cms.EDFilter("HLTBool",result=cms.bool(False)))

    return process
