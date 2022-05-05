///## NEW release
#include "L1Trigger/TrackFindingTracklet/interface/PurgeDuplicate.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"
#include "L1Trigger/TrackFindingTracklet/interface/CleanTrackMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Track.h"

#ifdef USEHYBRID
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"
#include "L1Trigger/TrackFindingTMTT/interface/KFParamsComb.h"
#include "L1Trigger/TrackFindingTracklet/interface/HybridFit.h"

#include "L1Trigger/TrackFindingTracklet/interface/TrackletCalculatorDisplaced.h"
#endif

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <unordered_set>
#include <algorithm>

using namespace std;
using namespace trklet;

PurgeDuplicate::PurgeDuplicate(std::string name, Settings const& settings, Globals* global)
    : ProcessBase(name, settings, global) {}

void PurgeDuplicate::addOutput(MemoryBase* memory, std::string output) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding output to " << memory->getName() << " to output "
                                 << output;
  }
  unordered_set<string> outputs = {"trackout",
                                   "trackout1",
                                   "trackout2",
                                   "trackout3",
                                   "trackout4",
                                   "trackout5",
                                   "trackout6",
                                   "trackout7",
                                   "trackout8",
                                   "trackout9",
                                   "trackout10",
                                   "trackout11"};
  if (outputs.find(output) != outputs.end()) {
    auto* tmp = dynamic_cast<CleanTrackMemory*>(memory);
    assert(tmp != nullptr);
    outputtracklets_.push_back(tmp);
    return;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " could not find output: " << output;
}

void PurgeDuplicate::addInput(MemoryBase* memory, std::string input) {
  if (settings_.writetrace()) {
    edm::LogVerbatim("Tracklet") << "In " << name_ << " adding input from " << memory->getName() << " to input "
                                 << input;
  }
  unordered_set<string> inputs = {"trackin",
                                  "trackin1",
                                  "trackin2",
                                  "trackin3",
                                  "trackin4",
                                  "trackin5",
                                  "trackin6",
                                  "trackin7",
                                  "trackin8",
                                  "trackin9",
                                  "trackin10",
                                  "trackin11",
                                  "trackin12"};
  if (inputs.find(input) != inputs.end()) {
    auto* tmp = dynamic_cast<TrackFitMemory*>(memory);
    assert(tmp != nullptr);
    inputtrackfits_.push_back(tmp);
    return;
  }
  throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__ << " could not find input: " << input;
}

void PurgeDuplicate::execute(std::vector<Track>& outputtracks_, unsigned int iSector) {
  inputtracklets_.clear();
  inputtracks_.clear();

  inputstubidslists_.clear();
  inputstublists_.clear();
  myinputstublists_.clear();
  mergedstubidslists_.clear();

  if (settings_.removalType() != "merge") {
    for (auto& inputtrackfit : inputtrackfits_) {
      if (inputtrackfit->nTracks() == 0)
        continue;
      for (unsigned int j = 0; j < inputtrackfit->nTracks(); j++) {
        Track* aTrack = inputtrackfit->getTrack(j)->getTrack();
        aTrack->setSector(iSector);
        inputtracks_.push_back(aTrack);
      }
    }
    if (inputtracks_.empty())
      return;
  }

  unsigned int numTrk = inputtracks_.size();

  ////////////////////
  // Hybrid Removal //
  ////////////////////
#ifdef USEHYBRID

  if (settings_.removalType() == "merge") {
    std::vector<std::pair<int, bool>> trackInfo;  // Track seed & duplicate flag
    // Vector to store the relative rank of the track candidate for merging, based on seed type
    std::vector<int> seedRank;

    // Get vectors from TrackFit and save them
    // inputtracklets: Tracklet objects from the FitTrack (not actually fit yet)
    // inputstublists: L1Stubs for that track
    // inputstubidslists: Stub stubIDs for that 3rack
    // mergedstubidslists: the same as inputstubidslists, but will be used during duplicate removal
    for (unsigned int i = 0; i < inputtrackfits_.size(); i++) {
      if (inputtrackfits_[i]->nStublists() == 0)
        continue;
      if (inputtrackfits_[i]->nStublists() != inputtrackfits_[i]->nTracks())
        throw "Number of stublists and tracks don't match up!";

      int countSeedsPreMerge = 0; 
      for (unsigned int j = 0; j < inputtrackfits_[i]->nStublists(); j++) {
        Tracklet* aTrack = inputtrackfits_[i]->getTrack(j);
        inputtracklets_.push_back(inputtrackfits_[i]->getTrack(j));

        std::vector<const Stub*> stublist = inputtrackfits_[i]->getStublist(j);

        countSeedsPreMerge = 0;
        for (auto& st : stublist){
          int stubLayer = (findLayerDisk(st)).first;
          int stubDisk  = (findLayerDisk(st)).second;
          if ( isSeedingStub(aTrack->seedIndex(), stubLayer, stubDisk)) countSeedsPreMerge++;    
        }
        inputstublists_.push_back(stublist);
        myinputstublists_.push_back(stublist);

        std::vector<std::pair<int, int>> stubidslist = inputtrackfits_[i]->getStubidslist(j);
        inputstubidslists_.push_back(stubidslist);
        mergedstubidslists_.push_back(stubidslist);

        // Encoding: L1L2=0, L2L3=1, L3L4=2, L5L6=3, D1D2=4, D3D4=5, L1D1=6, L2D1=7
        // Best Guess:          L1L2 > L1D1 > L2L3 > L2D1 > D1D2 > L3L4 > L5L6 > D3D4
        // Best Rank:           L1L2 > L3L4 > D3D4 > D1D2 > L2L3 > L2D1 > L5L6 > L1D1
        // Rank-Informed Guess: L1L2 > L3L4 > L1D1 > L2L3 > L2D1 > D1D2 > L5L6 > D3D4
        unsigned int curSeed = aTrack->seedIndex();
        if (curSeed == 0) {
          seedRank.push_back(1);
        } else if (curSeed == 2) {
          seedRank.push_back(2);
        } else if (curSeed == 5) {
          seedRank.push_back(3);
        } else if (curSeed == 4) {
          seedRank.push_back(4);
        } else if (curSeed == 1) {
          seedRank.push_back(5);
        } else if (curSeed == 7) {
          seedRank.push_back(6);
        } else if (curSeed == 3) {
          seedRank.push_back(7);
        } else if (curSeed == 6) {
          seedRank.push_back(8);
        } else if (settings_.extended()) {
          seedRank.push_back(9);
        } else {
          throw cms::Exception("LogError") << __FILE__ << " " << __LINE__ << " Seed " << curSeed
                                           << " not found in list, and settings->extended() not set.";
        }

        if (stublist.size() != stubidslist.size())
          throw "Number of stubs and stubids don't match up!";

        trackInfo.emplace_back(i, false);
      }
    }

    if (inputtracklets_.empty())
      return;
    unsigned int numStublists = inputstublists_.size();

    // Initialize all-false 2D array of tracks being duplicates to other tracks
    bool dupMap[numStublists][numStublists];  // Ends up symmetric
    for (unsigned int itrk = 0; itrk < numStublists; itrk++) {
      for (unsigned int jtrk = 0; jtrk < numStublists; jtrk++) {
        dupMap[itrk][jtrk] = false;
      }
    }

    // Find duplicates; Fill dupMap by looping over all pairs of "tracks"
    // numStublists-1 since last track has no other to compare to
    for (unsigned int itrk = 0; itrk < numStublists - 1; itrk++) {
      for (unsigned int jtrk = itrk + 1; jtrk < numStublists; jtrk++) {
        // Get primary track stubids
        const std::vector<std::pair<int, int>>& stubsTrk1 = inputstubidslists_[itrk];

        // Get and count secondary track stubids
        const std::vector<std::pair<int, int>>& stubsTrk2 = inputstubidslists_[jtrk];

        // Count number of Unique Regions (UR) that share stubs, and the number of UR that each track hits
        unsigned int nShareUR = 0;
        unsigned int nURStubTrk1 = 0;
        unsigned int nURStubTrk2 = 0;
        if (settings_.mergeComparison() == "CompareAll") {
          bool URArray[16];
          for (auto& i : URArray) {
            i = false;
          };
          for (const auto& st1 : stubsTrk1) {
            for (const auto& st2 : stubsTrk2) {
              if (st1.first == st2.first && st1.second == st2.second) {
                // Converts region encoded in st1->first to an index in the Unique Region (UR) array
                int i = st1.first;
                int reg = (i > 0 && i < 10) * (i - 1) + (i > 10) * (i - 5) - (i < 0) * i;
                if (!URArray[reg]) {
                  nShareUR++;
                  URArray[reg] = true;
                }
              }
            }
          }
        } else if (settings_.mergeComparison() == "CompareBest") {
          std::vector<const Stub*> fullStubslistsTrk1 = inputstublists_[itrk];
          std::vector<const Stub*> fullStubslistsTrk2 = inputstublists_[jtrk];

          // Arrays to store the index of the best stub in each region
          int URStubidsTrk1[16];
          int URStubidsTrk2[16];
          for (int i = 0; i < 16; i++) {
            URStubidsTrk1[i] = -1;
            URStubidsTrk2[i] = -1;
          }
          // For each stub on the first track, find the stub with the best residual and store its index in the URStubidsTrk1 array
          for (unsigned int stcount = 0; stcount < stubsTrk1.size(); stcount++) {
            int i = stubsTrk1[stcount].first;
            int reg = (i > 0 && i < 10) * (i - 1) + (i > 10) * (i - 5) - (i < 0) * i;
            double nres = getPhiRes(inputtracklets_[itrk], fullStubslistsTrk1[stcount]);
            double ores = 0;
            if (URStubidsTrk1[reg] != -1)
              ores = getPhiRes(inputtracklets_[itrk], fullStubslistsTrk1[URStubidsTrk1[reg]]);
            if (URStubidsTrk1[reg] == -1 || nres < ores) {
              URStubidsTrk1[reg] = stcount;
            }
          }
          // For each stub on the second track, find the stub with the best residual and store its index in the URStubidsTrk2 array
          for (unsigned int stcount = 0; stcount < stubsTrk2.size(); stcount++) {
            int i = stubsTrk2[stcount].first;
            int reg = (i > 0 && i < 10) * (i - 1) + (i > 10) * (i - 5) - (i < 0) * i;
            double nres = getPhiRes(inputtracklets_[jtrk], fullStubslistsTrk2[stcount]);
            double ores = 0;
            if (URStubidsTrk2[reg] != -1)
              ores = getPhiRes(inputtracklets_[jtrk], fullStubslistsTrk2[URStubidsTrk2[reg]]);
            if (URStubidsTrk2[reg] == -1 || nres < ores) {
              URStubidsTrk2[reg] = stcount;
            }
          }
          // For all 16 regions (6 layers and 10 disks), count the number of regions who's best stub on both tracks are the same
          for (int i = 0; i < 16; i++) {
            int t1i = URStubidsTrk1[i];
            int t2i = URStubidsTrk2[i];
            if (t1i != -1 && t2i != -1 && stubsTrk1[t1i].first == stubsTrk2[t2i].first &&
                stubsTrk1[t1i].second == stubsTrk2[t2i].second)
              nShareUR++;
          }
          // Calculate the number of unique regions hit by each track, so that this number can be used in calculating the number of independent
          // stubs on a track (not enabled/used by default)
          for (int i = 0; i < 16; i++) {
            if (URStubidsTrk1[i] != -1)
              nURStubTrk1++;
            if (URStubidsTrk2[i] != -1)
              nURStubTrk2++;
          }
        }

        // Fill duplicate map
        if (nShareUR >= settings_.minIndStubs()) {  // For number of shared stub merge condition
          dupMap[itrk][jtrk] = true;
          dupMap[jtrk][itrk] = true;
        }
      }
    }

    // invent stub coordinate before the merging
    bool newApproach = false;
    bool newApproachBeforeMerging = true;
    int mergingType = 4;  // 0 = original, 1 = original with inverted order, 
                          // 2 = do not merge seeds, 3 = 1 but set bend to be the same as for preferred track seeding stub 
                          // 4 to be used with newApproachBeforeMerging
    if (newApproachBeforeMerging){
      for (unsigned int itrk = 0; itrk < numStublists; itrk++) {
        Tracklet* tracklet = inputtracklets_[itrk];
        int theSeedIndex = tracklet->seedIndex() ;
      
        std::vector<const Stub*> originalStubsList = inputstublists_[itrk];
        std::vector<const Stub*> newStubList;

        for (unsigned int stubit = 0; stubit < originalStubsList.size(); stubit++) {
          const Stub* thisStub = originalStubsList[stubit];
          if ( isSeedingStub(tracklet->seedIndex(), (findLayerDisk(thisStub)).first, (findLayerDisk(thisStub)).second)) {
            // get a vector containing r, z, phi
            std::vector<double> inv_r_z_phi = get_invented_coords_displ(iSector, thisStub, tracklet );
//             std::vector<double> inv_r_z_phi = get_invented_coords_displ(iSector, thisStub, tracklet );
            
            double stub_x_invent = inv_r_z_phi[0] * std::cos(inv_r_z_phi[2]);
            double stub_y_invent = inv_r_z_phi[0] * std::sin(inv_r_z_phi[2]);
            double stub_z_invent = inv_r_z_phi[1];
            
            Stub* invent_stub_ptr = new Stub(*thisStub) ;

            const L1TStub* L1stub = thisStub->l1tstub();
            L1TStub invent_L1stub ( L1stub->DTClink(),
                                    L1stub->region(),
                                    L1stub->layerdisk(),
                                    L1stub->stubword(),
                                    L1stub->isPSmodule(),
                                    L1stub->isFlipped(),
                                    stub_x_invent,
                                    stub_y_invent,
                                    stub_z_invent,
                                    L1stub->bend(),
                                    L1stub->strip(),
                                    L1stub->tps(),
                                    L1stub->ttStubRef()	
                                  );
                               
            invent_stub_ptr->setl1tstub(new L1TStub(invent_L1stub));
            invent_stub_ptr->l1tstub()->setAllStubIndex(L1stub->allStubIndex());
            invent_stub_ptr->l1tstub()->setUniqueIndex(L1stub->uniqueIndex());
            
            newStubList.push_back(invent_stub_ptr);  

            // to enable comparison output file                  
            std::cout << invent_stub_ptr->isBarrel()<< "\t" << l1tinfo(&invent_L1stub, "invent").c_str() << "\t" <<  l1tinfo(L1stub, "original").c_str() 
                      << "\tdr\t"    << abs(invent_L1stub.r()-L1stub->r())  
                      << "\tdz\t"    << abs(invent_L1stub.z()-L1stub->z()) 
                      << "\tdphi\t"  << abs(invent_L1stub.phi()-L1stub->phi()) 
                      << "\tisPS\t"  << L1stub->isPSmodule() << std::endl;
          }
          else{
            newStubList.push_back(thisStub);
          }
        }
        myinputstublists_[itrk] = newStubList;
      }
      
    }





    // Merge duplicate tracks
    for (unsigned int itrk = 0; itrk < numStublists - 1; itrk++) {
      for (unsigned int jtrk = itrk + 1; jtrk < numStublists; jtrk++) {
        // Merge a track with its first duplicate found.
        if (dupMap[itrk][jtrk]) {
          // Set preferred track based on seed rank
          int preftrk;
          int rejetrk;
          if (seedRank[itrk] < seedRank[jtrk]) {
            preftrk = itrk;
            rejetrk = jtrk;
          } else {
            preftrk = jtrk;
            rejetrk = itrk;
          }

          // Get a merged stub list
          // original merging
          if (mergingType == 0){
            std::vector<const Stub*> newStubList;
            std::vector<const Stub*> stubsTrk1 = inputstublists_[rejetrk];
            std::vector<const Stub*> stubsTrk2 = inputstublists_[preftrk];
            newStubList = stubsTrk1;
            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()) {
                newStubList.push_back(stubsTrk2[stub2it]);
              }
            }
            //   Overwrite stublist of preferred track with merged list
            inputstublists_[preftrk] = newStubList;
  
            std::vector<std::pair<int, int>> newStubidsList;
            std::vector<std::pair<int, int>> stubidsTrk1 = mergedstubidslists_[rejetrk];
            std::vector<std::pair<int, int>> stubidsTrk2 = mergedstubidslists_[preftrk];
            newStubidsList = stubidsTrk1;
  
            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()){
                newStubidsList.push_back(stubidsTrk2[stub2it]);
              }
            }
            // Overwrite stubidslist of preferred track with merged list
            mergedstubidslists_[preftrk] = newStubidsList;
          }
          else if (mergingType == 1){
          
            // reverse order original merging
            std::vector<const Stub*> newStubList;
            std::vector<const Stub*> stubsTrk1 = inputstublists_[preftrk];
            std::vector<const Stub*> stubsTrk2 = inputstublists_[rejetrk];
            newStubList = stubsTrk1;
            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()) {
                newStubList.push_back(stubsTrk2[stub2it]);
              }
            }
            //   Overwrite stublist of preferred track with merged list
            inputstublists_[preftrk] = newStubList;
  
            std::vector<std::pair<int, int>> newStubidsList;
            std::vector<std::pair<int, int>> stubidsTrk1 = mergedstubidslists_[preftrk];
            std::vector<std::pair<int, int>> stubidsTrk2 = mergedstubidslists_[rejetrk];
            newStubidsList = stubidsTrk1;
  
            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()){
                newStubidsList.push_back(stubidsTrk2[stub2it]);
              }
            }
            // Overwrite stubidslist of preferred track with merged list
            mergedstubidslists_[preftrk] = newStubidsList;
          }
          else if (mergingType == 2){

            std::vector<const Stub*> newStubList;
            std::vector<const Stub*> stubsTrk1 = inputstublists_[preftrk];
            std::vector<const Stub*> stubsTrk2 = inputstublists_[rejetrk];
            newStubList = stubsTrk1;
 
            std::vector<std::pair<int, int>> newStubidsList;
            std::vector<std::pair<int, int>> stubidsTrk1 = mergedstubidslists_[rejetrk];
            std::vector<std::pair<int, int>> stubidsTrk2 = mergedstubidslists_[preftrk];
            newStubidsList = stubidsTrk1;
 
            Tracklet* PrefTracklet = inputtracklets_[preftrk];
            Tracklet* RejTracklet  = inputtracklets_[rejetrk];

            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()){
//               std::cout << "---  unique stub  --> "  ;
                if ( isSeedingStub(PrefTracklet->seedIndex(), (findLayerDisk(stubsTrk2[stub2it])).first, (findLayerDisk(stubsTrk2[stub2it])).second)) {
//                   std::cout << "---  not merging this stub  --> " << std::endl; 
//                   std::cout << l1tinfo( stubsTrk2[stub2it]->l1tstub(), "not merged ").c_str() << "\n" << std::endl;// <<  l1tinfo(L1stub, "original").c_str()  << std::endl;
//                  
//                   for (unsigned int seedstub1it = 0; seedstub1it < stubsTrk1.size(); seedstub1it++) {
//                       if ( isSeedingStub(PrefTracklet->seedIndex(), (findLayerDisk(stubsTrk1[seedstub1it])).first, (findLayerDisk(stubsTrk1[seedstub1it])).second)) {
//                           if (stubsTrk1[seedstub1it]->layerdisk() == stubsTrk2[stub2it]->layerdisk()){
//                               std::cout << l1tinfo( stubsTrk1[seedstub1it]->l1tstub(), "pref stub: ").c_str() << "\n" << std::endl;// <<  l1tinfo(L1stub, "original").c_str()  << std::endl;
//                           }
//                       }
//                   }
                  continue;                
                }  
                newStubList.push_back(stubsTrk2[stub2it]);
                newStubidsList.push_back(stubidsTrk2[stub2it]);
              }
            }
            inputstublists_[preftrk] = newStubList;
            mergedstubidslists_[preftrk] = newStubidsList;
          }
          else if (mergingType == 3){ //original merging but change bend of seeding stub to the one of the preferred track
            std::vector<const Stub*> newStubList;
            std::vector<const Stub*> stubsTrk1 = inputstublists_[preftrk];
            std::vector<const Stub*> stubsTrk2 = inputstublists_[rejetrk];
            Tracklet* PrefTracklet = inputtracklets_[preftrk];
            newStubList = stubsTrk1;
            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()) { // if it's new stub from reje track
                
                if ( isSeedingStub(PrefTracklet->seedIndex(), (findLayerDisk(stubsTrk2[stub2it])).first, (findLayerDisk(stubsTrk2[stub2it])).second)){ // and is from seeding layer

                  Stub* clone_stub_ptr = new Stub(*stubsTrk2[stub2it]) ; 
                  const L1TStub* tmp_stub_l1tstub = clone_stub_ptr->l1tstub();
                  for (unsigned int seedstub1it = 0; seedstub1it < stubsTrk1.size(); seedstub1it++) {
                    if ( isSeedingStub(PrefTracklet->seedIndex(), (findLayerDisk(stubsTrk1[seedstub1it])).first, (findLayerDisk(stubsTrk1[seedstub1it])).second) && \
                      stubsTrk1[seedstub1it] -> layerdisk() == stubsTrk2[stub2it]->layerdisk() ) {
                                
                      L1TStub clone_L1stub (    tmp_stub_l1tstub->DTClink(),
                                                tmp_stub_l1tstub->region(),
                                                tmp_stub_l1tstub->layerdisk(),
                                                tmp_stub_l1tstub->stubword(),
                                                tmp_stub_l1tstub->isPSmodule(),
                                                tmp_stub_l1tstub->isFlipped(),
                                                tmp_stub_l1tstub->x(),
                                                tmp_stub_l1tstub->y(),
                                                tmp_stub_l1tstub->z(),
                                                stubsTrk1[seedstub1it]->l1tstub()->bend(),
                                                tmp_stub_l1tstub->strip(),
                                                tmp_stub_l1tstub->tps()	,
                                                tmp_stub_l1tstub->ttStubRef()	
                                              );
                                 
                      clone_stub_ptr->setl1tstub(new L1TStub(clone_L1stub));

                      newStubList.push_back(clone_stub_ptr);
                      break;
                    }
                  } // end loop on pref track stubs
                   
                }
                else {
                  newStubList.push_back(stubsTrk2[stub2it]);
                }
              }
            }
            //   Overwrite stublist of preferred track with merged list
            inputstublists_[preftrk] = newStubList;
  
            std::vector<std::pair<int, int>> newStubidsList;
            std::vector<std::pair<int, int>> stubidsTrk1 = mergedstubidslists_[preftrk];
            std::vector<std::pair<int, int>> stubidsTrk2 = mergedstubidslists_[rejetrk];
            newStubidsList = stubidsTrk1;
  
            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()){
                newStubidsList.push_back(stubidsTrk2[stub2it]);
              }
            }
            // Overwrite stubidslist of preferred track with merged list
            mergedstubidslists_[preftrk] = newStubidsList;
          }
          else if (mergingType == 4){
            std::vector<const Stub*> newStubList;
            std::vector<const Stub*> stubsTrk1 = myinputstublists_[rejetrk];
            std::vector<const Stub*> stubsTrk2 = myinputstublists_[preftrk];
            newStubList = stubsTrk1;
            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()) {
                newStubList.push_back(stubsTrk2[stub2it]);
              }
            }
            //   Overwrite stublist of preferred track with merged list
            myinputstublists_[preftrk] = newStubList;
  
            std::vector<std::pair<int, int>> newStubidsList;
            std::vector<std::pair<int, int>> stubidsTrk1 = mergedstubidslists_[rejetrk];
            std::vector<std::pair<int, int>> stubidsTrk2 = mergedstubidslists_[preftrk];
            newStubidsList = stubidsTrk1;
  
            for (unsigned int stub2it = 0; stub2it < stubsTrk2.size(); stub2it++) {
              if (find(stubsTrk1.begin(), stubsTrk1.end(), stubsTrk2[stub2it]) == stubsTrk1.end()){
                newStubidsList.push_back(stubidsTrk2[stub2it]);
              }
            }
            // Overwrite stubidslist of preferred track with merged list
            mergedstubidslists_[preftrk] = newStubidsList;
          }
          
//           std::cout << "---- preferred track (" << preftrk << ") stubs after this merging:  ";
//           for (auto& st : inputstublists_[preftrk]){
//               int stubLayer = (findLayerDisk(st)).first;
//               int stubDisk  = (findLayerDisk(st)).second;
//               std::cout << stubLayer <<  "/" <<  stubDisk  << "  s? " << isSeedingStub(PrefTracklet->seedIndex(), stubLayer, stubDisk) << "; ";
//           }
//           std::cout << std::endl;


          // Mark that rejected track has been merged into another track
          trackInfo[rejetrk].second = true;
        }
      }
    }

    // Make the final track objects, fit with KF, and send to output
    for (unsigned int itrk = 0; itrk < numStublists; itrk++) {
      bool duplicateTrack = trackInfo[itrk].second;
      if (not duplicateTrack) { // Don't waste CPU by calling KF for duplicates

        Tracklet* tracklet = inputtracklets_[itrk];

        int theSeedIndex = tracklet->seedIndex() ;
        std::vector<const Stub*> alltrackstublist = inputstublists_[itrk];
        if (mergingType == 4)
          alltrackstublist = myinputstublists_[itrk];
        std::vector<const Stub*> trackstublist ;

        // Encoding: L1L2=0, L2L3=1, L3L4=2, L5L6=3, D1D2=4, D3D4=5, L1D1=6, L2D1=7
        if (newApproach){
            std::vector<const Stub*> matchedstublist ;
            std::vector<const Stub*> seedingstublist ;
            std::vector<const Stub*> newseedingstublist ;
            
            // find matched stubs and copy them to the collection of stubs for fitting
            for (auto& st : alltrackstublist){
      
                int stubLayer = (findLayerDisk(st)).first;
                int stubDisk  = (findLayerDisk(st)).second;
                if (!tracklet->match(st->layerdisk()) && isSeedingStub(tracklet->seedIndex(), stubLayer, stubDisk))
                    seedingstublist.push_back(st);
                else{
                    matchedstublist.push_back(st);
                }
            }
            for (auto& st : seedingstublist){
            
                const L1TStub* L1stub = st->l1tstub();
                // return a vector containing r, z, phi
                std::vector<double> inv_r_z_phi = get_invented_coords(iSector, st, tracklet );
                
                double stub_x_invent = inv_r_z_phi[0] * std::cos(inv_r_z_phi[2]);
                double stub_y_invent = inv_r_z_phi[0] * std::sin(inv_r_z_phi[2]);
                double stub_z_invent = inv_r_z_phi[1];

                Stub* invent_stub_ptr = new Stub(*st) ;
                L1TStub invent_L1stub ( L1stub->DTClink(),
                                        L1stub->region(),
                                        L1stub->layerdisk(),
                                        L1stub->stubword(),
                                        L1stub->isPSmodule(),
                                        L1stub->isFlipped(),
                                        stub_x_invent,
                                        stub_y_invent,
                                        stub_z_invent,
                                        L1stub->bend(),
                                        L1stub->strip(),
                                        L1stub->tps(),
                                         L1stub->ttStubRef()
                                      );
                                   
                invent_stub_ptr->setl1tstub(new L1TStub(invent_L1stub));
                invent_stub_ptr->l1tstub() -> setAllStubIndex(L1stub->allStubIndex());
                invent_stub_ptr->l1tstub()->setUniqueIndex(L1stub->uniqueIndex());
                
                newseedingstublist.push_back(invent_stub_ptr);  
                // to enable comparison output file                  
                std::cout << invent_stub_ptr->isBarrel()<< "\t" << l1tinfo(&invent_L1stub, "invent").c_str() << "\t" <<  l1tinfo(L1stub, "original").c_str() 
                        << "\tdr\t"    << abs(invent_L1stub.r()-L1stub->r())  
                        << "\tdz\t"    << abs(invent_L1stub.z()-L1stub->z()) 
                        <<  "\tdphi\t" << abs(invent_L1stub.phi()-L1stub->phi()) 
                        << "\tisPS\t"  << L1stub->isPSmodule() << std::endl;
            }
            // sum the matched and seeding lists
            trackstublist.insert(trackstublist.end(), newseedingstublist.begin(), newseedingstublist.end());
            trackstublist.insert(trackstublist.end(), matchedstublist.begin(), matchedstublist.end());
        }    
        else{
            trackstublist.insert(trackstublist.end(), alltrackstublist.begin(), alltrackstublist.end());
        }
          
        
        
        
//         std::vector<const Stub*> trackstublist = inputstublists_[itrk];

        // Run KF track fit 
        HybridFit hybridFitter(iSector, settings_, globals_);
        hybridFitter.Fit(tracklet, trackstublist);

        // If the track was accepted (and thus fit), add to output
        if (tracklet->fit()) {
          // Add fitted Track to output (later converted to TTTrack)
          Track* outtrack = tracklet->getTrack();
          outtrack->setSector(iSector);
          // Also store fitted track as more detailed Tracklet object.
          outputtracklets_[trackInfo[itrk].first]->addTrack(tracklet);

          // Add all tracks to standalone root file output
          outtrack->setStubIDpremerge(inputstubidslists_[itrk]);
          outtrack->setStubIDprefit(mergedstubidslists_[itrk]);
          outputtracks_.push_back(*outtrack);
        }
      }
    }
  }
#endif

  //////////////////
  // Grid removal //
  //////////////////
  if (settings_.removalType() == "grid") {
    // Sort tracks by ichisq/DoF so that removal will keep the lower ichisq/DoF track
    std::sort(inputtracks_.begin(), inputtracks_.end(), [](const Track* lhs, const Track* rhs) {
      return lhs->ichisq() / lhs->stubID().size() < rhs->ichisq() / rhs->stubID().size();
    });
    bool grid[35][40] = {{false}};

    for (unsigned int itrk = 0; itrk < numTrk; itrk++) {
      if (inputtracks_[itrk]->duplicate())
        edm::LogPrint("Tracklet") << "WARNING: Track already tagged as duplicate!!";

      double phiBin = (inputtracks_[itrk]->phi0(settings_) - 2 * M_PI / 27 * iSector) / (2 * M_PI / 9 / 50) + 9;
      phiBin = std::max(phiBin, 0.);
      phiBin = std::min(phiBin, 34.);

      double ptBin = 1 / inputtracks_[itrk]->pt(settings_) * 40 + 20;
      ptBin = std::max(ptBin, 0.);
      ptBin = std::min(ptBin, 39.);

      if (grid[(int)phiBin][(int)ptBin])
        inputtracks_[itrk]->setDuplicate(true);
      grid[(int)phiBin][(int)ptBin] = true;

      double phiTest = inputtracks_[itrk]->phi0(settings_) - 2 * M_PI / 27 * iSector;
      if (phiTest < -2 * M_PI / 27)
        edm::LogVerbatim("Tracklet") << "track phi too small!";
      if (phiTest > 2 * 2 * M_PI / 27)
        edm::LogVerbatim("Tracklet") << "track phi too big!";
    }
  }  // end grid removal

  //////////////////////////
  // ichi + nstub removal //
  //////////////////////////
  if (settings_.removalType() == "ichi" || settings_.removalType() == "nstub") {
    for (unsigned int itrk = 0; itrk < numTrk - 1; itrk++) {  // numTrk-1 since last track has no other to compare to

      // If primary track is a duplicate, it cannot veto any...move on
      if (inputtracks_[itrk]->duplicate() == 1)
        continue;

      unsigned int nStubP = 0;
      vector<unsigned int> nStubS(numTrk);
      vector<unsigned int> nShare(numTrk);
      // Get and count primary stubs
      std::map<int, int> stubsTrk1 = inputtracks_[itrk]->stubID();
      nStubP = stubsTrk1.size();

      for (unsigned int jtrk = itrk + 1; jtrk < numTrk; jtrk++) {
        // Skip duplicate tracks
        if (inputtracks_[jtrk]->duplicate() == 1)
          continue;

        // Get and count secondary stubs
        std::map<int, int> stubsTrk2 = inputtracks_[jtrk]->stubID();
        nStubS[jtrk] = stubsTrk2.size();

        // Count shared stubs
        for (auto& st : stubsTrk1) {
          if (stubsTrk2.find(st.first) != stubsTrk2.end()) {
            if (st.second == stubsTrk2[st.first])
              nShare[jtrk]++;
          }
        }
      }

      // Tag duplicates
      for (unsigned int jtrk = itrk + 1; jtrk < numTrk; jtrk++) {
        // Skip duplicate tracks
        if (inputtracks_[jtrk]->duplicate() == 1)
          continue;

        // Chi2 duplicate removal
        if (settings_.removalType() == "ichi") {
          if ((nStubP - nShare[jtrk] < settings_.minIndStubs()) ||
              (nStubS[jtrk] - nShare[jtrk] < settings_.minIndStubs())) {
            if ((int)inputtracks_[itrk]->ichisq() / (2 * inputtracks_[itrk]->stubID().size() - 4) >
                (int)inputtracks_[jtrk]->ichisq() / (2 * inputtracks_[itrk]->stubID().size() - 4)) {
              inputtracks_[itrk]->setDuplicate(true);
            } else if ((int)inputtracks_[itrk]->ichisq() / (2 * inputtracks_[itrk]->stubID().size() - 4) <=
                       (int)inputtracks_[jtrk]->ichisq() / (2 * inputtracks_[itrk]->stubID().size() - 4)) {
              inputtracks_[jtrk]->setDuplicate(true);
            } else {
              edm::LogVerbatim("Tracklet") << "Error: Didn't tag either track in duplicate pair.";
            }
          }
        }  // end ichi removal

        // nStub duplicate removal
        if (settings_.removalType() == "nstub") {
          if ((nStubP - nShare[jtrk] < settings_.minIndStubs()) && (nStubP < nStubS[jtrk])) {
            inputtracks_[itrk]->setDuplicate(true);
          } else if ((nStubS[jtrk] - nShare[jtrk] < settings_.minIndStubs()) && (nStubS[jtrk] <= nStubP)) {
            inputtracks_[jtrk]->setDuplicate(true);
          } else {
            edm::LogVerbatim("Tracklet") << "Error: Didn't tag either track in duplicate pair.";
          }
        }  // end nstub removal

      }  // end tag duplicates

    }  // end loop over primary track

  }  // end ichi + nstub removal

  //Add tracks to output
  if (settings_.removalType() != "merge") {
    for (unsigned int i = 0; i < inputtrackfits_.size(); i++) {
      for (unsigned int j = 0; j < inputtrackfits_[i]->nTracks(); j++) {
        if (inputtrackfits_[i]->getTrack(j)->getTrack()->duplicate() == 0) {
          if (settings_.writeMonitorData("Seeds")) {
            ofstream fout("seeds.txt", ofstream::app);
            fout << __FILE__ << ":" << __LINE__ << " " << name_ << "_" << iSector << " "
                 << inputtrackfits_[i]->getTrack(j)->getISeed() << endl;
            fout.close();
          }
          outputtracklets_[i]->addTrack(inputtrackfits_[i]->getTrack(j));
        }
        //For root file:
        outputtracks_.push_back(*inputtrackfits_[i]->getTrack(j)->getTrack());
      }
    }
  }
}

double PurgeDuplicate::getPhiRes(Tracklet* curTracklet, const Stub* curStub) {
  double phiproj;
  double stubphi;
  double phires;
  // Get phi position of stub
  stubphi = curStub->l1tstub()->phi();
  // Get region that the stub is in (Layer 1->6, Disk 1->5)
  int Layer = curStub->layerdisk() + 1;
  if (Layer > N_LAYER) {
    Layer = 0;
  }
  int Disk = curStub->layerdisk() - (N_LAYER - 1);
  if (Disk < 0) {
    Disk = 0;
  }
  // Get phi projection of tracklet
  int seedindex = curTracklet->seedIndex();
  // If this stub is a seed stub, set projection=phi, so that res=0
  if ((seedindex == 0 && (Layer == 1 || Layer == 2)) || (seedindex == 1 && (Layer == 2 || Layer == 3)) ||
      (seedindex == 2 && (Layer == 3 || Layer == 4)) || (seedindex == 3 && (Layer == 5 || Layer == 6)) ||
      (seedindex == 4 && (abs(Disk) == 1 || abs(Disk) == 2)) ||
      (seedindex == 5 && (abs(Disk) == 3 || abs(Disk) == 4)) || (seedindex == 6 && (Layer == 1 || abs(Disk) == 1)) ||
      (seedindex == 7 && (Layer == 2 || abs(Disk) == 1)) ||
      (seedindex == 8 && (Layer == 2 || Layer == 3 || Layer == 4)) ||
      (seedindex == 9 && (Layer == 4 || Layer == 5 || Layer == 6)) ||
      (seedindex == 10 && (Layer == 2 || Layer == 3 || abs(Disk) == 1)) ||
      (seedindex == 11 && (Layer == 2 || abs(Disk) == 1 || abs(Disk) == 2))) {
    phiproj = stubphi;
    // Otherwise, get projection of tracklet
  } else {
    phiproj = curTracklet->proj(curStub->layerdisk()).phiproj();
  }
  // Calculate residual
  phires = std::abs(stubphi - phiproj);
  return phires;
}


// Encoding: L1L2=0, L2L3=1, L3L4=2, L5L6=3, D1D2=4, D3D4=5, L1D1=6, L2D1=7
bool PurgeDuplicate::isSeedingStub(int seedindex, int Layer, int Disk){
  if ((seedindex == 0 && (Layer == 1 || Layer == 2)) || 
      (seedindex == 1 && (Layer == 2 || Layer == 3)) ||
      (seedindex == 2 && (Layer == 3 || Layer == 4)) || 
      (seedindex == 3 && (Layer == 5 || Layer == 6)) ||
      (seedindex == 4 && (abs(Disk) == 1 || abs(Disk) == 2)) ||
      (seedindex == 5 && (abs(Disk) == 3 || abs(Disk) == 4)) || 
      (seedindex == 6 && (Layer == 1 || abs(Disk) == 1)) ||
      (seedindex == 7 && (Layer == 2 || abs(Disk) == 1)) ||
      (seedindex == 8 && (Layer == 2 || Layer == 3 || Layer == 4)) ||
      (seedindex == 9 && (Layer == 4 || Layer == 5 || Layer == 6)) ||
      (seedindex == 10 && (Layer == 2 || Layer == 3 || abs(Disk) == 1)) ||
      (seedindex == 11 && (Layer == 2 || abs(Disk) == 1 || abs(Disk) == 2))) 
      return true;
  
  return false;
}

std::pair<int,int> PurgeDuplicate::findLayerDisk(const Stub* st){

    std::pair<int,int> layer_disk;
    layer_disk.first = st->layerdisk() + 1;
    if (layer_disk.first > N_LAYER) {
        layer_disk.first = 0;
    }
    layer_disk.second = st->layerdisk() - (N_LAYER - 1);
    if (layer_disk.second < 0) {
        layer_disk.second = 0;
    }
    return layer_disk;
}


std::string PurgeDuplicate::l1tinfo(const L1TStub* L1stub, std::string str=""){
    std::string thestr = Form("\t %s stub info:  r/z/phi:\t%f\t%f\t%f\t%d\t%f\t%d", str.c_str(), L1stub->r(), L1stub->z(), L1stub->phi(), L1stub->iphi(), L1stub->bend(),  L1stub->layerdisk()   );
    return thestr;                        
}


std::vector<double> PurgeDuplicate::get_invented_coords(unsigned int iSector, const Stub* st, Tracklet* tracklet){     

  int stubLayer = (findLayerDisk(st)).first;
  int stubDisk = (findLayerDisk(st)).second;
  const L1TStub* L1stub = st->l1tstub();
  
  double stub_phi  = -99;
  double stub_z    = -99;
  double stub_r    = -99;

  double tracklet_rinv = tracklet->rinv();

  if (st->isBarrel()){
      stub_r = settings_.rmean(stubLayer-1);
      stub_phi  = tracklet->phi0() - std::asin(stub_r*tracklet_rinv/2);
      stub_phi  = stub_phi + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG();
      stub_phi  = reco::reduceRange(stub_phi);
      stub_z  = tracklet->z0() + 2 * tracklet->t() * 1 / tracklet_rinv * std::asin(stub_r*tracklet_rinv/2);
  }
  else{
      stub_z = settings_.zmean(stubDisk-1)*tracklet->disk()/abs(tracklet->disk());
      stub_phi  = tracklet->phi0() - (stub_z - tracklet->z0())*tracklet_rinv/2/tracklet->t();
      stub_phi  = stub_phi + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG();
      stub_phi  = reco::reduceRange(stub_phi);
      stub_r  = 2/tracklet_rinv * std::sin ((stub_z - tracklet->z0())*tracklet_rinv/2/tracklet->t());
  }

  std::vector invented_coords{stub_r, stub_z, stub_phi};
  return invented_coords;
}


std::vector<double> PurgeDuplicate::get_invented_coords_displ(unsigned int iSector, const Stub* st, Tracklet* tracklet){     

  int stubLayer = (findLayerDisk(st)).first;
  int stubDisk = (findLayerDisk(st)).second;
  const L1TStub* L1stub = st->l1tstub();
  
  double stub_phi  = -99;
  double stub_z    = -99;
  double stub_r    = -99;
  double r_star    = -99;
  double eps       = -99;
  double coeff_1 = -99.;
  double r_over_2rho = -99.;
  double d0_over_r = -99.;
  double sin_val = -99.;
  double beta = -99.;
  double r_square = -99.;

  double tracklet_rinv = tracklet->rinv();
  double rho = 1/tracklet->rinv();
  double rho_minus_d0 = rho + tracklet->d0(); // should be -, but otherwise does not work
  //   if (rho < 0) rho_minus_d0 = -rho_minus_d0;

  // exact helix
  if (st->isBarrel() && L1stub->isPSmodule()){
//       stub_r = settings_.rmean(2);
      stub_r = settings_.rmean(stubLayer-1);
      
      sin_val = (stub_r*stub_r + rho_minus_d0*rho_minus_d0 - rho*rho) / (2 * stub_r * rho_minus_d0) ;
      stub_phi = tracklet->phi0() - std::asin(sin_val);
      stub_phi  = stub_phi + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG();
      stub_phi  = reco::reduceRange(stub_phi);

      beta =  std::acos((rho*rho + rho_minus_d0*rho_minus_d0 - stub_r*stub_r) / (2 * rho * rho_minus_d0) );
      stub_z = tracklet->z0() + tracklet->t() * std::abs(rho * beta);
//      if (abs(rho)>4000) {
//           std::cout << "rho sara: " << rho << std::endl;
//           std::cout << "projecting to : "  << stub_r << std::endl;
//           std::cout << "r0 sara: "   << rho_minus_d0 << std::endl;
//           std::cout << "d0 sara: "   << tracklet->d0() << std::endl;
//           std::cout << "phi0 sara: " << tracklet->phi0()<< std::endl;
//           std::cout << "phi before: "  << tracklet->phi0() - std::asin(sin_val) << std::endl;
//           std::cout << "phi intermediate: "  << tracklet->phi0() - std::asin(sin_val) + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG()<< std::endl;
//           std::cout << "phi sara: "  << stub_phi<< std::endl;
// //           std::cout << "beta sara: "  << beta<< std::endl;
// //           std::cout << "z0 sara: "  << tracklet->z0() << std::endl;
// //           std::cout << "t sara: "  << tracklet->t() << std::endl;
//       }
  }
  else if (!st->isBarrel() && L1stub->isPSmodule()){
//       stub_z = settings_.zmean(1)*tracklet->disk()/abs(tracklet->disk());
//       rho = std::abs(rho); //not clear from formulas though
//       if (tracklet->t() < 0) stub_z = -stub_z;
      stub_z = settings_.zmean(stubDisk-1)*tracklet->disk()/abs(tracklet->disk());
      beta = (stub_z - tracklet->z0()) / (tracklet->t() * std::abs(rho)); // maybe rho should be abs value
      r_square = -2 * rho * rho_minus_d0 * std::cos(beta) + rho*rho + rho_minus_d0*rho_minus_d0;
      stub_r = sqrt(r_square);
      
      sin_val = (stub_r*stub_r + rho_minus_d0*rho_minus_d0 - rho*rho) / (2 * stub_r * rho_minus_d0) ;
      stub_phi = tracklet->phi0() - std::asin(sin_val);
      stub_phi  = stub_phi + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG();
      stub_phi  = reco::reduceRange(stub_phi);

//      if (abs(rho)>4000) {
//           std::cout << "rho sara: " << rho << std::endl;
//           std::cout << "t sara: " << tracklet->t() << std::endl;
//           std::cout << "projecting to z: "  << stub_z << std::endl;
//           std::cout << "beta sara : "  << stub_z << " * " << tracklet->z0() << " / " << tracklet->t() << " * " << rho << " = " << beta << std::endl;
// //           std::cout << "beta sara: "  << beta << std::endl;
//           std::cout << "d0 sara: "   << tracklet->d0() << std::endl;
// //           std::cout << "r0 sara: "  << rho_minus_d0 << std::endl;
//           std::cout << "r square sara: "  << r_square << std::endl;
//           std::cout << "phi val before "  << tracklet->phi0() - std::asin(sin_val) << std::endl;
// //           std::cout << "r0 sara: "   << rho_minus_d0 << std::endl;
// //           std::cout << "d0 sara: "   << tracklet->d0() << std::endl;
// //           std::cout << "phi0 sara: " << tracklet->phi0()<< std::endl;
// //           std::cout << "phi before: "  << tracklet->phi0() - std::asin(sin_val) << std::endl;
// //           std::cout << "phi intermediate: "  << tracklet->phi0() - std::asin(sin_val) + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG()<< std::endl;
// //           std::cout << "phi sara: "  << stub_phi<< std::endl;
// // //           std::cout << "beta sara: "  << beta<< std::endl;
// // //           std::cout << "z0 sara: "  << tracklet->z0() << std::endl;
// // //           std::cout << "t sara: "  << tracklet->t() << std::endl;
//       }

  }
//       stub_z = settings_.zmean(stubDisk-1)*tracklet->disk()/abs(tracklet->disk());
//       r_star = (stub_z - tracklet->z0()) / tracklet->t();
//       eps    =  pow((tracklet->d0() * tracklet->t()/stub_z), 2) + 1/6*pow( (r_star * tracklet_rinv/2) , 2);
//       stub_r = r_star * (1-eps);
// //      stub_r = L1stub->r();
//      
//       coeff_1  = tracklet->d0()*tracklet->t()/stub_z;
//       stub_phi = tracklet->phi0() - stub_r*tracklet_rinv/2;
//       stub_phi = stub_phi + coeff_1 * (1 + tracklet->z0()/stub_z) * (1 + tracklet->d0()*tracklet_rinv/2) * (1 + eps);
//       stub_phi = stub_phi + 1/6 * pow((-r_star*tracklet_rinv/2 + coeff_1), 3);
//       stub_phi  = stub_phi + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG();
//       stub_phi  = reco::reduceRange(stub_phi);
// 
//   }
//   else{
//       stub_r = L1stub->r();
//       stub_z = L1stub->z();
//       stub_phi = L1stub->phi();
//   }

  // approx helix
//   if (st->isBarrel() && L1stub->isPSmodule()){
//   
//       stub_r = settings_.rmean(stubLayer-1);
// //       double phiproj[1], zproj[1], phider[1], zder[1];
// //       for (unsigned int i = 0; i < 1; i++) {
// //         exactproj(stub_r, tracklet_rinv, tracklet->phi0(), tracklet->d0(), 
// //                                                        tracklet->t(), tracklet->z0(), 
// //                                                        tracklet->d0()+1/tracklet_rinv, 
// //                                                        phiproj[i], zproj[i], phider[i], zder[i]);
// //       }
//       
//       r_over_2rho = stub_r * tracklet_rinv/2;
//       d0_over_r = tracklet->d0() /stub_r;
// //       stub_z = L1stub->z();
//       stub_z  = tracklet->z0() + tracklet->t() * stub_r * (1 + pow(tracklet->d0()/stub_r,2) + 1/6. *(pow(stub_r/2*tracklet_rinv, 2)));
//       stub_phi = tracklet->phi0() - r_over_2rho + d0_over_r + \
//                  d0_over_r * tracklet->d0()*tracklet_rinv/2 - \
//                  r_over_2rho*tracklet->d0() * tracklet_rinv + \
//                  1/6*(pow(-r_over_2rho + d0_over_r, 3));
//       
//       stub_phi  = stub_phi + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG();
//       stub_phi  = reco::reduceRange(stub_phi);
//       
//       if (stub_phi*L1stub->phi() < 0){
//         std::cout << "warning!!! \t r:" << L1stub->r() << " vs " << stub_r << "\t" <<
//                                            L1stub->z() << " vs " << stub_z << "\t" <<
//                                            L1stub->phi() << " vs " << stub_phi << 
//         std::endl;
//       }
// //         std::cout << "compare: \t r:"   << stub_r<< " vs " << stub_r << "\tphi  " <<
// //                                            phiproj[0] << " vs " << stub_phi << "\tz  " <<
// //                                            zproj[0] << " vs " << stub_z << 
// //         std::endl;
//       
//   }
//   else if (!st->isBarrel() && L1stub->isPSmodule()){
//       stub_z = settings_.zmean(stubDisk-1)*tracklet->disk()/abs(tracklet->disk());
//       r_star = (stub_z - tracklet->z0()) / tracklet->t();
//       eps    =  pow((tracklet->d0() * tracklet->t()/stub_z), 2) + 1/6*pow( (r_star * tracklet_rinv/2) , 2);
//       stub_r = r_star * (1-eps);
// //      stub_r = L1stub->r();
//      
//       coeff_1  = tracklet->d0()*tracklet->t()/stub_z;
//       stub_phi = tracklet->phi0() - stub_r*tracklet_rinv/2;
//       stub_phi = stub_phi + coeff_1 * (1 + tracklet->z0()/stub_z) * (1 + tracklet->d0()*tracklet_rinv/2) * (1 + eps);
//       stub_phi = stub_phi + 1/6 * pow((-r_star*tracklet_rinv/2 + coeff_1), 3);
//       stub_phi  = stub_phi + iSector * settings_.dphisector() - 0.5 * settings_.dphisectorHG();
//       stub_phi  = reco::reduceRange(stub_phi);
// 
//   }
  else{
      stub_r = L1stub->r();
      stub_z = L1stub->z();
      stub_phi = L1stub->phi();
  }
  std::vector invented_coords{stub_r, stub_z, stub_phi};
  return invented_coords;
}


void PurgeDuplicate::exactproj(double rproj,
                                            double rinv,
                                            double phi0,
                                            double d0,
                                            double t,
                                            double z0,
                                            double r0,
                                            double& phiproj,
                                            double& zproj,
                                            double& phider,
                                            double& zder) {
  double rho = 1 / rinv;
  if (rho < 0) {
    r0 = -r0;
  }
  phiproj = phi0 - asin((rproj * rproj + r0 * r0 - rho * rho) / (2 * rproj * r0));
  double beta = acos((rho * rho + r0 * r0 - rproj * rproj) / (2 * r0 * rho));
  zproj = z0 + t * std::abs(rho * beta);

  //not exact, but close
  phider = -0.5 * rinv / sqrt(1 - pow(0.5 * rproj * rinv, 2)) + d0 / (rproj * rproj);
  zder = t / sqrt(1 - pow(0.5 * rproj * rinv, 2));

//   if (settings_.debugTracklet())
//     edm::LogVerbatim("Tracklet") << "exact proj layer at " << rproj << " : " << phiproj << " " << zproj;
}
