#ifndef DataFormats_EgammaReco_ElectronSeed_h
#define DataFormats_EgammaReco_ElectronSeed_h

//********************************************************************
//
// A verson of reco::ElectronSeed which can have N hits as part of the 
// 2017 upgrade of E/gamma pixel matching for the phaseI pixels
//
// While it is technically named ElectronSeed, it is effectively a new class
// However to simplify things, the name ElectronSeed was kept
//
// Noticed that h/e values never seem to used anywhere and they are a 
// mild pain to propagate in the new framework so they were removed
//
// infinities are used to mark invalid unset values to maintain 
// compatibilty with the orginal ElectronSeed class
//
// author: S. Harper (RAL), 2017
//
//*********************************************************************


#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectoryState/interface/TrackCharge.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Ref.h"

#include <vector>
#include <limits>

namespace reco
{
  
  class ElectronSeed : public TrajectorySeed {
  public :
    struct PMVars {
      float dRZPos;
      float dRZNeg;
      float dPhiPos;
      float dPhiNeg;
      int detId; //this is already stored as the hit is stored in traj seed but a useful sanity check
      int layerOrDiskNr;//redundant as stored in detId but its a huge pain to hence why its saved here

      PMVars();
      void setDPhi(float pos,float neg){dPhiPos=pos;dPhiNeg=neg;}
      void setDRZ(float pos,float neg){dRZPos=pos;dRZNeg=neg;}
      void setDet(int iDetId,int iLayerOrDiskNr){detId=iDetId;layerOrDiskNr=iLayerOrDiskNr;}

    };
    
    
    typedef edm::OwnVector<TrackingRecHit> RecHitContainer ;
    typedef edm::RefToBase<CaloCluster> CaloClusterRef ;
    typedef edm::Ref<TrackCollection> CtfTrackRef ;
    static std::string const & name()
    {
      static std::string const name_("ElectronSeed") ;
      return name_;
    }
    
    //! Construction of base attributes
    ElectronSeed() ;
    ElectronSeed( const TrajectorySeed & ) ;
    ElectronSeed( PTrajectoryStateOnDet & pts, RecHitContainer & rh,  PropagationDirection & dir ) ;
    ElectronSeed * clone() const { return new ElectronSeed(*this) ; }
    virtual ~ElectronSeed()=default;

    //! Set additional info
    void setCtfTrack( const CtfTrackRef & ) ;
    void setCaloCluster( const CaloClusterRef& clus){caloCluster_=clus;isEcalDriven_=true;}
    void addHitInfo(const PMVars& hitVars){hitInfo_.push_back(hitVars);}
    void setNrLayersAlongTraj(int val){nrLayersAlongTraj_=val;}
    //! Accessors
    const CtfTrackRef& ctfTrack() const { return ctfTrack_ ; }
    const CaloClusterRef& caloCluster() const { return caloCluster_ ; }
   
    //! Utility
    TrackCharge getCharge() const { return startingState().parameters().charge() ; }

    bool isEcalDriven() const { return isEcalDriven_ ; }
    bool isTrackerDriven() const { return isTrackerDriven_ ; }

    const std::vector<PMVars>& hitInfo()const{return hitInfo_;}
    float dPhiNeg(size_t hitNr)const{return getVal(hitNr,&PMVars::dPhiNeg);}
    float dPhiPos(size_t hitNr)const{return getVal(hitNr,&PMVars::dPhiPos);}
    float dPhiBest(size_t hitNr)const{return bestVal(dPhiNeg(hitNr),dPhiPos(hitNr));}
    float dRZPos(size_t hitNr)const{return getVal(hitNr,&PMVars::dRZPos);}
    float dRZNeg(size_t hitNr)const{return getVal(hitNr,&PMVars::dRZNeg);}
    float dRZBest(size_t hitNr)const{return bestVal(dRZNeg(hitNr),dRZPos(hitNr));}
    int detId(size_t hitNr)const{return hitNr<hitInfo_.size() ? hitInfo_[hitNr].detId : 0;}
    int subDet(size_t hitNr)const{return DetId(detId(hitNr)).subdetId();}
    int layerOrDiskNr(size_t hitNr)const{return getVal(hitNr,&PMVars::layerOrDiskNr);}
    int nrLayersAlongTraj()const{return nrLayersAlongTraj_;}

    //redundant, backwards compatible function names
    //to be cleaned up asap
    //no new code should use them
    //they were created as time is short and there is less risk having
    //the functions here rather than adapting all the function call to them in other
    //CMSSW code 
    float dPhi1()const{return dPhiNeg(0);}
    float dPhi1Pos()const{return dPhiPos(0);}
    float dPhi2()const{return dPhiNeg(1);}
    float dPhi2Pos()const{return dPhiPos(1);}
    float dRz1()const{return dRZNeg(0);}
    float dRz1Pos()const{return dRZPos(0);}   
    float dRz2()const{return dRZNeg(1);}
    float dRz2Pos()const{return dRZPos(1);}   
    int subDet1()const{return subDet(0);}
    int subDet2()const{return subDet(1);}
    int hitsMask()const;
    void setNegAttributes(float dRZ2=std::numeric_limits<float>::infinity(),
			  float dPhi2=std::numeric_limits<float>::infinity(),
			  float dRZ1=std::numeric_limits<float>::infinity(),
			  float dPhi1=std::numeric_limits<float>::infinity());
    void setPosAttributes(float dRZ2=std::numeric_limits<float>::infinity(),
			  float dPhi2=std::numeric_limits<float>::infinity(),
			  float dRZ1=std::numeric_limits<float>::infinity(),
			  float dPhi1=std::numeric_limits<float>::infinity());
    
    

  private:
    static float bestVal(float val1,float val2){return std::abs(val1)<std::abs(val2) ? val1 : val2;}
    template<typename T>
    T getVal(size_t hitNr,T PMVars::*val)const{
      return hitNr<hitInfo_.size() ? hitInfo_[hitNr].*val : std::numeric_limits<T>::infinity();
    }
    
  private:

    CtfTrackRef ctfTrack_ ;
    CaloClusterRef caloCluster_ ;
    std::vector<PMVars> hitInfo_;
    int nrLayersAlongTraj_;
    
    bool isEcalDriven_ ;
    bool isTrackerDriven_ ;

  };
}

#endif
