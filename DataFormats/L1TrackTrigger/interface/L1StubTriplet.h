/*! \class   L1StubTriplet
 *
 */

#ifndef L1STUBTRIPLET_H
#define L1STUBTRIPLET_H

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

class L1StubTriplet {
public:
  /// Constructors
  L1StubTriplet();

  /// Destructor
  ~L1StubTriplet(); 

  float getStubRapprox(int i) const { return rapprox_.at(i); }
  void setStubRapprox(int i, float rapprox) {rapprox_[i] = rapprox; }

  float getStubZapprox(int i) const { return zapprox_.at(i); }
  void setStubZapprox(int i, float zapprox) {zapprox_[i] = zapprox; }

  float getStubBend(int i) const { return bend_.at(i); }
  void setStubBend(int i, float bend) {bend_[i] = bend; }

  int getStubRZbin(int i) const { return rzbin_.at(i); }
  void setStubRZbin(int i, int rzbin) {rzbin_[i] = rzbin; }

  int getStubIndex(int i) const { return index_.at(i); }
  void setStubIndex(int i, int idx) {index_[i] = idx; }

  int getStubLayerdisk(int i) const { return layerdisk_.at(i); }
  void setStubLayerdisk(int i, int layerdisk) {layerdisk_[i] = layerdisk; }

  unsigned int getSector() const { return sector_; }
  void setSector(unsigned int i) { sector_ = i; }

  int getRegion() const { return region_; }
  void setRegion(int i) { region_ = i; }

  int getTPDUnit() const { return tpdunit_; }
  void setTPDUnit(int i) { tpdunit_ = i; }
  
//   int rzbin = (outervmstub.vmbits().value() & (settings_->NLONGVMBINS() - 1));

private:

  /// Data members
  std::vector<float> rapprox_;
  std::vector<float> zapprox_;
  std::vector<float> bend_;
  std::vector<int> rzbin_;
  std::vector<int> index_;
  std::vector<int> layerdisk_;

  unsigned int sector_;
  int region_;
  int tpdunit_;
  
};  /// Close class
#endif
