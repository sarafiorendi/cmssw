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
  
// 
//   /// Detector element
//   DetId getDetId() const { return theDetId; }
//   void setDetId(DetId aDetId) { theDetId = aDetId; }
//   unsigned int getStackMember() const { return theStackMember; }
//   void setStackMember(unsigned int aStackMember) { theStackMember = aStackMember; }
// 
//   /// Rows and columns to get rid of Digi collection
//   std::vector<int> findRows() const;
//   std::vector<int> findCols() const;
//   void setCoordinates(std::vector<int> a, std::vector<int> b) {
//     theRows = a;
//     theCols = b;
//   }
//   std::vector<int> getRows() const { return theRows; }
//   std::vector<int> getCols() const { return theCols; }

//   /// Single hit coordinates
//   /// Information
//   std::string print(unsigned int i = 0) const;

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

/*! \brief   Implementation of methods
 *  \details Here, in the header file, the methods which do not depend
 *           on the specific type <T> that can fit the template.
 *           Other methods, with type-specific features, are implemented
 *           in the source file.
 */

/// Default Constructor
/// NOTE: to be used with setSomething(...) methods
// L1StubTriplet::L1StubTriplet() {
//   /// Set default data members
// //   theHits.clear();
//   innerStubRapprox_ = 0.;
// }

// /// Another Constructor
// template <typename T>
// TTCluster<T>::TTCluster(std::vector<T> aHits, DetId aDetId, unsigned int aStackMember, bool storeLocal) {
//   /// Set data members
//   this->setHits(aHits);
//   this->setDetId(aDetId);
//   this->setStackMember(aStackMember);
// 
//   theRows.clear();
//   theCols.clear();
//   if (storeLocal) {
//     this->setCoordinates(this->findRows(), this->findCols());
//   }
// }

/// Destructor
// L1StubTriplet::~L1StubTriplet() {}

/// Information
// template <typename T>
// std::string TTCluster<T>::print(unsigned int i) const {
//   std::string padding("");
//   for (unsigned int j = 0; j != i; ++j) {
//     padding += "\t";
//   }
// 
//   std::stringstream output;
//   output << padding << "TTCluster:\n";
//   padding += '\t';
//   output << padding << "DetId: " << theDetId.rawId() << '\n';
//   output << padding << "member: " << theStackMember << ", cluster size: " << theHits.size() << '\n';
//   return output.str();
// }
// 
// template <typename T>
// std::ostream& operator<<(std::ostream& os, const TTCluster<T>& aTTCluster) {
//   return (os << aTTCluster.print());
// }

#endif
