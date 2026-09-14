#ifndef EventFilter_Phase2TrackerRawToDigi_TrackerBlock_H
#define EventFilter_Phase2TrackerRawToDigi_TrackerBlock_H

#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2DAQFormatSpecification.h"
#include <vector>
#include <iostream>
#include <bitset>

using namespace Phase2DAQFormatSpecification;

// Base class to be then specialised for tracker header and trailer
class TrackerBlock {
public:
  explicit TrackerBlock(size_t nLines) : values_(nLines, 0) {}
  virtual ~TrackerBlock() = default;

  void setValue(std::vector<uint32_t>& newValues) {
    values_ = newValues;
    is2S_ = ((values_[0] >> (N_BITS_PER_WORD - MODULE_TYPE_BITS)) & ((1u << MODULE_TYPE_BITS) - 1)) == MODULE_TYPE_2S;
    setSpecificValue();
  }

  bool is2S() const { return is2S_; }

  void printValues() const {
    for (size_t i = 0; i < values_.size(); ++i)
      printValue(i);
  }

  void printValue(size_t i) const {
    std::cout << blockName() << "[" << i << "]: " << values_[i] << "   " << std::bitset<32>(values_[i]) << std::endl;
  }

protected:
  std::vector<uint32_t> values_;
  bool is2S_{false};

  virtual std::string blockName() const { return "TrackerBlock"; }
  virtual void setSpecificValue() {} // overridden by sub-classes

};

/**
 * @brief: Tracker Header Object, decodes Tracker Header described in the link below:
 * https://docs.google.com/spreadsheets/d/1RHZFqeHCoJhRaAfaKEO1Gx6U6c1Y3tRGhL_aSbZQROY/edit?gid=256168213#gid=256168213
 */
class TrackerHeader : public TrackerBlock {
public:
  TrackerHeader() : TrackerBlock(HEADER_N_LINES) {}

  explicit TrackerHeader(const std::vector<uint32_t>& words) : TrackerBlock(HEADER_N_LINES) {
    std::vector<uint32_t> copy(words);
    setValue(copy);
  }

  uint8_t hasExtendedData() const { return hasED_; }
  ////  FIXME remove numbers 
  uint32_t getBoardType() const { return (values_[0] >> 24) & 0xFF; }      // bits 31-24 (8 bits)
  uint32_t getVersionMajor() const { return (values_[0] >> 21) & 0x7; }    // bits 23-21 (3 bits)
  uint32_t getVersionMinor() const { return (values_[0] >> 16) & 0x1F; }   // bits 20-16 (5 bits)
  uint32_t getMode() const { return (values_[0] >> 13) & 0x7; }            // bits 15-13 (3 bits)
  uint32_t getED() const { return (values_[0] >> 12) & 0x1; }              // bit 12 (1 bit)
  uint32_t getBoardID() const { return (values_[0] >> 4) & 0xFF; }         // bits 11-4 (8 bits)
  uint32_t getDAQpathCoreID() const { return values_[0] & 0xF; }           // bits 3-0 (4 bits)

  void printFields() const {
    printf(
        "Board Type: %02X, Ver Major: %03d, Ver Minor: %04d, Mode: %03d, ED: %01d, Board ID: %02X, DAQpath Core "
        "ID: %01X\n",
        getBoardType(),
        getVersionMajor(),
        getVersionMinor(),
        getMode(),
        getED(),
        getBoardID(),
        getDAQpathCoreID());
  }

protected:
  std::string blockName() const override { return "TrackerHeader"; }

  void setSpecificValue() override {
    hasED_ = static_cast<uint8_t>(getED());
  }

private: 
  uint8_t hasED_{0}; // indicates if extended data are present

};



class TrackerTrailer : public TrackerBlock {
public:
  TrackerTrailer() : TrackerBlock(TRAILER_N_LINES) {}

  uint8_t afterExtendedData() const { return endED_; }

protected:
  std::string blockName() const override { return "TrackerTrailer"; }

  void setSpecificValue() override {
    endED_ = (values_[0]) & ((1u << TRAILER_ENDED_BITS) - 1);
  }

private:
  // indicates if trailer is at end of a normal or extended data
  uint8_t endED_{0};  
  
};

#endif