#ifndef CARDDATA_H
#define CARDDATA_H

#include <stdint.h>
#include <iostream>
#include <stdlib.h>

namespace{
  constexpr int BOGUS_INT = std::numeric_limits<int>::max();
  //std::numeric_limits<double>::quiet_NaN()
  constexpr uint64_t BOGUS_UINT64 = std::numeric_limits<uint64_t>::max();
  constexpr uint32_t BOGUS_UINT32 = std::numeric_limits<uint32_t>::max();
  constexpr uint16_t BOGUS_UINT16 = std::numeric_limits<uint16_t>::max();
}

class CardData{
  public:
  uint64_t LastSync;                   // ? clock ticks of last GPS sync pulse
  int SequenceID;                      // unique minibuffer id. should account for event, trigger, minibuffer!
  int StartTimeSec;                    // unix epoch timestamp at StartCount ADC counts
  int StartTimeNSec;                   // above, ns part
  uint64_t StartCount;                 // ADC counts from powerup when ... it gets initialized? 
  static int TriggerNumber;            // number of minibuffers. 40 for Hefty mode. 
  //uint64_t* triggerCounts;           // --
  std::vector<ULong64_t> TriggerCounts; // ADC ticks of the start of each minibuffer in the Full buffer?
  //uint32_t* Rates;                   // --
  std::vector<UInt_t> Rates;            // Avg pulse rate? units? size is number of channels.
  int CardID;                          // card position in VME crate
  static int Channels;                 // num channels on card
  static int BufferSize;               // datapoints per channel in a full buffer.
  static int Eventsize;                // samples per minibuffer per channel divided by 4.
  static int FullBufferSize;           // Channels * BufferSize. corresponds to 80us. 
  //uint16_t* Data;                    // Concatenated data from all channels on card. Size FullBufferSize.
  std::vector<uint16_t> Data;
  // With a fullbuffer of 160k and 4 channels, there's 40k samples per channel, 
  // and with 40 minibuffers per fullbuffer, that's 1k samples per minibuffer per channel.
  // At 2ns per sample that's 2000ns = 2us per minibuffer.
  
  CardData(){ Reset(); }
  void Reset();
  
};

int CardData::TriggerNumber = 40;
int CardData::Channels = 4;
int CardData::BufferSize = 40000;      // UNUSED: calculated in utilityfuncs using trigger window
int CardData::Eventsize = 250;         // UNUSED: as above
int CardData::FullBufferSize = (CardData::Channels * CardData::BufferSize); // UNUSED: as above

void CardData::Reset(){
  LastSync = BOGUS_UINT64;
  SequenceID = BOGUS_INT;
  StartTimeSec = BOGUS_INT;
  StartTimeNSec = BOGUS_INT;
  StartCount = BOGUS_UINT64;
  CardID = BOGUS_INT;
  TriggerCounts.assign(TriggerNumber,BOGUS_UINT64);   // start counts of each minibuffer
  Rates.assign(Channels,BOGUS_UINT32);                // avg pulse rate on each PMT
  Data.assign(FullBufferSize,BOGUS_UINT16);           // 160k datapoints
}

#endif
