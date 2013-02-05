#ifndef ADPCM_LMS_H
#define ADPCM_LMS_H

/****************************************/
// このソースコードは以下の文献で提案されたADPCMを実装してある
// ADPCM WITH ADAPTIVE PRE- AND POST-FILTERING FOR DELAY-FREE AUDIO CODING
// 変数や配列の名前は論文を参照
/***************************************/


#include <itpp/base/vec.h>
#include "./adpcm_delayfree.h"


class adpcmLMS 
  : public adpcmDelayFree
{
protected:
  virtual double predictor(const double xHat, const double eHat);
  
public:
  // adpcmDelayFree();	      // ------------ デフォルトコンストラクタ
  adpcmLMS(int q = 4, double max = 1.5, int d = 30, double lAT = 0.9, double lRT = 0.1):
    adpcmDelayFree(q,max,d,lAT,lRT)
  {}
  // ~adpcmDelayFree();------------ デフォルトデストラクタ
  
  virtual void encoding(const itpp::vec &inputData, itpp::vec &encodedData);
  virtual void encoding(const itpp::vec &inputData, itpp::ivec &outputData);
  virtual void decoding(const itpp::vec &encodedData, itpp::vec &decodedData);
  virtual void decoding(const itpp::ivec &encodedData, itpp::vec &decodedData);

};


#endif
