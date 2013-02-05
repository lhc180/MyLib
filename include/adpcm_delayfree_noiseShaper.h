#ifndef ADPCM_DELAYFREE_NOISESHAPER_H
#define ADPCM_DELAYFREE_NOISESHAPER_H

/****************************************/
// このソースコードは以下の文献で提案されたADPCMを実装してある
// DELAY-FREE AUDIO CODING BASED ON ADPCM AND ERROR FEEDBACK
// noise shaperを継承を使って実装する
// 変数や配列の名前は論文を参照
//
// GAL algorithmによる更新式は以下の論文を参照した
// Simulation and Performance Analysis of Adaptive Filtering Algorithms in Noise Cancellation
/***************************************/
/****************to do list*************/
// noiseShaperを作る
/***************************************/

#include <itpp/base/vec.h>
#include "./adpcm_delayfree.h"


class adpcmWithNoiseShaper: public adpcmDelayFree  
{
private:
    
protected:
  itpp::vec hs;			// noise shaper内のフィルタ係数
  virtual double noiseShaper(const double u, const double x, const double xHat);

public:
  // adpcmWithNoiseShaper();	      // ------------ デフォルトコンストラクタ
  
  adpcmWithNoiseShaper(int q = 4, double max = 1.5, int d = 50, double lAT = 0.9, double lRT = 0.1):
    adpcmDelayFree(q, max, d, lAT, lRT), hs(9)
  {
    hs[0] = 0.0;
    hs[1] = 0.5;
    hs[2] = 0.04;
    hs[3] = 0.0245;
    hs[4] = 6.5*pow(10,-4);
    hs[5] = 0.274;
    hs[6] = -0.1;
    hs[7] = -2.5*pow(10,-3);
    hs[8] = 0.218;
  }
  // ~adpcmDelayFree();------------ デフォルトデストラクタ

  virtual void encoding(const itpp::vec &input, itpp::vec &encoded);
  virtual void encoding(const itpp::vec &input, itpp::ivec &encoded);
};
  
  


#endif
