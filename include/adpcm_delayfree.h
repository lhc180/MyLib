#ifndef ADPCM_DELAYFREE_H
#define ADPCM_DELAYFREE_H

/****************************************/
// このソースコードは以下の文献で提案されたADPCMを実装してある
// DELAY-FREE AUDIO CODING BASED ON ADPCM AND ERROR FEEDBACK
// 変数や配列の名前は論文を参照
/***************************************/

#include "./quantizer.h"
#include <itpp/base/vec.h>

const double _QUANTIZE_MAX_AMP = 1.5;

class adpcmDelayFree{

protected:
  // int quantizeSize;		// 量子化サイズ  
  // itpp::vec quantizeVec;		// 
  cQuantizer Quantizer;
  int delay;			// 最大遅延
  double lambdaAT,lambdaRT;		// ローパスフィルタの係数
  bool start;				// 各プライベート関数を初期化するかどうか
  itpp::vec kappa;
  // virtual double quantize(const double in);
  // virtual int quantizeReturnIndex(const double in); // quantizeVecのインデックスを返す
  virtual double envelopeEstimation(const double v, const double eHat); 
  virtual double predictor(const double xHat);

public:
  // adpcmDelayFree();	      // ------------ デフォルトコンストラクタ
  adpcmDelayFree(int q = 3, double max = 1.5, int d = 50, double lAT = 0.9, double lRT = 0.1):
    // quantizeSize(q), quantizeVec(static_cast<int>(pow(2,q))), 
    Quantizer(q,max),
    delay(d), lambdaAT(lAT), lambdaRT(lRT), start(true), kappa(delay+1)
  {
    // 以下のコメントは論文の量子化器
    // if(q == 3){
//       quantizeVec[0] = 0.15;
//       quantizeVec[1] = 0.55;
//       quantizeVec[2] = 1.05;
//       quantizeVec[3] = 1.80;
//       quantizeVec[7] = -1*quantizeVec[0];
//       quantizeVec[6] = -1*quantizeVec[1];
//       quantizeVec[5] = -1*quantizeVec[2];
//       quantizeVec[4] = -1*quantizeVec[3];
//     }
//     else if(q == 4){
//       quantizeVec[0] = 0.1;
//       quantizeVec[1] = 0.3;
//       quantizeVec[2] = 0.52;
//       quantizeVec[3] = 0.788;
//       quantizeVec[4] = 1.0;
//       quantizeVec[5] = 1.313;
//       quantizeVec[6] = 1.738;
//       quantizeVec[7] = 2.70;
//       int j = 0;
//       for(int i = quantizeVec.size()-1; i >= 8; i--,j++){
// 	quantizeVec[i] = -1*quantizeVec[j];
//       }
   //  }
//     else{
//       std::cerr << "Error: Quantize size is not supported.\n";
//       exit(1);
//     }
  
  }

  virtual ~adpcmDelayFree(){};	// ------------ デフォルトデストラクタ
  
  virtual void encoding(const itpp::vec &inputData, itpp::vec &encodedData);
  // 量子化されたデータを出力
  virtual void encoding(const itpp::vec &inputData, itpp::ivec &encodedData);
  
  virtual void decoding(const itpp::vec &encodedData, itpp::vec &decodedData);
  // 量子化されたデータから復号
  virtual void decoding(const itpp::ivec &encodedData, itpp::vec &decodingFromQuantizedData);

};


#endif
