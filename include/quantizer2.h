#ifndef QUANTIZER_H
#define QUANTIZER_H

/***********************************/
// スカラー量子化器
//
// 最大振幅値を2^(quantizeSize-1)で割った数が量子化幅
/***********************************/
/*************to do list***********/
// vecでivecに対応する各振幅を保持しておき、double型に対して最も近いvecの値を返す機能を付ける
/**********************************/

#include <iostream>
#include <itpp/base/vec.h>

class cQuantizer{
protected:
  int quantizeSize;
  double maxAmp;
  int divNum;			// 分割数
  itpp::vec quantizeVec;
public:
  cQuantizer(int q = 3, double max = 1.0, bool midrise = false)
  {
    
    set(q, max, midrise);
    
  }
  // デフォルトのコピーコンストラクタ
  // デフォルトデストラクタ

  // 初期化
  // midrise = falseならindex = 0は振幅値0
  void set(int q = 3, double max = 1.0, bool midrise = false);

  double quantize(const double input);
  std::vector<double> quantize(const std::vector<double> &input);
  
  int quantizeReturnIndex(const double input);
  std::vector<int> quantizeReturnIndex(const std::vector<double> &input);
  
  double dequantize(const int index);
  std::vector<double> dequantize(const std::vector<int> &index);

  int getQuantizeSize();
  
};

#endif
