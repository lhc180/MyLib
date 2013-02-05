#ifndef PQMF_H
#define PQMF_H

/**********************************/
// PQMF
/**********************************/
/*********to do list***************/
// 符号化後のデータをitpp::matに対応させる
// clear関数を実装する
/*********************************/

#include<vector>
#include<assert.h>
#include <itpp/base/vec.h>
#include "./mymatrix.h"

const unsigned FILTER_LENGTH = 512;
const unsigned SUBBAND_NUM = 32;


class cPQMF{
private:
  double x[FILTER_LENGTH],y[1024]; // 分析フィルタバッファと合成フィルタのバッファ
  double m[SUBBAND_NUM][64],n[SUBBAND_NUM][64],c[FILTER_LENGTH];
  bool start;			// 合成フィルタがスタートかどうか
  
public:
  cPQMF();			// コンストラクタ
  
				// デフォルトデストラクタ

  // コピーコンストラクタはデフォルト

  void analyze(const std::vector<double> &input, std::vector<double> &output);
  void analyze(const itpp::vec &input, itpp::vec &output);
  void analyze(const std::vector<double> &input, mylib::vec_2D &output_2D);


//   void fastAnalysis(const std::vector<double> &input, std::vector<double> &output);

  // 481サンプルの遅延があるので注意
  void synthesize(const std::vector<double> &input, std::vector<double> &output);
  void synthesize(const itpp::vec &input, itpp::vec &output);
  void synthesize(const mylib::vec_2D &input_2D, std::vector<double> &output);
  
//   void fastSynthesis(const std::vector<double> &input, std::vector<double> &output);
  
};



#endif
