#ifndef MY_LDPC_H
#define MY_LDPC_H

/************************************************/
// Gallager法によって検査行列を作る
// 行列演算はitppのものを使う
// 通信路はAWGNを想定
// 符号化変調に対応
// ただし、復号の際に変調方式等の情報を渡す
// 
// H2が実際に使う検査行列
// f関数は高速化のため近似する
/************************************************/

/***************to do list *********************/
// 行処理と列処理はもっと高速化できそう。
// Hmatの行数などのメンバ変数は無くても良いかも。
// 符号化変調用のdecoding関数も多重定義で作る。
// 和田山先生のやり方でやる。
// メンバ関数にlambdaの計算を入れる。
/************************************************/

#include <iostream>
#include <itpp/base/gf2mat.h>
#include <itpp/base/vec.h>
#include <itpp/base/specmat.h>
#include <itpp/comm/modulator.h>
#include <vector>
#include <cstdlib>

typedef struct tIndex1{
    std::vector<int> indexVec;
} index1;

class cLDPC{
protected:
  virtual int signOfNumber(double input); // inputの符号を返す
  virtual double fFunction(double x);	     // Gallagerのf関数
  
  virtual void rowsProcessing(itpp::mat &alpha, itpp::mat &beta,
			      const itpp::vec &llrVec,
			      const std::vector<index1> &rowsIndex); // 復号過程における行処理
  
  virtual void colsProcessing(itpp::mat &alpha, itpp::mat &beta,
			      const std::vector<index1> &colsIndex); // 列処理

  virtual itpp::vec calcLLR(const itpp::Modulator_2D &mod,
			    const itpp::cvec &receivedVec,
			    const itpp::bvec &estimatedVec,
			    const double N0);

  // This is for MLC and MSD.
  virtual itpp::vec calcLLR(const std::vector< itpp::Modulator_2D > &vecMod,
			    const itpp::cvec &receivedVec,
			    const double N0);


  virtual void estimateCode(const itpp::mat &alpha, 
			    itpp::bvec &decoded,
			    const itpp::vec &llrVec,
			    const std::vector<index1> &colsIndex); // 一時推定語を求める

  virtual bool checkParity(itpp::bvec &decoded); // パリティ検査

  
protected:
  unsigned HRowSize;		// 検査行列の行数
  unsigned HColSize;		// 検査行列の列数(符号長)
  unsigned HRowWeight;		// 検査行列の行重み
  unsigned HColWeight;		// 検査行列の列重み
  itpp::GF2mat Hmat;		// Gallager法による検査行列
  itpp::GF2mat Gmat;		// 生成行列
  itpp::ivec perm;		// Hmatの列入れ替え情報を保持
  bool bSetDone;

public:
  cLDPC(): bSetDone(false) { }			// デフォルトコンストラクタ

  cLDPC(unsigned n, unsigned rowWeight = 6, unsigned colWeight = 3):
    HRowSize(n/rowWeight*colWeight), HColSize(n), HRowWeight(rowWeight), HColWeight(colWeight) 
  {
    // 列数が行重みの倍数でなかったら失敗
    if(HColSize%HRowWeight != 0 || HColSize < HRowSize){
      std::cout << "You must choose correct code length.\n";
      exit(1);
    }

    setEncoder();

  } // コンストラクタ

  // ~ldpc() --- デフォルトデストラクタ

  //   void initialize(unsigned k, unsigned n, unsigned rowWeight = 6, unsigned colWeight = 3);
  virtual void setEncoder(unsigned n, unsigned rowWeight = 6, unsigned colWeight =3);

  virtual void setEncoder();		// パリティ検査行列の生成から
  // 生成行列の生成までを行う

  virtual void makeParityCheckMatrix();	// パリティ検査行列を作る
  
  virtual void printParityCheckMatrix(); // パリティ検査行列を出力

  virtual void makeGeneratorMatrix();	// 生成行列を作る

  virtual void checkParityCondition(); 	// GHt=0かどうかチェックする

  virtual unsigned getCodeLength();	// 符号長(Hの列数)を返す

  virtual unsigned getCheckSymbolLength(); // 検査記号数(Hの行数)を返す

  virtual unsigned getInputLength();	// 符号化される情報ベクトルの長さ
  
  virtual double getRate();		// 符号化率

  virtual void encode(const itpp::bvec &input, itpp::bvec &coded); // 符号化

  // iteration回数を返す
  virtual unsigned decode(const itpp::Modulator_2D &mod,
			  const itpp::cvec &symbol,
			  itpp::bvec &decodedBits, 
			  const double N0,
			  const unsigned loopMax = 100); // 対数領域sum-product復号法
  // loopMaxは最大反復回数

  // これはMLC、MSD用
  virtual unsigned decode(const std::vector< itpp::Modulator_2D > &vecMod,
			  const itpp::cvec &symbol,
			  itpp::bvec &decodedCodes, // 符号語であり、情報ビットではない
			  const double N0,
			  const unsigned loopMax = 100);

};




#endif
