#ifndef LDPC_SPARSE_H
#define LDPC_SPARSE_H
/*************************************/
// ボツになりそう
/*************************************/


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
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <vector>
#include <cstdlib>


typedef struct tIndex1{
    std::vector<int> indexVec;
} index1;

class ldpc{
private:
  int signOfNumber(double input); // inputの符号を返す
  double fFunction(double x);	     // Gallagerのf関数
  
  void rowsProcessing(itpp::sparse_mat &alpha, itpp::sparse_mat &beta,
		      const itpp::vec &llrVec,
		      const std::vector<index1> &rowsIndex); // 復号過程における行処理
  
  void colsProcessing(itpp::sparse_mat &alpha, itpp::sparse_mat &beta,
		      const std::vector<index1> &colsIndex); // 列処理

  template< typename type >
  itpp::vec calcLLR(const itpp::Modulator< type > &mod,
		    const itpp::Vec< type > &receivedVec,
		    const itpp::bvec &estimatedVec,
		    const double N0);

  void estimateCode(const itpp::sparse_mat &alpha, 
		    itpp::bvec &decoded,
		    const itpp::vec &llrVec,
		    const std::vector<index1> &colsIndex); // 一時推定語を求める

  bool checkParity(itpp::bvec &decoded); // パリティ検査

  
protected:
  unsigned HRowSize;		// 検査行列の行数
  unsigned HColSize;		// 検査行列の列数(符号長)
  unsigned HRowWeight;		// 検査行列の行重み
  unsigned HColWeight;		// 検査行列の列重み
  itpp::GF2mat_sparse Hmat;		// Gallager法による検査行列
  itpp::GF2mat_sparse Gmat;		// 生成行列
  itpp::ivec perm;		// Hmatの列入れ替え情報を保持
  
public:
  ldpc();			// デフォルトコンストラクタ

  ldpc(unsigned n, unsigned rowWeight = 6, unsigned colWeight = 3):
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
  
  void setEncoder();		// パリティ検査行列の生成から
				// 生成行列の生成までを行う

  void makeParityCheckMatrix();	// パリティ検査行列を作る
  
  void printParityCheckMatrix(); // パリティ検査行列を出力

  void makeGeneratorMatrix();	// 生成行列を作る

  void checkParityCondition(); 	// GHt=0かどうかチェックする

  unsigned getCodeLength();	// 符号長(Hの列数)を返す

  unsigned getCheckSymbolLength(); // 検査記号数(Hの行数)を返す

  unsigned getInputLength();	// 符号化される情報ベクトルの長さ
  
  double getRate();		// 符号化率

  void encoding(const itpp::bvec &input, itpp::bvec &coded); // 符号化

  // iteration回数を返す
  template <typename type>
  unsigned decoding(const itpp::Modulator< type > &mod,
		    const itpp::Vec< type > &symbol,
		    itpp::bvec &decodedBits, 
		    const double N0,
		    const unsigned loopMax = 100); // 対数領域sum-product復号法
  // loopMaxは最大反復回数
};




#endif
