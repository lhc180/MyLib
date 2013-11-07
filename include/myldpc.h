#ifndef MYLDPC_H
#define MYLDPC_H

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
#include <mymatrix.h>
#include <vector>
#include <cstdlib>

namespace mylib{
  // マルチスレッドの戻り値用に使う
  struct DecoderParas {
    itpp::bvec decodedBits;
    int loop;
  };

  class Ldpc{
  protected:
    virtual void SetupMatrix(); // パリティ検査行列の生成から
				// 生成行列の生成までを行う

    virtual void MakeParityCheckMatrix(itpp::GF2mat &tHmat);	// パリティ検査行列を作る

    virtual void MakeGeneratorMatrix(itpp::GF2mat &tHmat, itpp::GF2mat &tGmat);	// 生成行列を作る

    virtual void CheckParityCondition(itpp::GF2mat &tHmat, itpp::GF2mat &tGmat); 	// GHt=0かどうかチェックする


    virtual double Ffunction(double x);	     // Gallagerのf関数
  
    virtual void RowsProcessing(itpp::mat* alpha, const itpp::mat &beta,
                                const itpp::vec &llrVec); // 復号過程における行処理
  
    virtual void ColsProcessing(const itpp::mat &alpha, itpp::mat* beta); // 列処理

    virtual itpp::vec CalcLLR(const itpp::Modulator_2D &mod,
                              const itpp::cvec &receivedVec,
                              const itpp::bvec &estimatedVec,
                              double N0);   

    virtual itpp::vec CalcLLRWithPads(const itpp::Modulator_2D& mod,
                                      const itpp::cvec& receivedVec,
                                      const itpp::bvec& estimatedVec,
                                      int numPads,
                                      double n0);
    

    
    
    virtual void EstimateCode(const itpp::mat &alpha, 
                              itpp::bvec* decoded,
                              const itpp::vec &llrVec); // 一時推定語を求める

    virtual bool CheckParity(const itpp::bvec &decoded); // パリティ検査

  
  protected:
    int hRowSize_;		// 検査行列の行数
    int hColSize_;		// 検査行列の列数(符号長)
    int hRowWeight_;		// 検査行列の行重み
    int hColWeight_;		// 検査行列の列重み
    int infoLength_;
    mylib::ivec_2D hMat_;	// Gallager法による検査行列
				// 1である要素の番号だけを格納する
    mylib::ivec_2D hMatTrans_; // Hmatの転置行列
    mylib::ivec_2D gMat_;	// 生成行列
				// 同様に1である要素の番号だけを格納する
    mylib::ivec_2D gMatTrans_; // Gmatのtranspose
    itpp::ivec perm_;		// Hmatの列入れ替え情報を保持
    bool setDone_;

  public:
    Ldpc(): setDone_(false) { }			// デフォルトコンストラクタ

    Ldpc(unsigned codeLength, unsigned rowWeight = 6, unsigned colWeight = 3)
    {
      assert(codeLength > 0 && rowWeight > 0 && colWeight > 0);
      Set(codeLength, rowWeight, colWeight);
    } // コンストラクタ

    virtual ~Ldpc() { } // --- デフォルトデストラクタ

    //   void initialize(unsigned k, unsigned n, unsigned rowWeight = 6, unsigned colWeight = 3);
    void Set(unsigned n, unsigned rowWeight = 6, unsigned colWeight =3);
  
    //  virtual void printParityCheckMatrix(); // パリティ検査行列を出力

  
    unsigned CodeLength()	// 符号長(Hの列数)を返す
    {
      return hColSize_;
    }

    unsigned SymbolLength() // 検査記号数(Hの行数)を返す
    {
      return hRowSize_;
    }

    unsigned InfoLength()	// 符号化される情報ベクトルの長さ
    {
      return infoLength_;
    }

    double CodeRate()		// 符号化率
    {
      return static_cast<double>(infoLength_)/static_cast<double>(hColSize_);
    }

    
    void Encode(const itpp::bvec &input, itpp::bvec &coded); // 符号化
    itpp::bvec Encode(const itpp::bvec& input);
    
    // iteration回数を返す

    int Decode(const itpp::Modulator_2D &mod,
                       const itpp::cvec &symbol,
                       itpp::bvec &decodedBits, 
                       double n0,
                       int loopMax = 100); // 対数領域sum-product復号法
    // loopMaxは最大反復回数
    itpp::bvec Decode(const itpp::Modulator_2D& mod,
                              const itpp::cvec& symbol,
                              double n0, int loopMax = 100)
    {
      itpp::bvec decoded;
      Decode(mod, symbol, decoded, n0, loopMax);
      return decoded;
    }

    // 受信器側でpadding bitsの数が分かっているとき
    int DecodeWithPadding0(const itpp::Modulator_2D& mod,
                                   const itpp::cvec& symbol,
                                   itpp::bvec& decodedBits,
                                   double N0,
                                   int numPad = 0,
                                   int loopMax = 100);

    
    
    
  };

  class LdpcForMlcMsd: public Ldpc
  {
  protected:
        // This is for MLC and MSD.
    virtual itpp::vec CalcLLR(const std::vector< itpp::Modulator_2D > &vecMod,
                              const itpp::cvec &receivedVec,
                              double N0);

    // This is for multithread MSD.
    virtual itpp::vec CalcLLRAtLevel(const itpp::Modulator_2D &mod,
                              const itpp::cvec &receivedVec,
                              double N0,
                              int level);

    
  public:
    LdpcForMlcMsd() { }
    LdpcForMlcMsd(unsigned codeLength, unsigned rowWeight = 6, unsigned colWeight = 3):
      Ldpc(codeLength, rowWeight, colWeight)
    { }
    
    virtual ~LdpcForMlcMsd() { }

        // これはMLC、MSD用
    int Decode(const std::vector< itpp::Modulator_2D > &vecMod,
                       const itpp::cvec &symbol,
                       itpp::bvec &decodedCodes, // 符号語であり、情報ビットではない
                       const double N0,
                       const int loopMax = 100);
    
    // 並列処理するMSD用
    DecoderParas DecodeAtLevel(const itpp::Modulator_2D& vecMod,
                                     const itpp::cvec &symbol,
                                     // itpp::bvec &decodedCodes,
                                     double N0,
                                       int level, // どのレベルの復号か．
                                     // BlockPartitioning8PSKなら0か
                                     // 1が並列処理できる
                                       // int *loop,
                                     int loopMax = 100);


  };
  
}



#endif
