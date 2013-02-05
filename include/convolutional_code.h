#ifndef CONVOLUTIONAL_CODE_H
#define CONVOLUTIONAL_CODE_H

/***************************************/
// This code includes class definition of 
// Convolutional code, punctured CC and RCPC codes. 
//
// 
/****************************************/
// 今のところ符号化率は1/2しか使えない

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>

class cConvolutionalCode{
protected:
  bool init;			// 符号の初期化が完了しているかどうか
  bool termination;		// 終端処理がしてあるかどうか
  // std::vector<int> signal;	// 符号化前の信号
//   std::vector<int> codedSignal;	// 符号化後の信号
  int constraintLength;		// 拘束長
  int inverseRate;		// 符号化率の逆数
  int registerSize;		// constraintLength-1がメモリー数
  int stateNum;			// 状態数
  std::vector<unsigned> generatorPolynomial; // 生成多項式
  unsigned en_shiftRegister;                 // シフトレジスタの中身はビット列で表す 
  
  struct tState{
    int nextState;
    std::vector<int> outputSequence;
  } state;

  virtual void checkInitialization();	// 符号の初期化が完了しているかどうかチェックする

  

  virtual int getNextState(int lastState, int bit);
  virtual void getEncodedSequence(int currentState, int input, std::vector<int> &sequence);
  virtual int get1bit(int currentState);

public:
  cConvolutionalCode():init(false) {}		// デフォルトコンストラクタ
  cConvolutionalCode(const int constraint, const int invRate);
  virtual ~cConvolutionalCode(){ }	// デフォルトデストラクタ

  virtual void initialize();		// 初期化関数

  // 畳み込み符号のステートなどをセットする
  virtual void setCode(const int constraint, const int invRate);

  void encode(const std::vector<int> &input, std::vector<int> &output); // 終端処理せず符号化
  void encode_term(const std::vector<int> &input, std::vector<int> &output); // 終端処理して符号化

  void encode(const itpp::bvec &input, itpp::bvec &output);
  void encode_term(const itpp::bvec &input, itpp::bvec &output);

  // hard decision decoding
  void decode_hd(const std::vector<int> &input, std::vector<int> &output);
  
  // soft decision decoding
  void decode_sd(const std::vector<double> &input, std::vector<int> &output);

  void decode_sd(const itpp::vec &input, itpp::bvec &output);
  
  friend class cRCPCCode;
};

// +++++++++++++ cPuncturedConvolutionalCode ++++++++++++++++++++++
class cPuncturedConvolutionalCode
  : public cConvolutionalCode
{
protected:
  int puncturingPeriod;
  int denominatorCodingRate;	// puncturingPeriod/denominatorCodingRateがパンクチャ後の符号化率
  std::vector<int> puncturingMatrix;
  
  void setPuncturingMatrix();

  itpp::bvec puncturing(const itpp::bvec &input);
  itpp::vec fillingUp(const itpp::vec &input);

public:
  cPuncturedConvolutionalCode() { }
  cPuncturedConvolutionalCode(const int constraint,
                              const int invRate,
                              const int pp,
                              const int denom);
  // デフォルトデストラクタ
  void setCode(const int constraint, const int invRate, const int pp, const int denom);

  void encode(const itpp::bvec &input, itpp::bvec &output);
  void encode_term(const itpp::bvec &input, itpp::bvec &output);

  // Only soft decision.
  void decode_sd(const itpp::vec &input, itpp::bvec &output);

  friend class cRCPCCode;
};

/****************************************/
// RCPC符号のクラス
// cPuncturedConvolutionalCodeのフレンドクラス
//
// cPuncturedConvolutionalCodeのen_shiftRegisterを保存して次の符号化に使う。
/*****************************************/

class cRCPCCode{
private:
  int constraint_;
  int puncturingPeriod_;
  std::vector<int> denomRate_;
  std::vector<int> codeLength_; // 各レイヤーの符号長を格納
  bool init_;

public:
  cRCPCCode(): constraint_(3), puncturingPeriod_(1), denomRate_(0), codeLength_(0), init_(false){ }
  cRCPCCode(int constraint, int puncturingPeriod, const std::vector<int> &denomRate)
  {
    set(constraint, puncturingPeriod, denomRate);
  }
  
  void set(int constraint, int puncturingPeriod, const std::vector<int> &denomRate)
  {
    assert(denomRate.size() != 0);

    constraint_ = constraint;
    puncturingPeriod_ = puncturingPeriod;
    denomRate_ = denomRate;
    codeLength_.resize(denomRate_.size());
    
    init_ = true;
  }

  itpp::bvec encode(const std::vector<itpp::bvec> &input);
  
  std::vector<itpp::bvec> decode(const itpp::vec &output);
};



#endif
