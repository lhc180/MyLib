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

  virtual void initialize();		// 初期化関数

  virtual int getNextState(int lastState, int bit);
  virtual void getEncodedSequence(int currentState, int input, std::vector<int> &sequence);
  virtual int get1bit(int currentState);

public:
  cConvolutionalCode():init(false) {}		// デフォルトコンストラクタ
  cConvolutionalCode(const int constraint, const int invRate);
  virtual ~cConvolutionalCode(){ }	// デフォルトデストラクタ

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
  // cPuncturedConvolutionalCode() { }
  cPuncturedConvolutionalCode(int constraint,
                              int invRate,
                              int pp,
                              int denom);
  // デフォルトデストラクタ
  void setCode(const int constraint, const int invRate, const int pp, const int denom);

  void encode(const itpp::bvec &input, itpp::bvec &output);
  void encode_term(const itpp::bvec &input, itpp::bvec &output);

  // Only soft decision.
  void decode_sd(const itpp::vec &input, itpp::bvec &output);

};




#endif
