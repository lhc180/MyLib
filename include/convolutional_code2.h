#ifndef CONVOLUTIONAL_CODE_H
#define CONVOLUTIONAL_CODE_H

/***************************************/
// This code includes class definition of 
// Convolutional code, punctured CC and RCPC codes. 
//
// 
/****************************************/
// ���̂Ƃ��땄��������1/2�����g���Ȃ�

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>

class cConvolutionalCode{
protected:
  bool init;			// �����̏��������������Ă��邩�ǂ���
  bool termination;		// �I�[���������Ă��邩�ǂ���
  // std::vector<int> signal;	// �������O�̐M��
//   std::vector<int> codedSignal;	// ��������̐M��
  int constraintLength;		// �S����
  int inverseRate;		// ���������̋t��
  int registerSize;		// constraintLength-1���������[��
  int stateNum;			// ��Ԑ�
  std::vector<unsigned> generatorPolynomial; // ����������
  unsigned en_shiftRegister;                 // �V�t�g���W�X�^�̒��g�̓r�b�g��ŕ\�� 
  
  struct tState{
    int nextState;
    std::vector<int> outputSequence;
  } state;

  virtual void checkInitialization();	// �����̏��������������Ă��邩�ǂ����`�F�b�N����

  virtual void initialize();		// �������֐�

  virtual int getNextState(int lastState, int bit);
  virtual void getEncodedSequence(int currentState, int input, std::vector<int> &sequence);
  virtual int get1bit(int currentState);

public:
  cConvolutionalCode():init(false) {}		// �f�t�H���g�R���X�g���N�^
  cConvolutionalCode(const int constraint, const int invRate);
  virtual ~cConvolutionalCode(){ }	// �f�t�H���g�f�X�g���N�^

  // ��ݍ��ݕ����̃X�e�[�g�Ȃǂ��Z�b�g����
  virtual void setCode(const int constraint, const int invRate);

  void encode(const std::vector<int> &input, std::vector<int> &output); // �I�[��������������
  void encode_term(const std::vector<int> &input, std::vector<int> &output); // �I�[�������ĕ�����

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
  int denominatorCodingRate;	// puncturingPeriod/denominatorCodingRate���p���N�`����̕�������
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
  // �f�t�H���g�f�X�g���N�^
  void setCode(const int constraint, const int invRate, const int pp, const int denom);

  void encode(const itpp::bvec &input, itpp::bvec &output);
  void encode_term(const itpp::bvec &input, itpp::bvec &output);

  // Only soft decision.
  void decode_sd(const itpp::vec &input, itpp::bvec &output);

};




#endif
