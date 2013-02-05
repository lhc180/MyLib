#ifndef MLC_MSD_H
#define MLC_MSD_H

#include <vector>
#include "../include/myldpc.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
/*************************************************/
// �����Multilevel Coding, Multistage Decoding�̃N���X�ł���
// ���������@�͍��̂Ƃ���LDPC�݂̂�z�肵�Ă���
// ���ꂼ��ʁX�̃N���X�ŁA��������
//
/***********************************************/
/***************** to do list ******************/
// error concealment�ɑΉ�������
// 
/**********************************************/


/***************** class cMLC_MSD *********************/
class cMLC_MSD{
protected:
  std::vector<cLDPC> vecLDPC;
  itpp::Modulator_2D Modulator;
  int nLevel;			// bitsPerSymbol�ƈꏏ
  unsigned nCodeLength;
  double avrCodeRate;		// �V���{�����[�g
  bool bSetDone;

public:
  cMLC_MSD() : nLevel(0), nCodeLength(0), avrCodeRate(0), bSetDone(false)
  { }
  cMLC_MSD(std::vector<cLDPC> &vLDPC, itpp::Modulator_2D &Mod);

  void set(std::vector<cLDPC> &vLDPC, itpp::Modulator_2D &Mod);
  
  void encode(const std::vector<itpp::bvec> &vBits, itpp::cvec &symbols);

  itpp::cvec encode(const std::vector<itpp::bvec> &vBits)
  {
    itpp::cvec symbols;
    encode(vBits, symbols);
    return symbols;
  }
  
  // �e���x���̃C�e���[�V�����񐔂�Ԃ�
  itpp::ivec decode(const itpp::cvec &symbols,
		    std::vector< itpp::bvec > &vDecodedBits,
		    double N0,
		    unsigned loopMax);

  std::vector<itpp::bvec> decode(const itpp::cvec &symbols,
                                 double N0,
                                 unsigned loopMax)
  {
    std::vector< itpp::bvec > vDecodedBits(Modulator.bits_per_symbol());
    decode(symbols, vDecodedBits, N0, loopMax);
    return vDecodedBits;
  }
  
  
  double symbolRate(){
    return avrCodeRate;
  }
};

/*************** end of cMLC_MSD ********************/

// MLC�p��LDPC���Z�b�g����
// ���ϕ���������Ԃ�
double setLDPC_MLC(std::vector< cLDPC > &vecLDPC,
		   unsigned nCodeLength,
		   itpp::ivec &vecRowWeight,
		   itpp::ivec &vecColWeight);

// �Z�b�g�p�[�e�B�V���j���O�}�b�s���O�ɂ���
void setModulatorNatural(itpp::Modulator_2D &Mod);

// Proposals for MLC and MSD below.
// These adjust lengths of info bits such that higher level info bits move to
// lower level.
std::vector< itpp::bvec > adjustInfoLength_forMLC(const std::vector< itpp::bvec > &infoBits,
                                                  const std::vector< int > &infoLengths_LDPC);

std::vector< itpp::bvec > adjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo,
                                                   const std::vector< int > &infoLengths);

// end of proposals.
                                                   
// Divide input into some levels.
template< typename kind >
mylib::Vector_2D< kind > divideIntoLevels(const std::vector< kind > &input, const std::vector< int > &numSamples)
{
  int sum = 0;
  for(int i = 0; i < numSamples.size(); i++)
    {
      sum += numSamples[i];
    }
  assert(input.size() == sum);

  mylib::Vector_2D< kind > output(numSamples.size());
  int inputIndex = 0;
  for(int level = 0; level < numSamples.size(); level++)
    {
      output(level).resize(numSamples[level]);
      for(int i = 0; i < numSamples[level]; i++)
        {
          output(level, i) = input[inputIndex];
          inputIndex++;
        }
    }

  return output;
}

#endif
