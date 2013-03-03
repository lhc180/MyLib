#ifndef MLC_MSD_H
#define MLC_MSD_H

#include <vector>
#include "../include/myldpc.h"
#include "../include/myutl.h"
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
namespace mylib{
  class MlcMsd{
  protected:
    std::vector<Ldpc> vecLDPC;
    itpp::Modulator_2D Modulator;
    int nLevel;			// bitsPerSymbol�ƈꏏ
    unsigned nCodeLength;
    double avrCodeRate;		// �V���{�����[�g
    bool bSetDone;

  public:
    MlcMsd() : nLevel(0), nCodeLength(0), avrCodeRate(0), bSetDone(false)
    { }
    MlcMsd(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod);

    void Set(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod);
  
    void Encode(const std::vector<itpp::bvec> &vBits, itpp::cvec &symbols);

    itpp::cvec Encode(const std::vector<itpp::bvec> &vBits)
    {
      itpp::cvec symbols;
      Encode(vBits, symbols);
      return symbols;
    }
  
    // �e���x���̃C�e���[�V�����񐔂�Ԃ�
    itpp::ivec Decode(const itpp::cvec &symbols,
                      std::vector< itpp::bvec > &vDecodedBits,
                      double N0,
                      unsigned loopMax);

    std::vector<itpp::bvec> Decode(const itpp::cvec &symbols,
                                   double N0,
                                   unsigned loopMax)
    {
      std::vector< itpp::bvec > vDecodedBits(Modulator.bits_per_symbol());
      Decode(symbols, vDecodedBits, N0, loopMax);
      return vDecodedBits;
    }
  
  
    double SymbolRate(){
      return avrCodeRate;
    }
  };

  /*************** end of cMLC_MSD ********************/

  // MLC�p��LDPC���Z�b�g����
  // ���ϕ���������Ԃ�
  double SetLDPC_MLC(std::vector< Ldpc > &vecLDPC,
                     unsigned nCodeLength,
                     itpp::ivec &vecRowWeight,
                     itpp::ivec &vecColWeight);

  // �Z�b�g�p�[�e�B�V���j���O�}�b�s���O�ɂ���
  void SetModulatorNatural(itpp::Modulator_2D &Mod);

  // Proposals for MLC and MSD below.
  // These adjust lengths of info bits such that higher level info bits move to
  // lower level.
  std::vector< itpp::bvec > AdjustInfoLength_forMLC(const std::vector< itpp::bvec > &infoBits,
                                                    const std::vector< int > &infoLengths_LDPC);

  std::vector< itpp::bvec > AdjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo,
                                                     const std::vector< int > &infoLengths);

  // end of proposals.

  /************************************************************************************
   * DivideIntoLevels -- MLC�p�ɓ��͌n��𕡐��̃��x���ɕ�������B
   * 
   * Arguments:
   *   input -- ���͌n��
   *   numSamples -- �e���x���̏��r�b�g�����i�[
   *
   * Return Value:
   *   Vector_2D< kind > output -- ���x���������ꂽ���n�񂪊e�s�Ɋ��蓖�Ă��Ă���
   ************************************************************************************/
  template< typename kind >
  inline mylib::Vector_2D< kind > DivideIntoLevels(const std::vector< kind > &input, const std::vector< int > &numSamples)
  {
    int sum = mylib::Sum(numSamples);
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
} // end of namespace mylib

#endif
