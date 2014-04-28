#ifndef CAPACITY_H
#define CAPACITY_H
/**************************************************/
// ����̓L���p�V�e�B�v�Z�p�̃N���X�ł���
/*************************************************/
#include <itpp/itbase.h>
#include <itpp/itcomm.h>

namespace mylib{
  class Capacity{
  private:
    itpp::Modulator_2D modulator_;
    int bitsPerSymbol_;
    int numSymbols_;
    itpp::cvec symbols_;
    bool setMod_;

  public:
    Capacity(): bitsPerSymbol_(0), numSymbols_(0), symbols_(0), setMod_(false) 
    { }
    Capacity(const itpp::Modulator_2D &Mod)
    {
      SetModulator(Mod);
    }

    // ~Capacity() --- �f�t�H���g�f�X�g���N�^

    void SetModulator(const itpp::Modulator_2D &Mod );

    // nTrans --- the number of transmitted symbol, nTrial --- trial times 
    double GetCapacity(double N0, int nTrans = 1000, int nTrial = 100);
    
  };
  
}

#endif
