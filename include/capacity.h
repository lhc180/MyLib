#ifndef CAPACITY_H
#define CAPACITY_H
/**************************************************/
// これはキャパシティ計算用のクラスである
/*************************************************/
#include <itpp/itbase.h>
#include <itpp/itcomm.h>

class cCapacity{
private:
  itpp::Modulator_2D Modulator;
  int nBitsPerSymbol;
  int nSymbols;
  itpp::cvec vecSymbols;
  bool bSetMod;

public:
  cCapacity(): nBitsPerSymbol(0), nSymbols(0), vecSymbols(0), bSetMod(false) 
  { }
  cCapacity(itpp::Modulator_2D &Mod)
  {
    setModulator(Mod);
  }

  // ~cCapacity() --- デフォルトデストラクタ

  void setModulator(itpp::Modulator_2D &Mod );

  // nTrans --- the number of transmitted symbol, nTrial --- trial times 
  double getCapacity(double EsN0dB, int nTrans = 1000, int nTrial = 100);
};


#endif
