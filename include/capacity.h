#ifndef CAPACITY_H
#define CAPACITY_H
/**************************************************/
// これはキャパシティ計算用のクラスである
/*************************************************/
#include <itpp/itbase.h>
#include <itpp/itcomm.h>

namespace mylib{
  class Capacity{
  private:
    itpp::Modulator_2D modulator_;
    int bitsPerSymbol_;
    itpp::cvec symbols_;
    int numSymbols_;
    int nTrans_;
    int nTrial_;

  public:
    Capacity(const itpp::Modulator_2D &mod, int nTrans = 1000, int nTrial = 100):
      modulator_(mod), bitsPerSymbol_(mod.bits_per_symbol()), symbols_(mod.get_symbols()), 
      numSymbols_(symbols_.size()), nTrans_(nTrans), nTrial_(nTrial)
    { }

    virtual ~Capacity()
    { }

    // nTrans --- the number of transmitted symbol, nTrial --- trial times 
    double operator()(double N0);
    
  };

  class CodingExponent
  {
  private:
    const itpp::Modulator_2D modulator_;
    const int bitsPerSymbol_;
    const itpp::cvec symbols_;
    const int numSymbols_;
    const int codeLength_;
    const int nTrans_;
    const int nTrial_;

    double CutoffRate(double rho, double n0) const;
    
  public:
    CodingExponent(const itpp::Modulator_2D &mod, int codeLength,
                   int nTrans = 1000, int nTrial = 100):
      modulator_(mod), bitsPerSymbol_(mod.bits_per_symbol()), symbols_(mod.get_symbols()),
      numSymbols_(symbols_.size()), codeLength_(codeLength),
      nTrans_(nTrans), nTrial_(nTrial)
    { }
  
    virtual ~CodingExponent()
    { }

    double CodeRate(double n0, double targetErrorRate, double prevR = 0) const;

    double ErrorRate(double n0, double codeRate) const;
    
  };
  
}

#endif
