#include <iostream>
#include <cmath>
#include <cassert>
#include <complex>
#include "../include/capacity.h"

void cCapacity::setModulator(itpp::Modulator_2D &Mod)
{
  Modulator = Mod;
  nBitsPerSymbol = Modulator.bits_per_symbol();
  vecSymbols = Modulator.get_symbols();
  
  // std::cout << "## vecSymbols = " << vecSymbols << std::endl;
  
  nSymbols = vecSymbols.size();

  bSetMod = true;
}
  
// 条件付き確率
inline double probability(std::complex<double> y, std::complex<double> x, double N0)
{
  // std::cout << "abs(y-x) = " << abs(y-x) << std::endl;

  return exp(-pow(std::abs(y-x),2)/N0);
}

// キャパシティを返す
double cCapacity::getCapacity(double EsN0dB, int nTrans, int nTrial)
{
  assert(bSetMod);

  
  double Es = 1.0;
  double EsN0 = pow(10.0, EsN0dB/10.0);
  double N0 = Es/EsN0;

  itpp::AWGN_Channel awgn(N0);

  

  double avr = 0.0;
  for(int trial = 0; trial < nTrial; trial++){

    itpp::bvec bits = itpp::randb(nBitsPerSymbol*nTrans);

    // nTrial個のシンボル
    itpp::cvec symbol = Modulator.modulate_bits(bits);
    assert(symbol.size() == nTrans);

    itpp::cvec received = awgn(symbol);  

    // 受信信号の数
    double sum2 = 0.0;		// Ex用 
    for(int y = 0; y < nTrans; y++){ 
      // コンスタレーションの数
    
      // for(int x = 0; x < nSymbols; x++){
  
      double denominator = probability(received(y), symbol(y), N0);

      // logの内側のΣ用
      double sum1 = 0.0;	// log用
      for(int x_dash = 0; x_dash < nSymbols; x_dash++){
	sum1 += probability(received(y), vecSymbols(x_dash), N0);
      }	// for x_dash
      
      sum2 += log2(sum1/denominator);
    
    } // for y
  
    avr += sum2;// static_cast<double>(nTrans);
  }
  avr /= static_cast<double>(nTrans);

  avr /= static_cast<double>(nTrial);
  
  return nBitsPerSymbol - avr;
}
  
