#include <iostream>
#include <cmath>
#include <cassert>
#include <complex>
#include "../include/capacity.h"

namespace mylib{
    
  // 条件付き確率
  inline double pdf_awgn(std::complex<double> y, std::complex<double> x, double N0)
  {
    // std::cout << "abs(y-x) = " << abs(y-x) << std::endl;

    return exp(-pow(std::abs(y-x),2)/N0);
  }

  // キャパシティを返す
  double Capacity::operator()(double N0)
  {
    itpp::AWGN_Channel awgn(N0);

    double avr = 0.0;
    for(int trial = 0; trial < nTrial_; trial++){

      itpp::bvec bits = itpp::randb(bitsPerSymbol_*nTrans_);

      // nTrial個のシンボル
      itpp::cvec symbol = modulator_.modulate_bits(bits);

      itpp::cvec received = awgn(symbol);  

      // 受信信号の数
      double sum2 = 0.0;		// Ex用 
      for(int y = 0; y < nTrans_; y++){ 
        // コンスタレーションの数
    
        // for(int x = 0; x < nSymbols; x++){
  
        double denominator = pdf_awgn(received(y), symbol(y), N0);

        // logの内側のΣ用
        double sum1 = 0.0;	// log用
        for(int x_dash = 0; x_dash < numSymbols_; x_dash++){
          sum1 += pdf_awgn(received(y), symbols_(x_dash), N0);
        }	// for x_dash
      
        sum2 += log2(sum1/denominator);
    
      } // for y
  
      avr += sum2;// static_cast<double>(nTrans);
    }
    avr /= static_cast<double>(nTrans_);

    avr /= static_cast<double>(nTrial_);
  
    return bitsPerSymbol_ - avr;
  }
  
} // namespace mylib

