/************************************************************************************
 * mymodulation.cpp
 *   
 * 色々なラベリング手法
 *
 * Contents:
 *   HybridPartitioningPSK
 *
 * Last Updated: <2014/03/24 23:49:05 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/

#include <cassert>
#include "../include/mymodulation.h"

namespace mylib {

  // ++++++++ 8PSK用のAsymmetric Constellation
  static itpp::cvec AsymmetricConstellation8PSK(double angle)
  {
    itpp::cvec symbols(8);
    
    double pi_2 = M_PI / 2.0;   // pi/2
    double epsilon = pi_2/10000.0; // 誤差補正用
    
    for (int i = 0; i < 4; ++i){
      std::complex< double > symbol = std::complex<double>(std::polar(1.0, pi_2*i));
      if (std::fabs(std::real(symbol)) < epsilon) {
        symbols(2*i) = std::complex<double>(0.0, std::imag(symbol));
      }
      else if (std::fabs(std::imag(symbol)) < epsilon) {
        symbols(2*i) = std::complex<double>(std::real(symbol), 0.0);
      }
      else {
        symbols(2*i) = symbol;
      }
      
      symbols[2*i+1] = std::polar(1.0, pi_2*i + angle); // 補正無し
    } // for i

    return symbols;
  }
    
  /****************** class NaturalPSK ************************/
  void SetPartitioningPSK::Init()
  {
    itpp::ivec newBitmap(M);
    
    for (int i = 0; i < M; ++i){
      // std::cout << "## i = " << i << std::endl;
      // std::cout << "## i_reversed = " << i_reversed << std::endl;      
      newBitmap[i] = itpp::reverse_int(k, i);
    } // for i

    itpp::cvec newSymbols = symbols;

    set(newSymbols, newBitmap);
  }
  /******************* end of NaturalPSK **************************/

  void ReverseSpPSK::Init()
  {
    itpp::ivec newBitmap(M);

    for (int i = 0; i < M; ++i){
      newBitmap[i] = i;
    } // for i

    itpp::cvec newSymbols = symbols;

    set(newSymbols, newBitmap);
  }

  void BlockPartitioningPSK::Init(double angle)
  {
    assert(M == 8);
    itpp::ivec newBits2Symbol = bits2symbols;
    itpp::cvec newSymbols = AsymmetricConstellation8PSK(angle);

    set(newSymbols, newBits2Symbol);
  }
  
  void HybridPartitioningPSK_1::Init(double angle)
  {
    assert(M == 8);             // まだ8PSKにしか対応していない
    itpp::ivec newBits2Symbol = "0 6 1 7 3 5 2 4";
    itpp::cvec newSymbols = AsymmetricConstellation8PSK(angle);
    
    set(newSymbols, newBits2Symbol);
  }

  void HybridPartitioningPSK_2::Init(double angle)
  {
    assert(M == 8);
    itpp::ivec newBits2Symbol = "0 1 4 5 3 2 7 6";
    itpp::cvec newSymbols = AsymmetricConstellation8PSK(angle);

    set(newSymbols, newBits2Symbol);
  }

  // void BlockPartitioning16APSK::Init()
  // {
    
  // }
  
}
