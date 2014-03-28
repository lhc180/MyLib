/************************************************************************************
 * mymodulation.cpp
 *   
 * 色々なラベリング手法
 *
 * Contents:
 *   HybridPartitioningPSK
 *
 * Last Updated: <2014/03/26 00:14:22 from Yoshitos-iMac.local by yoshito>
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

  static itpp::cvec AsymmetricConstellationQAM(int mAry, double dRatio)
  {
    int k = itpp::levels2bits(mAry);
    assert((itpp::pow2i(k) == mAry) && (itpp::is_even(k)));

    int num1side = itpp::round_i(std::sqrt(static_cast< double >(mAry)));
    double signalDistance = dRatio - 1.0;

    itpp::vec signalInSide(num1side/2);
    for (int i = 0; i < num1side/2; ++i){
      signalInSide[i] = 1 + signalDistance*i;
    } // for i
    double averageEnergy = itpp::sum_sqr(signalInSide);
    double scalingFactor = std::sqrt(averageEnergy);

    itpp::cvec symbols(mAry);

    // 正の虚数部分
    for (int i = 0; i < num1side/2; ++i){
      for (int j = 0; j < num1side/2; ++j){
        symbols[i*num1side + j] = std::complex<double>((1+(num1side/2-1-j)*signalDistance)/scalingFactor,
                                                       (1+(num1side/2-1-i)*signalDistance)/scalingFactor);
      } // for j
      for (int j = num1side/2; j < num1side; ++j){
        symbols[i*num1side + j] = std::complex<double>((-1-(j-num1side/2)*signalDistance)/scalingFactor,
                                                       (1+(num1side/2-1-i)*signalDistance)/scalingFactor);
      } // for j
    } // for i
    // 負の実数部分
    for (int i = num1side/2; i < num1side; ++i){
      for (int j = 0; j < num1side/2; ++j){
        symbols[i*num1side + j] = std::complex<double>((1+(num1side/2-1-j)*signalDistance)/scalingFactor,
                                                       (-1-(i-num1side/2)*signalDistance)/scalingFactor);
      } // for j
      for (int j = num1side/2; j < num1side; ++j){
        symbols[i*num1side + j] = std::complex<double>((-1-(j-num1side/2)*signalDistance)/scalingFactor,
                                                       (-1-(i-num1side/2)*signalDistance)/scalingFactor);
      } // for j
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

  void BlockPartitioningQAM::Init(int m, double dRatio)
  {
    assert(m == 16);          // まだ16QAMにしか対応していない

    itpp::ivec newBits2Symbol = "0 1 4 5 2 3 6 7 8 9 12 13 10 11 14 15";

    itpp::cvec newSymbols = AsymmetricConstellationQAM(m, dRatio);
    
    set(newSymbols, newBits2Symbol);
      
  }

  
  void SemiBlockPartitioningQAM::Init(int m, double dRatio)
  {
    assert(m == 16);

    itpp::ivec newBits2Symbol = "0 4 1 2 3 7 5 6 15 11 14 13 12 8 10 9";

    itpp::cvec newSymbols = AsymmetricConstellationQAM(m, dRatio);
    
    set(newSymbols, newBits2Symbol);
  }
  
  // void BlockPartitioning16APSK::Init()
  // {
    
  // }
  
}
