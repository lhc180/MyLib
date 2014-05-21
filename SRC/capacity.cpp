#include <iostream>
#include <cmath>
#include <cassert>
#include <complex>
#include <boost/lexical_cast.hpp>
#include "../include/capacity.h"

namespace mylib{

  inline double PdfAWGN(std::complex<double> y, 
                        std::complex<double> x, double N0)
  {
    return 1.0/sqrt(M_PI*N0)*exp(-pow(std::abs(y-x),2)/N0);
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

        double denominator = PdfAWGN(received(y), symbol(y), N0);

        // logの内側のΣ用
        double sum1 = 0.0;	// log用
        for(int x_dash = 0; x_dash < numSymbols_; x_dash++){
          sum1 += PdfAWGN(received(y), symbols_(x_dash), N0);
        } // for x_dash
      
        sum2 += log2(sum1/denominator);
    
      } // for y
  
      avr += sum2;// static_cast<double>(nTrans);
    }
    avr /= static_cast<double>(nTrans_);

    avr /= static_cast<double>(nTrial_);
  
    return bitsPerSymbol_ - avr;
  }

/************************************************************************************
 * CodingExponent 
 ************************************************************************************/
  double CodingExponent::CutoffRate(double rho, double n0) const
  {
    itpp::AWGN_Channel awgn(n0);
    
    double avr = 0.0;
    // このループがE_{x,y}のxにあたる
    for (int trial = 0; trial < nTrial_; ++trial){
      itpp::bvec bits = itpp::randb(bitsPerSymbol_*nTrans_);
      itpp::cvec transSymbol = modulator_.modulate_bits(bits);
      // itpp::cvec transSymbol = itpp::zeros_c(bitsPerSymbol_*nTrans_);
      itpp::cvec received = awgn(transSymbol);

      double sum2 = 0.0;
      for (int y = 0; y < nTrans_; ++y){

        // double denominator = PdfAWGN(received[y], transSymbol[y], 4*n0);
        
        // logの内側のSigma
        double sum1 = 0;
        double denomSum = 0;
        for (int x_dash = 0; x_dash < numSymbols_; ++x_dash){
          double fY = PdfAWGN(received[y], symbols_[x_dash], n0);
          sum1 += pow(fY, 1.0/(1.0+rho));
          denomSum += PdfAWGN(received[y], symbols_[x_dash], n0);
        } // for x_dash
        
        double numerator = pow(sum1, 1.0+rho);
        double denominator = denomSum / numSymbols_;
        
        // sum2 += log2(numerator/denominator);
        sum2 += (numerator/denominator)/static_cast< double >(nTrans_);
        
      } // for y

      avr += sum2/static_cast< double >(nTrial_);
    } // for trial

      // avr /= static_cast< double >(nTrans_);
    
    double logAvr = log2(avr);
    // std::cout << "##avr = " << avr << std::endl;

    // std::cout << "## cutoff rate = " <<  (1.0+rho)*static_cast< double >(bitsPerSymbol_) - logAvr << std::endl;
    
    return (1.0+rho)*static_cast< double >(bitsPerSymbol_) - logAvr;

  }

  /************************************************************************************
   * CodeRate -- 符号化率を求める
   * 
   * Arguments:
   *   n0 -- 雑音の分散
   *   targetErrorRate -- 達成したいFER
   *   prevR -- この値よりも小さい符号化率は返さない
   *
   * Return Value:
   *   double -- 符号化率
   ************************************************************************************/
  double CodingExponent::CodeRate(double n0, double targetErrorRate, double prevR) const
  {
    static const itpp::vec rho("0:0.01:1");
    std::string codeRateString = boost::lexical_cast< std::string >(prevR) + ":0.001:" +
      boost::lexical_cast< std::string >(bitsPerSymbol_);
    const itpp::vec rate(codeRateString.c_str());

    double targetCodingExponent = -1.0/static_cast< double >(codeLength_)*log2(targetErrorRate);

    itpp::vec E0_rho(rho.size());
    for (int rho_i = 0; rho_i < rho.size(); ++rho_i){
      E0_rho[rho_i] = CutoffRate(rho[rho_i], n0);
    } // for rho_i
    
    double finalR = prevR;
    for (int rate_i = 0; rate_i < rate.size(); ++rate_i){
      double Er_max = std::numeric_limits< double >::min();

      for (int rho_i = 0; rho_i < rho.size(); ++rho_i){
        double Er_temp = E0_rho[rho_i] - rho[rho_i]*rate[rate_i];
        if (Er_temp > Er_max){
          Er_max = Er_temp;
        } // if 
      } // for rho_i

      if (targetCodingExponent <= Er_max){
        finalR = rate[rate_i];
      } // if
      else{
        break;
      } // else 
    } // for r_i

    return finalR;
  }

  double CodingExponent::ErrorRate(double n0, double codeRate) const
  {
    static const itpp::vec rho("0:0.01:1");
    double finalErrorRate = 1.0;
    double errorRate = 1.0;

    itpp::vec E0_rho(rho.size());
    for (int rho_i = 0; rho_i < rho.size(); ++rho_i){
      E0_rho[rho_i] = CutoffRate(rho[rho_i], n0);
    } // for rho_i
    
    while (errorRate >= 1e-6){
      double codingExponent = -1.0/static_cast< double >(codeLength_)*log2(errorRate);
      double Er_max = std::numeric_limits< double >::min();

      for (int rho_i = 0; rho_i < rho.size(); ++rho_i){
        double Er_temp = E0_rho[rho_i] - rho[rho_i]*codeRate;
        if (Er_temp > Er_max){
          Er_max = Er_temp;
        } // if 
      } // for rho_i

      if (codingExponent <= Er_max){
        finalErrorRate = errorRate;
        errorRate /= 2.0;
      } // if
      else{
        break;
      } // else
    } // while 

    return finalErrorRate;
  }
  
  
} // namespace mylib

