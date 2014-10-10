/************************************************************************************
 * Path Loss Modelの実装部分
 *   
 * 
 *
 * Contents:
 *   
 *
 * Last Updated: <2014/10/08 20:35:00 from WatanabeYoshito-no-iMac.local by yoshito>
 ************************************************************************************/
#include <iostream>
#include <cmath>
#include "../include/path_loss_model.h"

namespace mylib {

  static const double LIGHT_SPEED = 3e8;
  
  // Simplified Path-Loss Model
  void SimplifiedPathLossModel::Set(double transPower, double d0, double gamma, double frequency)
  {
    transPower_ = transPower;
    d0_ = d0;
    gamma_ = gamma;
    double lambda = LIGHT_SPEED / frequency;

    gain_ = pow(lambda/(4*M_PI*d0_),2);
  }

  double SimplifiedPathLossModel::ReceivedPower(double distance) const
  {
    return transPower_*gain_*pow(d0_/distance,gamma_);
  }

  double SimplifiedPathLossModel::DistanceFromReceivedPower(double pr) const
  {
    return d0_*pow(transPower_ * gain_/pr, 1/gamma_);
  }


  // Free-Space Path-Loss Model
  void FreeSpacePathLossModel::Set(double transPower, double frequency, double gain)
  {
    transPower_ = transPower;
    double lambda = LIGHT_SPEED / frequency;

    pathLoss_ = gain*pow(lambda/4*M_PI,2);
  }

  
  double FreeSpacePathLossModel::ReceivedPowerFromDistance(double distance) const
  {
    return transPower_*pathLoss_/pow(distance, 2);
  }

  double FreeSpacePathLossModel::DistanceFromReceivedPower(double pr) const
  {
    return sqrt(pathLoss_*transPower_/pr);
  }

  /************************************************************************************
   * ShadowFadingModel 
   * 
   * シャドウイングのクラス
   ************************************************************************************/
  double ShadowFadingModel::ReceivedPowerAtDistance(double distance) const
  {
    itpp::Normal_RNG rng(0.0, sigma_dB_*sigma_dB_);
    
    double snr = splm_.ReceivedPower(distance)/splm_.TransPower();
    double snr_dB = itpp::dB(snr);
    
    double receivedSNR_dB = snr_dB + rng();
    double receivedSNR = itpp::inv_dB(receivedSNR_dB);

    return receivedSNR * splm_.TransPower();
  }
  
  double ShadowFadingModel::CoverageRate(double Pmin, double coverageDistance) const
  {
    double a = (Pmin - splm_.ReceivedPower(coverageDistance))/sigma_dB_;
    double b = 10*splm_.Gamma()*std::log10(std::exp(1))/sigma_dB_;

    double first = itpp::Qfunc(a);
    double second = std::exp((2-2*a*b)/(b*b))*itpp::Qfunc((2-a*b)/b);

    return first + second;
  }
  
}
