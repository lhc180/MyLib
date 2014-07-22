/************************************************************************************
 * Path Loss Modelの実装部分
 *   
 * 
 *
 * Contents:
 *   
 *
 * Last Updated: <2014/07/22 21:22:56 from WatanabeYoshito-no-iMac.local by yoshito>
 ************************************************************************************/
#include <iostream>
#include <cmath>
#include "../include/path_loss_model.h"

namespace mylib {

  static const double LIGHT_SPEED = 3e8;
  
  // Simplified Path-Loss Model
  void SimplifiedPathLossModel::Set(double transPoewr, double d0, double gamma, double frequency)
  {
    transPoewr_ = transPoewr;
    d0_ = d0;
    gamma_ = gamma;
    double lambda = LIGHT_SPEED / frequency;

    gain_ = pow(lambda/(4*M_PI*d0_),2);
  }

  double SimplifiedPathLossModel::ReceivedPowerFromDistance(double distance) const
  {
    return transPoewr_*gain_*pow(d0_/distance,gamma_);
  }

  double SimplifiedPathLossModel::DistanceFromReceivedPower(double pr) const
  {
    return d0_*pow(transPoewr_ * gain_/pr, 1/gamma_);
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
}
