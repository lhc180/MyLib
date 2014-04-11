/************************************************************************************
 * Path Loss Modelの実装部分
 *   
 * 
 *
 * Contents:
 *   
 *
 * Last Updated: <2014/03/31 13:30:05 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/
#include <iostream>
#include <cmath>
#include "../include/path_loss_model.h"

namespace mylib {

  static const double LIGHT_SPEED = 3e8;
  
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
}
