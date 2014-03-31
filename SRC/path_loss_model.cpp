/************************************************************************************
 * Path Loss Modelの実装部分
 *   
 * 
 *
 * Contents:
 *   
 *
 * Last Updated: <2014/03/30 02:27:53 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/
#include <iostream>
#include <cmath>
#include "../include/path_loss_model.h"

namespace mylib {

  static const double LIGHT_SPEED = 3e8;
  
  void SimplifiedPathLossModel::Set(double d0, double gamma, double frequency)
  {
    d0_ = d0;
    gamma_ = gamma;
    double lambda = LIGHT_SPEED / frequency;

    gain_ = pow(lambda/(4*M_PI*d0_),2);
  }

  double SimplifiedPathLossModel::ReceivedPowerFromDistance(double distance) const
  {
    return gain_*pow(d0_/distance,gamma_);
  }

  double SimplifiedPathLossModel::DistanceFromReceivedPower(double pr) const
  {
    return d0_*pow(gain_/pr, 1/gamma_);
  }
}
