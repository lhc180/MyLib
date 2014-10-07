/************************************************************************************
 * Path Loss Model
 *   
 * Class of Simplified Path Loss Model (Wireless Communications by Goldsmith, pp.46)
 * 距離と受信パワーの相互変換
 * 
 * Contents:
 *   ReceivedPower
 *
 * Last Updated: <2014/08/29 18:48:48 from WatanabeYoshito-no-iMac.local by yoshito>
 ************************************************************************************/

#ifndef PATH_LOSS_MODEL_H
#define PATH_LOSS_MODEL_H

#include <itpp/itcomm.h>
#include <itpp/itbase.h>

namespace mylib {
  
  class SimplifiedPathLossModel
  {
  private:
    double transPower_;
    double d0_;
    double gamma_;
    double gain_;               // 式(2.41)より求める

  protected:
        
  public:
    SimplifiedPathLossModel():transPower_(0), d0_(0), gamma_(0), gain_(0){ }
    SimplifiedPathLossModel(double transPoewr, double d0, double gamma, double frequency)
    {
      Set(transPoewr, d0, gamma, frequency);
    }
    virtual ~SimplifiedPathLossModel() { }

    void Set(double transPoewr, double d0, double gamma, double frequency);

    double ReceivedPowerFromDistance(double distance) const;
    double DistanceFromReceivedPower(double pr) const;

    double TransPower() const
    {
      return transPower_;
    }
  };

  class FreeSpacePathLossModel
  {
  private:
    double transPower_;
    double pathLoss_;
    
  public:
    FreeSpacePathLossModel();
    FreeSpacePathLossModel(double transPower, double frequency, double gain = 1.0);
    virtual ~FreeSpacePathLossModel()
    { }

    void Set(double transPower, double frequency, double gain = 1.0);
    double ReceivedPowerFromDistance(double distance) const;
    double DistanceFromReceivedPower(double pr) const;
    
  };

  class ShadowFadingModel
  {
  private:
    SimplifiedPathLossModel splm_;
    double sigma_dB_;
    
  public:
    ShadowFadingModel(const SimplifiedPathLossModel &splm, double sigma_dB):
      splm_(splm), sigma_dB_(sigma_dB)
    {
      itpp::RNG_randomize();
    }
    virtual ~ShadowFadingModel()
    { }

    double ReceivedPowerAtDistance(double distance) const;
  };

  
  
}

#endif
