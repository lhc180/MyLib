/************************************************************************************
 * Path Loss Model
 *   
 * Class of Simplified Path Loss Model (Wireless Communications by Goldsmith, pp.46)
 * 距離と受信パワーの相互変換
 * 
 * Contents:
 *   ReceivedPower
 *
 * Last Updated: <2014/03/28 22:37:14 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/

#ifndef PATH_LOSS_MODEL_H
#define PATH_LOSS_MODEL_H

namespace mylib {
  class SimplifiedPathLossModel
  {
  private:
    double d0_;
    double gamma_;
    double gain_;               // 式(2.41)より求める

  protected:
        
  public:
    SimplifiedPathLossModel():d0_(0), gamma_(0), gain_(0){ }
    SimplifiedPathLossModel(double d0, double gamma, double frequency)
    {
      Set(d0, gamma, frequency);
    }
    virtual ~SimplifiedPathLossModel() { }

    void Set(double d0, double gamma, double frequency);

    double ReceivedPowerFromDistance(double distance) const;
    double DistanceFromReceivedPower(double pr) const;
  };

}

#endif
