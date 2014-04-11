/************************************************************************************
 * Path Loss Model
 *   
 * Class of Simplified Path Loss Model (Wireless Communications by Goldsmith, pp.46)
 * 距離と受信パワーの相互変換
 * 
 * Contents:
 *   ReceivedPower
 *
 * Last Updated: <2014/03/31 13:26:58 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/

#ifndef PATH_LOSS_MODEL_H
#define PATH_LOSS_MODEL_H

namespace mylib {

  // 送信電力を1として計算
  class SimplifiedPathLossModel
  {
  private:
    double transPoewr_;
    double d0_;
    double gamma_;
    double gain_;               // 式(2.41)より求める

  protected:
        
  public:
    SimplifiedPathLossModel():transPoewr_(0), d0_(0), gamma_(0), gain_(0){ }
    SimplifiedPathLossModel(double transPoewr, double d0, double gamma, double frequency)
    {
      Set(transPoewr, d0, gamma, frequency);
    }
    virtual ~SimplifiedPathLossModel() { }

    void Set(double transPoewr, double d0, double gamma, double frequency);

    double ReceivedPowerFromDistance(double distance) const;
    double DistanceFromReceivedPower(double pr) const;
  };

}

#endif
