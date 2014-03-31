/************************************************************************************
 * InvFunc
 *   
 * セットされたベクトルの逆関数を求める
 *
 * Contents:
 *   class InvFunc
 *
 * Last Updated: <2014/03/29 17:58:23 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/

#ifndef INV_FUNC_H
#define INV_FUNC_H

#include <cassert>
#include <itpp/itbase.h>

namespace mylib {

  // 単純上昇する関数を想定
  class InvFunc
  {
  private:
    itpp::vec x_;               // x軸
    itpp::vec y_;               // y軸
    int size_;
    double minY_, maxY_;
    
  public:
    InvFunc(){ }
    
    InvFunc(const itpp::vec& x, const itpp::vec& y)
    {
      Set(x, y);
    }
    virtual ~InvFunc() { }

    void Set(const itpp::vec& x, const itpp::vec& y)
    {
      x_ = x;
      y_ = y;
      assert(x_.size() == y_.size());
      size_ = x_.size();
      minY_ = itpp::min(y_);
      maxY_ = itpp::max(y_);
    }
    
    // yに最も近づくx軸の値を変えす
    double X(double y) const
    {
      return x_[IndexX(y)];
    }

    int IndexX(double y) const
    {
      assert(y >= minY_ && y <= maxY_);

      double nearestY = minY_;
      for (int i = 0; i < size_; ++i){
        if (y_[i] > y){
          return i-1;
        } // if
        else{
          nearestY = y_[i];
        } // else 
      } // for i
      return size_ - 1;
    }

  };

  
  
} // namespace mylib

#endif
