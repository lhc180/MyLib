#ifndef MYCOMM_H
#define MYCOMM_H

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <iostream>
#include <ctime>

#include "myldpc.h"
#include "mlc_msd.h"
#include "capacity.h"
#include "convolutional_code.h"
#include "rcpc_theoretical.h"
#include "mybinvec.h"
#include "mybinmat.h"
#include "length_adjuster.h"
#include "mycontainer_mlc.h"

enum MODULATOR_TYPE{TYPE_BPSK, TYPE_QPSK, TYPE_8PSK, TYPE_16PSK, TYPE_16QAM};

namespace mylib{
  double erfc_inv(double x);      // 逆誤差補関数
  void makeInterleaver(itpp::ivec &interleaver);
  itpp::bvec GenErrors(const itpp::bvec &input, double errorRate);

  itpp::bvec std_ivec2bvec(int length, const std::vector<int> &input);
  itpp::bvec itpp_ivec2bvec(int length, const itpp::ivec &input);
  std::vector<int> std_bvec2ivec(int length, const itpp::bvec &input);
  itpp::ivec itpp_bvec2ivec(int length, const itpp::bvec &input);


  // +++++++++++++++++++ 与えられたBERからある変調方式のEbN0を返す +++++++++++++++++++
  class BERtoEbN0
  {
  protected:
    MODULATOR_TYPE modulatorType_;
    bool setDone_;

  public:
    BERtoEbN0(): setDone_(false){ }
    BERtoEbN0(MODULATOR_TYPE modType)
    {
      setModulatorType(modType);
    }

    void setModulatorType(MODULATOR_TYPE modType)
    {
      modulatorType_ = modType;
    
      setDone_ = true;
    }
            
    double calcEbN0(double ber)
    {
      double EbN0;
    
      if(modulatorType_ == TYPE_BPSK){
        EbN0 = pow(mylib::erfc_inv(2.0*ber),2);
      }
      else{
        std::cout << "The modulator type can not be available for now.\n";
        exit(1);
      }
    
      return EbN0;
    }
  };

  // +++++++++++++++++++ 与えられたレートで1を返す +++++++++++++++++++++
  class Rand01{
  private: 
    double rate_;
  
  public:
    Rand01(double rate = 0){
      setRate(rate);
    }
  
    void setRate(double rate = 0){
      assert(rate >= 0 && rate <= 1);
      rate_ = rate;
    }

    itpp::bin get(){
      double rnd;
      rnd = static_cast<double>(rand()+1)/(static_cast<double>(RAND_MAX)+1);
    
      itpp::bin out(0);
      if(rnd <= rate_){
        out = 1;
      }

      return out;
    }
  
    double getRate(){
      return rate_;
    }
  };

  // +++++++++++++++++ 特定数の1と残りは0で埋めた配列を返す +++++++++++++
  // +++++++++++++++++ 毎回ランダムでインターリーブされる  +++++++++++++
  class RandomArrange1vec
  {
  protected:
    int allNum_;                     // 要素数
    int num1_;
    bool setDone_;

  public:
    RandomArrange1vec() : setDone_(false){ }
    RandomArrange1vec(int allNum, int num1)
    {
      set(allNum, num1);
    }

    void set(int allNum, int num1){
      assert(allNum >= num1);
      allNum_ = allNum;
      num1_ = num1;
      setDone_ = true;
    }
  
    itpp::bvec get()
    {
      assert(setDone_);
    
      itpp::bvec temp(allNum_);
      temp.zeros();

      for(int i = 0; i < num1_; i++){
        temp[i] = 1;
      }

      itpp::Sequence_Interleaver<itpp::bin> interleaver(allNum_);
      interleaver.randomize_interleaver_sequence();

      itpp::bvec output = interleaver.interleave(temp);

      assert(temp.size() == output.size());
    
      return output;
    }
  };

  // +++++++++++++++ イベント確率を求めるためのクラス ++++++++++++++++
  // +++++++++++++++ 今のところitpp::binにしか対応してない +++++++++++
  class BinCounter
  {
  protected:
    unsigned trialsCounter_;
    unsigned count1_;
    
  public:
    BinCounter()
    {
      reset();
    }

    // default copy constructor
    // default destructor
  
    void reset()
    {
      trialsCounter_ = 0;
      count1_ = 0;
    }
  
    void update(const itpp::bin &input)
    {
      if(input == itpp::bin(1)) {
        count1_++;
      }
      trialsCounter_++;
    }

    void update(const itpp::bvec &input)
    {
      for (int i = 0, size = input.size(); i < size; ++i){
        update(input[i]);
      } // for i
    }
    
    double GetRate(const itpp::bin &input)
    {
      if(input == itpp::bin(1))
        {
          return static_cast<double>(count1_)/static_cast<double>(trialsCounter_);
        }
      else
        {
          unsigned count0 = trialsCounter_ - count1_;
          return static_cast<double>(count0)/static_cast<double>(trialsCounter_);
        }
    }
      
  };

  // インターリーバを作る関数
  inline void makeInterleaver(itpp::ivec &interleaver)
  {
    itpp::ivec temp(interleaver.size());

    for(int i = 0; i < temp.size(); i++){
      temp[i] = i;
    }
  
  
    for(int i = 0; i < interleaver.size(); i++){
      int num = rand()%temp.size();
      interleaver[i] = temp[num];
      temp.del(num);
    }
  }

  inline double erfc_inv(double x) 
  { 
    double Z,W,WI,F,D,Z2; 
  
    if (x>1.0) { 
      return -erfc_inv(2.0-x); 
    } 
    Z=1.0-x; 
    if (Z>0.85) { 
      W=sqrt(-log(x+x*Z)); 
      if (W>=2.5) { 
        if (W>=4.0) { 
          WI=1.0/W; 
          F=0.01078639*WI-0.1498384; 
          F=(F*WI-0.002028152)*WI; 
          D=WI-0.06888301; 
          D=D*WI+0.5211733; 
          D=D*WI+0.09952975; 
          return W+W*(0.1851159E-3+F/D); 
        } else { 
          F=0.06208963*W-0.3166501; 
          F=(F*W+0.3937021)*W; 
          D=W-2.962883; 
          D=D*W+4.666263; 
          D=D*W-6.266786; 
          return W+W*(-0.05668422+F/D); 
        } 
      } else { 
        F=0.05073975*W-0.2368201;
        F=(F*W-0.1314774)*W; 
        D=(W-7.586103)*W+21.98546; 
        D=D*W-44.27977; 
        return W+W*(-0.1146666+F/D); 
      } 
    } else { 
      Z2=Z*Z; 
      F=-1.187515+Z2; 
      F=-2.374996+Z2-0.05496261/F; 
      F=-3.293474+Z2-1.896513/F; 
      F=-0.1137730-0.5751703*Z2/F; 
      return Z+Z*F; 
    } 
  } 
     
  inline itpp::bvec GenErrors(const itpp::bvec &input, double errorProb)
  {
    itpp::bvec output(input.size());
    Rand01 Rand01(errorProb);

    for(int i = 0; i < input.size(); i++){
      output[i] = input[i] + Rand01.get();
    }

    return output;
  }
  
  inline itpp::bvec std_ivec2bvec(int length, const std::vector<int> &input)
  {
    itpp::bvec output(0);
  
    for(int i = 0; i < static_cast<int>(input.size()); i++){
      itpp::bvec t_bvec = itpp::dec2bin(length, input[i]);
      output = itpp::concat(output, t_bvec);
    }

    return output;

  }

  inline itpp::bvec itpp_ivec2bvec(int length, const itpp::ivec &input)
  {
    itpp::bvec output(0);
  
    for(int i = 0; i < input.size(); i++){
      itpp::bvec t_bvec = itpp::dec2bin(length, input[i]);
      output = itpp::concat(output, t_bvec);
    }

    return output;

  }

  inline std::vector<int> std_bvec2ivec(int length, const itpp::bvec &input)
  {
    std::vector<int> output(0);

    for(int start = 0; start < input.size(); start += length){
      itpp::bvec t_bvec = input.mid(start, length);
      int temp = itpp::bin2dec(t_bvec);
      output.push_back(temp);
    }

    return output;
  }

  inline itpp::ivec itpp_bvec2ivec(int length, const itpp::bvec &input)
  {
    itpp::ivec output(0);
    for(int start = 0; start < input.size(); start += length){
      itpp::bvec t_bvec = input.mid(start, length);
      int temp = itpp::bin2dec(t_bvec);
      output = itpp::concat(output, temp);
    }

    return output;
  }

}

#endif
