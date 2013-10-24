#ifndef QUANTIZER_H
#define QUANTIZER_H

/***********************************/
// スカラー量子化器
//
// 最大振幅値を2^(quantizeSize-1)で割った数が量子化幅
/***********************************/
/*************to do list***********/
// maxAmp_ in normalize should be changed for a positive value.
// only midtread is supported so far. If you add midrise, you have only to change toInt function.
/**********************************/

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <itpp/itbase.h>

namespace mylib{

  class Quantizer{
  protected:
    int quantizeSize_;
    double maxAmp_;
    double div_;

    //   int divNum_;			// 分割数
    //   itpp::vec quantizeVec_;
    bool midrise_;
  
    virtual int ToInt(double input)
    {
  
      double output;
      if(input >=0){
        if(input > maxAmp_ - div_){
          input = maxAmp_ - div_;
        }
        output = input/(maxAmp_);
        output *= pow(2, quantizeSize_ - 1);
      }
      else{
        if(input < -maxAmp_){
          input = -maxAmp_;
        }
        output = input/maxAmp_;
        output *= pow(2, quantizeSize_ - 1);
      }

      int output_i = itpp::round_i(output);

      return output_i;
    }
  
    virtual double ToDouble(int input)
    {
      double output;
      if(input >= 0){
        if(quantizeSize_ == 1){
          output = 0;
        }
        else{
          output = input / pow(2, quantizeSize_ - 1);
          output *= maxAmp_;
        }
      }
      else{
        output = input / pow(2, quantizeSize_ - 1);
        output *= maxAmp_;
      }

      return output;
    }

    // msb last
    virtual itpp::bvec Dec2bin(int input)
    {
      int bintemp = input;
      itpp::bvec output(quantizeSize_);

      for(int i = 0; i < quantizeSize_-1; i++){
        output[i] = itpp::bin((bintemp >> i) & 1);
      }

      // sign bit
      if(input >= 0){
        output[quantizeSize_-1] = 0;
      }
      else{
        output[quantizeSize_-1] = 1;
      }
     
      return output;
    }

    // itpp::bvec dec2bin(const std::vector<int> &input)
    //   {
    //     itpp::bvec output(0);

    //     for(unsigned i = 0; i < input.size(); i++){
    //       output = itpp::concat(output, dec2bin(input[i]));
    //     }
    
    //     return output;
    //   }

    //   itpp::bvec dec2bin(const itpp::ivec &input)
    //   {
    //     itpp::bvec output(0);

    //     for(unsigned i = 0; i < input.size(); i++){
    //       output = itpp::concat(output, dec2bin(input[i]));
    //     }
    
    //     return output;
    //   }
  

    virtual int Bin2dec(const itpp::bvec &input)
    {
      assert(input.size() == quantizeSize_);

      int output = 0;
    
      if(input[quantizeSize_-1]==0){
        for(int i = 0; i < quantizeSize_-1; i++){
          output += pow(2, i) * int(input[i]);
        }
      }
      else{
        for(int i = 0; i < quantizeSize_-1; i++){
          output += pow(2, i) * int(input[i]+itpp::bin(1));
        }
        output = -output;
        output--;
      }
      return output;
    }

   
  public:
    Quantizer(): quantizeSize_(1), maxAmp_(1), midrise_(false) { }
  
    Quantizer(int QuantizeSize, double maxAmp, bool midrise = false)
    {
      Set(QuantizeSize, maxAmp, midrise);
    }
    // デフォルトのコピーコンストラクタ
    // デフォルトデストラクタ

    // 初期化
    // midrise = falseならindex = 0は振幅値0
    virtual void Set(int q, double max, bool midrise = false)
    {
      if(q < 0){
        std::cerr << "Error: Quantize size must be positive.\n";
        exit(1);
      }

      assert(max > 0);

      assert(midrise == false);

      quantizeSize_ = q;
      maxAmp_ = max;
      // maxAmp_ is an absolute of negative max value.
      div_ = maxAmp_/pow(2, quantizeSize_ - 1);

    }


    // return value is the amplitude of a quantized sample.
    virtual double Quantize(double input);
    virtual std::vector<double> Quantize(const std::vector<double> &input);
  
    // return value is the binary signals of quantized samples.
    virtual itpp::bvec Quantize2bvec(double input);
    virtual itpp::bvec Quantize2bvec(const std::vector<double> &input);
    virtual itpp::bvec Quantize2bvec(const itpp::vec &input);
  
    // Dequantize
    virtual double Dequantize(const itpp::bvec &input);
    virtual std::vector<double> Dequantize2stdvec(const itpp::bvec &input);
    virtual itpp::vec Dequantize2itppvec(const itpp::bvec &input);

    virtual int QuantizeSize() const
    {
      return quantizeSize_;
    }
  
  };

  // UnsignedQ
  class UnsignedQ: public Quantizer
  {
  protected:
    virtual int ToInt(double input)
    {
    
      input = (input < 0 ) ? 0 : input;
      input = (input > maxAmp_) ? maxAmp_ : input;

      double output = input/static_cast<double>(maxAmp_);
      output *= (pow(2, quantizeSize_)-1);
      int output_i = itpp::round_i(output);
    
      return output_i;
    }

    virtual double ToDouble(int input)
    {
      double output = input / (pow(2, quantizeSize_)-1);
      output *= maxAmp_;
    
      return output;
    }

    virtual itpp::bvec Dec2bin(int input)
    {
      assert(input >= 0);
    
      int bintemp = input;
      itpp::bvec output(quantizeSize_);

      for(int i = 0; i < quantizeSize_; i++)
        {
          output[i] = itpp::bin((bintemp >> i) & 1);
        }

      return output;
    }

    virtual int Bin2dec(const itpp::bvec &input)
    {
      assert(input.size() == quantizeSize_);

      int output = 0;

      for(int i = 0; i < quantizeSize_; i++)
        {
          output += pow(2, i) * int(input[i]);
        }

      return output;
    }

  public:
    // cUnsignedQ(): cQuantizer()
    // {
    // }

    UnsignedQ(int QuantizeSize, double maxAmp)
    {
      Set(QuantizeSize, maxAmp);
    }

    virtual void Set(int q, double max)
    {
      if(q < 0){
        std::cerr << "Error: Quantize size must be positive.\n";
        exit(1);
      }

      assert(max > 0);

      quantizeSize_ = q;
      maxAmp_ = max;

      div_ = maxAmp_/(pow(2, quantizeSize_)-1);
    }

    virtual double Quantize(double input)
    {
      int int_ized = ToInt(input);
  
      double output = static_cast<double>(int_ized)/(pow(2, quantizeSize_)-1);

      return output * maxAmp_;
    }
  };

  class StepQ: public Quantizer
  {
  protected:
    // 四捨五入(もし振幅がmaxAmp_を超えてたら補正)
    virtual int ToInt(double input)
    {
      if(quantizeSize_ == 0){
        return 0;
      }
      else{
        if (input >= 0 && input > maxAmp_ - 1){
          input = maxAmp_ - 1;
        } // if input
        else if(input < 0 && input < -maxAmp_){
          input = -maxAmp_;
        }
        int output = itpp::round_i(input);
      
        return output;
      }
      
    }
    
  public:
    StepQ(): Quantizer()
    { }
    StepQ(int qSize)
    {
      Set(qSize);
    }

    int MaxAmp() const{
      return maxAmp_ - 1;
    }

    int MinAmp() const{
      return -maxAmp_;
    }
    
    virtual void Set(int qSize)
    {
      if (qSize < 0){
        std::cerr << "Error: Quantize size must be positive." << std::endl;
        exit(1);
      } // if qSize
      quantizeSize_ = qSize;
      maxAmp_ = pow(2, qSize-1);
      div_ = 1;
    }
    
    virtual double Quantize(double input)
    {
      return static_cast< double >(ToInt(input));
    }

    virtual std::vector<double> Quantize(const std::vector<double> &input)
    {
      std::vector< double > output(input.size());
      for (unsigned i = 0; i < input.size(); ++i){
        output[i] = Quantize(input[i]);
      } // for i
      return output;
    }
    
    // MSB first
    virtual itpp::bvec Quantize2bvec(double input)
    {
      int quantized = ToInt(input);
      itpp::bvec out_lsbFirst = Quantizer::Dec2bin(quantized);
      itpp::bvec output = itpp::reverse(out_lsbFirst);
      return output;
    }
    virtual itpp::bvec Quantize2bvec(const std::vector<double> &input)
    {
      return Quantizer::Quantize2bvec(input);
    }
    virtual itpp::bvec Quantize2bvec(const itpp::vec &input)
    {
      return Quantizer::Quantize2bvec(input);
    }
    
    // Assuming MSB first
    virtual double Dequantize(const itpp::bvec &input)
    {
      itpp::bvec reversed = itpp::reverse(input);
      int output = Quantizer::Bin2dec(reversed);
      return output;
    }

    virtual std::vector<double> Dequantize2stdvec(const itpp::bvec &input)
    {
      return Quantizer::Dequantize2stdvec(input);
    }

    virtual itpp::vec Dequantize2itppvec(const itpp::bvec &input)
    {
      return Quantizer::Dequantize2itppvec(input);
    }

  };

  int NecessaryBits(double max, double min);
  
  
}

#endif
