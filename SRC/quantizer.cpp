#include <cstdlib>
#include <assert.h>
#include "../include/quantizer.h"

/************to do list*************/
// quantizeはもっと高速化できる
/**********************************/

namespace mylib{

  inline double Quantizer::Quantize(const double input)
  {
    int int_ized = ToInt(input);
  
    double output;
    if(int_ized >= 0){
      if(quantizeSize_ == 1){
        output = 0;
      }
      else{
        output = static_cast<double>(int_ized)/(pow(2, quantizeSize_ - 1));
      }
    }
    else{
      output = static_cast<double>(int_ized)/(pow(2, quantizeSize_ - 1));
    }

    return output*maxAmp_;
  }

 
  inline std::vector<double> Quantizer::Quantize(const std::vector<double> &input)
  {
    std::vector<double> output(input.size());

    int i = 0;
    for(std::vector< double >::const_iterator curInput = input.begin();
        curInput != input.end(); ++curInput, ++i){
      output[i] = Quantize(*curInput);
    }

    return output;
  }

  inline itpp::bvec Quantizer::Quantize2bvec(const double input)
  {
    int output_i = ToInt(input);

    itpp::bvec output(0);

    output = Dec2bin(output_i);

    return output;
  }

  inline itpp::bvec Quantizer::Quantize2bvec(const std::vector<double> &input)
  {
    itpp::bvec output(0);

    for(int i = 0; i < static_cast<int>(input.size()); i++){
      output = itpp::concat(output, Quantize2bvec(input[i]));
    }

    return output;
  }

  inline itpp::bvec Quantizer::Quantize2bvec(const itpp::vec &input)
  {
    itpp::bvec output(0);

    for(int i = 0; i < static_cast<int>(input.size()); i++){
      output = itpp::concat(output, Quantize2bvec(input[i]));
    }

    return output;
  }

  inline double Quantizer::Dequantize(const itpp::bvec &input)
  {
    assert(input.size() == quantizeSize_);

    int output_i = bin2dec(input);

    double output = ToDouble(output_i);
  
    return output;
  }

  inline std::vector<double> Quantizer::Dequantize2stdvec(const itpp::bvec &input)
  {
    assert(input.size()%quantizeSize_ == 0);

    std::vector<double> output(input.size()/quantizeSize_);

    itpp::bvec temp;
    for(int i = 0; i < static_cast<int>(output.size()); i++){
      temp = input.mid(i*quantizeSize_, quantizeSize_);
      output[i] = Dequantize(temp);
    }

    return output;
  }

  inline itpp::vec Quantizer::Dequantize2itppvec(const itpp::bvec &input)
  {
    assert(input.size()%quantizeSize_ == 0);

    itpp::vec output(input.size()/quantizeSize_);

    itpp::bvec temp;
    for(int i = 0; i < static_cast<int>(output.size()); i++){
      temp = input.mid(i*quantizeSize_, quantizeSize_);
      output[i] = Dequantize(temp);
    }

    return output;
  }

  int NecessaryBits(double max, double min)
  {
    int qSize = 0;
    double amp = pow(2, qSize-1);
        
    while (amp-1 < max || -amp > min){
      ++qSize;
      amp = pow(2, qSize-1);
    } // while
    return qSize;
  }
  
}
