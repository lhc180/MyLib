#include <cstdlib>
#include <assert.h>
#include "../include/quantizer.h"

/************to do list*************/
// quantizeはもっと高速化できる
/**********************************/



double cQuantizer::quantize(const double input)
{
  int int_ized = toInt(input);
  
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

 
std::vector<double> cQuantizer::quantize(const std::vector<double> &input)
{
  std::vector<double> output(input.size());

  for(int i = 0; i < static_cast<int>(input.size()); i++){
    output[i] = quantize(input[i]);
  }

  return output;
}

itpp::bvec cQuantizer::quantize2bvec(const double input)
{
  int output_i = toInt(input);

  itpp::bvec output(0);

  output = dec2bin(output_i);

  return output;
}

itpp::bvec cQuantizer::quantize2bvec(const std::vector<double> &input)
{
  itpp::bvec output(0);

  for(int i = 0; i < static_cast<int>(input.size()); i++){
    output = itpp::concat(output, quantize2bvec(input[i]));
  }

  return output;
}

itpp::bvec cQuantizer::quantize2bvec(const itpp::vec &input)
{
  itpp::bvec output(0);

  for(int i = 0; i < static_cast<int>(input.size()); i++){
    output = itpp::concat(output, quantize2bvec(input[i]));
  }

  return output;
}

double cQuantizer::dequantize(const itpp::bvec &input)
{
  assert(input.size() == quantizeSize_);

  int output_i = bin2dec(input);

  double output = toDouble(output_i);
  
  return output;
}

std::vector<double> cQuantizer::dequantize2stdvec(const itpp::bvec &input)
{
  assert(input.size()%quantizeSize_ == 0);

  std::vector<double> output(input.size()/quantizeSize_);

  itpp::bvec temp;
  for(int i = 0; i < static_cast<int>(output.size()); i++){
    temp = input.mid(i*quantizeSize_, quantizeSize_);
    output[i] = dequantize(temp);
  }

  return output;
}

itpp::vec cQuantizer::dequantize2itppvec(const itpp::bvec &input)
{
  assert(input.size()%quantizeSize_ == 0);

  itpp::vec output(input.size()/quantizeSize_);

  itpp::bvec temp;
  for(int i = 0; i < static_cast<int>(output.size()); i++){
    temp = input.mid(i*quantizeSize_, quantizeSize_);
    output[i] = dequantize(temp);
  }

  return output;
}

