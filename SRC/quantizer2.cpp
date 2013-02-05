#include <cstdlib>
#include <assert.h>
#include "../include/quantizer.h"

/************to do list*************/
// quantizeはもっと高速化できる
/**********************************/

void cQuantizer::set(int q, double max, bool midrise)
{
  if(q < 0){
    std::cerr << "Error: Quantize size must be positive.\n";
    exit(1);
  }

  quantizeSize = q;
  maxAmp = max;
  divNum = static_cast<int>(pow(2,quantizeSize));
  quantizeVec.set_size(divNum);

  double divAmp = maxAmp / (static_cast<double>(divNum)/2.0);

  if(q == 0){
    quantizeVec[0] = 0.0;
  }
  else{
    if(midrise){
      double divAmp = 2.0*maxAmp/static_cast<double>(divNum-1);
      for(int i = 0; i < divNum/2; i++){
        quantizeVec[i] = divAmp/2.0 + divAmp*i;
      }
    
      int j = 1;
      for(int i = divNum-1; i >= divNum/2; i--,j++){
        quantizeVec[i] = -divAmp/2.0 - divAmp*j;
      }
    }
    else{
      for(int i = 0; i < divNum/2; i++){
        quantizeVec[i] = divAmp*i;
      }
      int j = 1;
      for(int i = divNum-1; i >= divNum/2; i--,j++){
        quantizeVec[i] = -1.0*divAmp*j;
      }
    }
    // std::cout << "## quantizeVec = " << quantizeVec << std::endl;
  }
}

double cQuantizer::quantize(const double input)
{
  int index = 0;
  double minimumError = 100.0;
  double error;
  
  for(int i = 0; i < quantizeVec.size(); i++){
    error = pow(input - quantizeVec[i],2);
    if(error < minimumError){
      index = i;
      minimumError = error;
    }
  }

  return quantizeVec[index];
}
  
std::vector<double> cQuantizer::quantize(const std::vector<double> &input)
{
  std::vector<double> output(input.size());

  for(int i = 0; i < static_cast<int>(input.size()); i++){
    output[i] = quantize(input[i]);
  }

  return output;
}

int cQuantizer::quantizeReturnIndex(const double input)
{
  int index = 0;
  double minimumError = 100.0;
  double error;
  
  for(int i = 0; i < quantizeVec.size(); i++){
    error = pow(input - quantizeVec[i],2);
    if(error < minimumError){
      index = i;
      minimumError = error;
    }
  }

  return index;
}

std::vector<int> cQuantizer::quantizeReturnIndex(const std::vector<double> &input)
{
  std::vector<int> output(input.size());

  for(int i = 0; i < static_cast<int>(input.size()); i++){
    output[i] = quantizeReturnIndex(input[i]);
  }

  return output;
}

double cQuantizer::dequantize(const int index)
{
  assert((index >= 0) && (index < divNum));

  return quantizeVec[index];
}

std::vector<double> cQuantizer::dequantize(const std::vector<int> &index)
{
  std::vector<double> output(index.size());

  for(int i = 0; i < static_cast<int>(index.size()); i++){
    output[i] = dequantize(index[i]);
  }

  return output;
}

int cQuantizer::getQuantizeSize()
{
  return quantizeSize;
}
