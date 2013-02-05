#ifndef UTILITY2_H
#define UTILITY2_H

#include<iostream>
#include<vector>
#include<fstream>
#include <cmath>
#include <cassert>

/***************************/
// itppが使えない環境用の各種変換用関数
// また、その他の変換用の関数など
//
// itppは使わない
/****************************/
// cppファイルを作って分割コンパイルしたが、上手くリンクできないみたいなので、
// とりあえずはヘッダファイルに関数の定義も書いておく

namespace mylib
{

  // int型のvectorをnormで正規化する
  void normalizeData(const std::vector<int> &input, std::vector<double> &output, unsigned norm);

  // fileにvectorを出力するが、withIndexがtrueならインデックスも付加する
  template<typename kind>
  void fprintVec(std::ofstream &file,
		  const std::vector<kind> &input,
		  bool withIndex = true);

  // vectorを0で初期化する
  template<typename type>
  void vecZeros(std::vector<type> &input);


  // 10進数から2進数に変換
  // MLBから出力
  std::vector<int> dec2bin(int length, int index);

  // 2進数から10進数に変換
  int bin2dec(const std::vector<int> &inbvec, bool msb_first = false);

  // inputのstart番目の要素からnumber分抜き出す
  template<typename kind>
  std::vector<kind> getMid(const std::vector<kind> &input, const int start, const int number);

  // forwardVecにbackwardVecを結合する
  template<typename kind>
  std::vector<kind> concat(const std::vector<kind> &forwardVec, const std::vector<kind> &backwardVec);

  template<typename type>
  type maxOf(type x1, type x2);

  template<typename type>
  type minOf(type x1, type x2);

  int DoubleToShort(double input);
  std::vector<int> DoubleToShort(const std::vector<double> &input);

  //小数点以下を四捨五入
  int round_uPoint(double input);

  // 符号
  int sgn(double in);
} // end of namespace mylib

template<typename type>
inline void mylib::vecZeros(std::vector<type> &input)
{
  for(int i = 0; i < static_cast<int>(input.size()); i++){
    input[i] = 0;
  }
}

inline void mylib::normalizeData(const std::vector<int> &input, std::vector<double> &output, unsigned norm)
{
  for(unsigned i = 0; i < input.size(); i++){
    output[i] = static_cast<double>(input[i])/static_cast<double>(norm);
  }
}

template<typename kind> 
void mylib::fprintVec(std::ofstream &file, 
		      const std::vector<kind> &input,
		      bool withIndex = true)
{
  if(withIndex){
    for(int i = 0; i < static_cast<int>(input.size()); i++){
      file << i << '\t' << input[i] << '\n';
    }
  }
  else{
    for(int i = 0; i < static_cast<int>(input.size()); i++){
      file << input[i] << '\n';
    }
  }
}

  // MLBから出力
inline std::vector<int> mylib::dec2bin(int length, int index)
{
  int i, bintemp = index;
  std::vector<int> temp(length);
 
  for (i = 0; i < length; i++) {
    temp[i] = bintemp & 1;
    bintemp = (bintemp >> 1);
  }
  return temp;
}

inline int mylib::bin2dec(const std::vector<int> &inbvec, bool msb_first)
{
  int i, temp = 0;
  int sizebvec = inbvec.size();
  if (msb_first) {
    for (i = 0; i < sizebvec; i++) {
      temp += pow(2,sizebvec - i - 1) * inbvec[i];
    }
  }
  else {
    for (i = 0; i < sizebvec; i++) {
      temp += pow(2,i) * inbvec[i];
    }
  }
  return temp;
}

template<typename kind>
inline std::vector<kind> mylib::getMid(const std::vector<kind> &input, const int start, const int number)
{
  assert((start >= 0) && (number >= 0));

  std::vector<kind> temp(number);

  for(int i = 0; i < number; i++){
    temp[i] = input[start+i];
  }

  return temp;
}

template<typename kind>
inline std::vector<kind> mylib::concat(const std::vector<kind> &forwardVec, const std::vector<kind> &backwardVec)
{
  std::vector<kind> temp = forwardVec;

  for(int i = 0; i < static_cast<int>(backwardVec.size()); i++){
    temp.push_back(backwardVec[i]);
  }

  return temp;
}

template<typename type>
inline type mylib::maxOf(type x1, type x2)
{
  if(x1 > x2){
    return x1;
  }
  else{
    return x2;
  }
}

template<typename type>
type mylib::minOf(type x1, type x2)
{
  if(x1 < x2){
    return x1;
  }
  else{
    return x2;
  }
}

inline int mylib::DoubleToShort(double input)
{

  int output;
  if(input>=0){
    output = mylib::round_uPoint(input*32767);
    while(output>32767)
      output -= 32767; 
  }
  else{
    output = mylib::round_uPoint(input*32768);
    while(output<-32768)
      output += 32768;
  }

  return output;
}

inline std::vector<int> mylib::DoubleToShort(const std::vector<double> &input)
{
  std::vector<int> output(input.size());
  for(int i = 0; i < static_cast<int>(input.size()); i++){
    if(input[i]>=0){
      output[i] = mylib::round_uPoint(input[i]*32767);
      while(output[i]>32767)
	output[i] -= 32767; 
    }
    else{
      output[i] = mylib::round_uPoint(input[i]*32768);
      while(output[i]<-32768)
	output[i] += 32768;
    }
  }

  return output;
}

inline int mylib::round_uPoint(double a)
{
  double out;

  out = floor(a + 0.5);
  
  return (int)out;
}
  
inline int mylib::sgn(double in)
{
  if(in >= 0.0)
    return 1;
  else
    return -1;
}


#endif
