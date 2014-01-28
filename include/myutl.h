#ifndef MYUTL_H
#define MYUTL_H

#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<cassert>
#include<cmath>
#include<itpp/itbase.h>
// #include"./mymatrix.h"
// #include"./mymatrix_utl.h"

/***************************/
// itppが使えない環境用の各種変換用関数
// また、その他の変換用の関数など
//
/****************************/

typedef unsigned char u_char;
typedef unsigned int  u_int;

namespace mylib
{
  const double INFTY = 1e30;

  // itppとstdの相互変換
  template < typename kind >
  itpp::Vec< kind > ToItppVec(const std::vector< kind >& input);

  template < typename kind >
  std::vector< kind > ToStdVector(const itpp::Vec< kind >& input);
  
  // int型のvectorをnormで正規化する
  inline void                normalizeData(const std::vector<int> &input, std::vector<double> &output, unsigned norm);
  inline std::vector<double> normalizeData(const std::vector<int> &input, unsigned norm);

  // double型のvectorをnormでdenormalizeする
  inline void             denormalizeData(const std::vector<double> &input, std::vector<int> &output, unsigned norm);
  inline std::vector<int> denormalizeData(const std::vector<double> &input, unsigned norm);

  // fileにvectorを出力するが、withIndexがtrueならインデックスも付加する
  template<typename kind>
  void fprintVec(std::ofstream            &file,
                 const std::vector<kind> &input,
                 bool                     withIndex = true);

  // start~endをdivで刻んだ値にする
  template<typename kind>
  std::vector<kind> setVec(kind start, kind end, kind div);

  // vectorを0で初期化する
  template<typename type>
  void vecZeros(std::vector<type> &input);

  // calculating average of vector
  template<typename type>
  double vecAverage(std::vector<type> &input);
  
  // 10進数から2進数に変換
  // MLBから出力
  std::vector<int> dec2bin(int length, int index);

  // 2進数から10進数に変換
  int bin2dec(const std::vector<int> &inbvec, bool msb_first = false);

  // inputのstart番目の要素からnumber分抜き出す
  // numberに負の数を入れた場合は、最後の要素まで抜き出す
  template<typename kind>
  std::vector<kind> Mid(const std::vector<kind> &input, const int start, const int number);

  // forwardVecにbackwardVecを結合する
  template<typename kind>
  std::vector<kind> Concat(const std::vector<kind> &forwardVec, const std::vector<kind> &backwardVec);
  template <typename kind>
  std::vector<kind> Concat(kind forward, const std::vector<kind>& backwardVec);
  
  
  template<typename type>
  type maxOf(type x1, type x2);

  template<typename type>
  type minOf(type x1, type x2);

  template <typename kind>
  kind Max(const std::vector< kind > &input);

  template <typename kind>
  kind Min(const std::vector< kind > &input);

  template <typename kind>
  kind Sum(const std::vector< kind > &input);
  
  //   int toShort(double input);
  //   std::vector<int> toShort(const std::vector<double> &input);
  
  template < typename kind >
  std::vector<int> ToInt(const std::vector<kind> &input);

  template< typename kind >
  double toDouble(kind input);

  template <typename kind>
  std::vector<double> toDouble(const std::vector<kind> &input);

  //小数点以下を四捨五入
  int round_uPoint(double input);

  // 符号
  int sgn(double in);

  // bSetがfalseなら標準入力からの入力をinputに入れる。
  // trueなら標準出力に出力する。
  void setString(std::string &input, bool &bSet);

  // inputFileNameの毎行ごとのデータをvectorで返す
  template<typename kind>
  std::vector< kind > fReadLines(const char* inputFileName, int lines);

  // ファイルの拡張子を取り出す
  std::string FileExtention(const char* fileName);

  // ファイルのサイズを返す
  size_t FileSize(const char* fileName);
  
  // 2のべき乗かどうかのチェック
  bool Radix2(int input);

  // 絶対値
  std::vector< double > Abs(const std::vector< double > &input);

  template < typename kind >
  inline std::vector< kind > LevelShift(const std::vector < kind > &input, kind level)
  {
    std::vector< kind > output(input.size());
    for(int i = 0; i < static_cast< int >(input.size()); i++)
      {
        output[i] = input[i] + level;
      }

    return output;
  
  }

  template<typename kind>
  inline std::vector< kind > Divide(const std::vector< kind > &input, const std::vector< kind > &qTable)
  {
    assert(input.size() == qTable.size());
  
    std::vector< kind > output(input.size());
    for(int i = 0; i < input.size(); i++)
      {
        output[i] = input[i] / qTable[i];
      }

    return output;
  }

  template<typename kind>
  inline std::vector< kind > Multiply(const std::vector< kind > &input, const std::vector< kind > &qTable)
  {
    assert(input.size() == qTable.size());
  
    std::vector< kind > output(input.size());
    for(int i = 0; i < input.size(); i++)
      {
        output[i] = input[i] * qTable[i];
      }

    return output;
  }

  template < typename kind >
  inline itpp::Vec< kind > ToItppVec(const std::vector< kind >& input)
  {
    itpp::Vec< kind > output(input.size());

    int i = 0;
    for (typename std::vector< kind >::const_iterator ite = input.begin(); ite != input.end(); ++ite){
      output[i] = *ite;
      ++i;
    } // for ite
    return output;
  }

  template < typename kind >
  inline std::vector< kind > ToStdVector(const itpp::Vec< kind >& input)
  {
    std::vector< kind > output(input.size());
    for (int i = 0, size = output.size(); i < size; ++i){
      output[i] = input[i];
    } // for i
    return output;
  }
  
} // end of namespace mylib

template<typename kind>
inline std::vector<kind> mylib::setVec(kind start, kind end, kind div)
{
  assert((start <= end));
  assert((div > 0));

  int size = (end - start)/div + 1;
  std::vector< kind > output(size);

  for(int i = 0; i < size; i++)
    {
      output[i] = start + i*div;
    }

  return output;
}

template<typename type>
inline void mylib::vecZeros(std::vector<type> &input)
{
  for(int i = 0; i < static_cast<int>(input.size()); i++){
    input[i] = 0;
  }
}

template<typename type>
inline double mylib::vecAverage(std::vector<type> &input)
{
  double ave  = 0.0;
  for(int i = 0; i < static_cast<int>(input.size()); i++){
    ave      += input[i];
  }
  ave /= static_cast<double>(input.size());

  return ave;
}

inline void mylib::normalizeData(const std::vector<int> &input, std::vector<double> &output, unsigned norm)
{
  output.resize(input.size());
  for(unsigned i = 0; i < input.size(); i++){
    output[i] = static_cast<double>(input[i])/static_cast<double>(norm);
  }
}

inline std::vector<double> mylib::normalizeData(const std::vector<int> &input, unsigned norm)
{
  std::vector<double> output(0);

  mylib::normalizeData(input, output, norm);

  return output;
}

inline void mylib::denormalizeData(const std::vector<double> &input, std::vector<int> &output, unsigned norm)
{
  output.resize(input.size());
  for(unsigned i = 0; i < input.size(); i++){
    output[i] = static_cast<int>(input[i] * norm);
  }
}

inline std::vector<int> mylib::denormalizeData(const std::vector<double> &input, unsigned norm)
{
  std::vector<int> output(0);

  mylib::denormalizeData(input, output, norm);

  return output;
}

template<typename kind> 
inline void mylib::fprintVec(std::ofstream           &file, 
		      const std::vector<kind> &input,
		      bool withIndex)
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

  // LSBから出力
inline std::vector<int> mylib::dec2bin(int length, int index)
{
  int              i, bintemp = index;
  std::vector<int> temp(length);
 
  for (i = 0; i < length; i++) {
    temp[i] = bintemp & 1;
    bintemp = (bintemp >> 1);
  }
  return temp;
}

inline int mylib::bin2dec(const std::vector<int> &inbvec, bool msb_first)
{
  int i, temp   = 0;
  int sizebvec  = inbvec.size();
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
inline std::vector<kind> mylib::Mid(const std::vector<kind> &input, const int start, const int number)
{
  if (number >= 0){
    return std::vector< kind >(input.begin()+start, input.begin()+ start+ number);
  } // if 
  else{
    return std::vector< kind >(input.begin()+start, input.end());
  } // else
}

template<typename kind>
inline std::vector<kind> mylib::Concat(const std::vector<kind> &forwardVec, const std::vector<kind> &backwardVec)
{
  std::vector<kind> temp = forwardVec;

  temp.insert(temp.end(), backwardVec.begin(), backwardVec.end());

  return temp;
}

template <typename kind>
inline std::vector<kind> mylib::Concat(kind forward, const std::vector<kind>& backwardVec)
{
  std::vector<kind> temp(1);
  temp[0] = forward;
  std::vector<kind> output = mylib::Concat(temp, backwardVec);

  return output;
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

template <typename kind>
inline kind mylib::Max(const std::vector< kind > &input)
{
  kind max = input[0];

  for (typename std::vector<kind>::const_iterator curInput = input.begin();
       curInput != input.end(); ++curInput){
    if (max < *curInput){
      max = *curInput;
    }
  } // for i

  return max;
}

template<typename type>
inline type mylib::minOf(type x1, type x2)
{
  if(x1 < x2){
    return x1;
  }
  else{
    return x2;
  }
}

template <typename kind>
inline kind mylib::Min(const std::vector< kind > &input)
{
  kind min = input[0];
  for (typename std::vector< kind >::const_iterator curInput = input.begin();
       curInput != input.end(); ++curInput){
    if (min > *curInput){
      min = *curInput;
    }
  } // for i

  return min;
}

template <typename kind>
inline kind mylib::Sum(const std::vector< kind > &input)
{
  kind sum = 0;
  for (typename std::vector< kind >::const_iterator curInput = input.begin();
       curInput != input.end(); ++curInput){
    sum += *curInput;
  } // for curInput

  return sum;
}

// inline int mylib::toShort(double input)
// {

//   int output;
//   if(input>=0){
//     output    = mylib::round_uPoint(input*32767);
//     while(output>32767)
//       output -= 32767; 
//   }
//   else{
//     output    = mylib::round_uPoint(input*32768);
//     while(output<-32768)
//       output += 32768;
//   }

//   return output;
// }

// inline std::vector<int> mylib::toShort(const std::vector<double> &input)
// {
//   std::vector<int> output(input.size());
//   for(int i = 0; i < static_cast<int>(input.size()); i++){
//     if(input[i]>=0){
//       output[i]  = mylib::round_uPoint(input[i]*32767);
//       while(output[i]>32767)
// 	output[i]  -= 32767; 
//     }
//     else{
//       output[i]  = mylib::round_uPoint(input[i]*32768);
//       while(output[i]<-32768)
// 	output[i]  += 32768;
//     }
//   }

//   return output;
// }


template < typename kind >
inline std::vector<int> mylib::ToInt(const std::vector<kind> &input)
{
  std::vector<int> output(input.size());

  for(int i = 0; i < static_cast<int>(input.size()); i++){
    output[i] = static_cast< int >(input[i]);
  }

  return output;
}

template< typename kind >
inline
double mylib::toDouble(kind input)
{
  return static_cast<double>(input);
}

template< typename kind >
inline
std::vector<double> mylib::toDouble(const std::vector<kind> &input)
{
  std::vector<double> output(input.size());

  for(int i = 0; i < static_cast<int>(input.size()); i++){
    output[i] = static_cast<double>(input[i]);
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

inline void mylib::setString(std::string &input, bool &bSet)
{
  if(!bSet){
    std::cin >> input;
  }
  else{
    std::cout << input << std::endl;
  }

  bSet = true;
}

template< typename kind > inline
std::vector< kind > mylib::fReadLines(const char* inputFileName, int lines)
{
  std::ifstream inputFile(inputFileName);
  std::vector< kind > output(lines);

  for(int i = 0; i < lines; i++)
    {
      assert(!inputFile.eof());
      inputFile >> output[i];
    }

  return output;
}

inline std::string mylib::FileExtention(const char* filename_c)
{
  std::string filename  = filename_c;
  int         stringLen = filename.length();

  for (int i = stringLen - 1; i >= 0; --i)
    {
      if (filename[i] == '.')
        {
          assert(i < stringLen-1);
          return filename.substr(i+1);
        }
    }                           // for i

  std::cout << "Error: File extention unrecognized." << std::endl;
  exit(1);
}

inline size_t mylib::FileSize(const char* filename)
{
  std::ifstream infile(filename);
  size_t fileSize = (size_t)infile.seekg(0, std::ios::end).tellg();
  infile.close();
  return fileSize;
}

inline bool mylib::Radix2(int input)
{
  assert(input > 0);

  if ((input & (input-1)) == 0){
    return true;
  }
  else{
    return false;
  }
}

inline std::vector< double > mylib::Abs(const std::vector< double > &input)
{
  std::vector< double > output(input.size());

  for (int i = 0; i < static_cast< int >(input.size()); ++i){
    output[i] = fabs(input[i]);
  } // for i

  return output;
}

#endif
