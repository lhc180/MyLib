#ifndef MYPROGJPEG_H
#define MYPROGJPEG_H

#include <vector>
#include "mymatrix.h"
#include "myjpeg.h"
/************************************************************************************
 * myprogjpeg.h
 *   
 * Progressive JPEGに対応するための色んなクラスを入れる予定。
 * フレームの並び替えに関するクラスとかエントロピー符号とか。
 *
 * Contents:
 *   FunctionName of ClassName
 *
 * Last Updated: <2013/07/13 18:58:19 from yoshitos-mac-mini.local by yoshito>
 ************************************************************************************/

namespace mylib
{
  // ジグザグシーケンスデータのvectorを低周波から順に並べていく。
  // すべての要素サイズが64であると仮定しているのでそれ以外だと上手くいかない
  template < typename kind > 
  inline itpp::Vec< kind > ForwardProgJpegArrange(const std::vector< itpp::Vec< kind > >& input)
  {
    int numFrame = input.size();
    int totalSize = DCTSIZE2 * numFrame;
    itpp::Vec< kind > output(totalSize);
    
    int index = 0;
    for (int elem = 0; elem < DCTSIZE2; ++elem){
      for (int frame = 0; frame < numFrame; ++frame){
        output[index] = input[frame][elem];
        ++index;
      } // for frame
    } // for elem

    return output;
  }

  template < typename kind > 
  inline itpp::Vec< kind > ForwardProgJpegArrange(const itpp::Vec< itpp::Vec< kind > >& input)
  {
    int numFrame = input.size();
    int totalSize = DCTSIZE2 * numFrame;
    itpp::Vec< kind > output(totalSize);
    
    int index = 0;
    for (int elem = 0; elem < DCTSIZE2; ++elem){
      for (int frame = 0; frame < numFrame; ++frame){
        output[index] = input[frame][elem];
        ++index;
      } // for frame
    } // for elem

    return output;
  }

  
  template < typename kind >
  inline std::vector< itpp::Vec< kind > > InverseProgJpegArrange(const itpp::Vec< kind >& input)
  {
    int totalSize = input.size();
    assert(totalSize % DCTSIZE2 == 0); // 64で割れなければいけない
    int numFrame = totalSize / DCTSIZE2;
    std::vector< itpp::Vec< kind > > output(numFrame, itpp::Vec< kind >(DCTSIZE2));
    
    for (int index = 0; index < totalSize; ++index){
      output[index % numFrame][index / numFrame] = input[index];
    } // for index
    
    return output;
  }
  
  template < typename kind >
  inline std::vector< itpp::Vec< kind > > InverseProgJpegArrange(const itpp::Vec< kind >& dcCoef,
                                                                 const itpp::Vec< kind >& acCoef)
  {
    itpp::Vec< kind > input = itpp::concat(dcCoef, acCoef);
    return InverseProgJpegArrange(input);
  }

  /************************************************************************************
   * ForwardDPCM -- DC成分だけDPCMにする
   * 
   * Arguments:
   *   input -- ジグザグシーケンスのブロックが複数入ったvector
   *
   * Return Value:
   *   output -- ジグザグシーケンスのブロックが複数入ったvector
   ************************************************************************************/
  template < typename kind >
  inline std::vector< itpp::Vec< kind > > ForwardDPCM(const std::vector< itpp::Vec< kind > >& input)
  {
    std::vector< itpp::Vec < kind > > output = input;

    kind preCoef = 0;
    for (int i = 0; i < input.size(); ++i){
      output[i][0] -= preCoef;
      preCoef = input[i][0];
    } // for i
    return output;
  }

  template < typename kind >
  inline itpp::Vec< itpp::Vec< kind > > ForwardDPCM(const itpp::Vec< itpp::Vec< kind > >& input)
  {
    itpp::Vec< itpp::Vec < kind > > output = input;

    kind preCoef = 0;
    for (int i = 0; i < input.size(); ++i){
      output[i][0] -= preCoef;
      preCoef = input[i][0];
    } // for i
    return output;
  }

  
  template < typename kind >
  inline std::vector< itpp::Vec< kind > > InverseDPCM(const std::vector< itpp::Vec< kind > >& input)
  {
    std::vector< itpp::Vec< kind > > output = input;

    kind preCoef = 0;
    for (int i = 0; i < input.size(); ++i){
      preCoef += output[i][0];
      output[i][0] = preCoef;
    } // for i
    return output;
  }

  template < typename kind >
  inline itpp::Vec< itpp::Vec< kind > > InverseDPCM(const itpp::Vec< itpp::Vec< kind > >& input)
  {
    itpp::Vec< itpp::Vec< kind > > output = input;

    kind preCoef = 0;
    for (int i = 0; i < input.size(); ++i){
      preCoef += output[i][0];
      output[i][0] = preCoef;
    } // for i
    return output;
  }

  
}


#endif
