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
 * Last Updated: <2013/07/04 15:07:15 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/

namespace mylib
{
  // ジグザグシーケンスデータのvectorを低周波から順に並べていく。
  // すべての要素サイズが64であると仮定しているのでそれ以外だと上手くいかない
  template < typename kind > 
  inline std::vector< kind > ForwardProgJpegArrange(const std::vector< std::vector< kind > >& input)
  {
    u_int numFrame = input.size();
    u_int totalSize = DCTSIZE2 * numFrame;
    std::vector< kind > output(totalSize);
    
    u_int index = 0;
    for (int elem = 0; elem < DCTSIZE2; ++elem){
      for (u_int frame = 0; frame < numFrame; ++frame){
        output[index] = input[frame][elem];
        ++index;
      } // for frame
    } // for elem

    return output;
  }

  template < typename kind >
  inline std::vector< std::vector< kind > > InverseProgJpegArrange(const std::vector< kind >& input)
  {
    u_int totalSize = input.size();
    assert(totalSize % DCTSIZE2 == 0); // 64で割れなければいけない
    u_int numFrame = totalSize / DCTSIZE2;
    std::vector< std::vector< kind > > output(numFrame, std::vector< kind >(DCTSIZE2, 0));
    
    for (u_int index = 0; index < totalSize; ++index){
      output[index % numFrame][index / numFrame] = input[index];
    } // for index
    
    return output;
  }

  template < typename kind >
  inline std::vector< std::vector< kind > > InverseProgJpegArrange(const std::vector< kind >& dcCoef, const std::vector< kind >& acCoef)
  {
    std::vector< kind > input = mylib::Concat(dcCoef, acCoef);
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
  inline std::vector< std::vector< kind > > ForwardDPCM(const std::vector< std::vector< kind > >& input)
  {
    std::vector< std::vector < kind > > output = input;

    kind preCoef = 0;
    for (int i = 0; i < input.size(); ++i){
      output[i][0] -= preCoef;
      preCoef = input[i][0];
    } // for i
    return output;
  }

  template < typename kind >
  inline std::vector< std::vector< kind > > InverseDPCM(const std::vector< std::vector< kind > >& input)
  {
    std::vector< std::vector< kind > > output = input;

    kind preCoef = 0;
    for (int i = 0; i < input.size(); ++i){
      preCoef += output[i][0];
      output[i][0] = preCoef;
    } // for i
    return output;
  }

  
}


#endif
