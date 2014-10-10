#ifndef MYINTERLEAVERS_H
#define MYINTERLEAVERS_H

#include <itpp/itbase.h>

namespace mylib {
  /************************************************************************************
   * RandomInterleaver -- ランダムインタリーバを生成する
   * 
   * Arguments:
   *   length -- インタリーバ長
   *
   * Return Value:
   *   itpp::ivec -- ランダムなインデクスの並び
   ************************************************************************************/
  inline itpp::ivec RandomInterleaver(int length)
  {
    itpp::ivec interleaver = itpp::sort_index(itpp::randu(length));
    
    return interleaver;
  }
  
  static bool S_RandomInterleaver(itpp::ivec* output, int size, int spread, int maxTrial = 500)
  {
    bool success = false;
    itpp::ivec interleaver = RandomInterleaver(size);

    for (int reverse = 0; reverse < maxTrial; ++reverse){
      std::cout << "reverse = " << reverse << '\r' << std::flush;
      int pivot;
      for (pivot = 0; pivot < size; ++pivot){
        int start_i = pivot - spread;
        if (start_i < 0){
          start_i = 0;
        } // if
        int swap_i;
        for (swap_i = pivot; swap_i < size; ++swap_i){
          int criteria = interleaver[swap_i];
          int i;
          for (i = start_i; i < pivot; ++i){
            if (std::abs(interleaver[i] - criteria) < spread){ // スプレッド値より下回ったら次のswap候補を探す
              // ここの判定に=付けるかどうかで若干特性変わる
              break;
            } // if 
          } // for i
          if (i == pivot){      // もしpivotの値まで比較が成功していたら
            std::swap(interleaver[pivot], interleaver[swap_i]);
            break;
          } // if i
        } // for swap_i
        if (swap_i == size){      // もしpivot以降にswapできる要素が無ければ
          break;
        } // if 
      } // for pivot
      if (pivot == size){
        success = true;
        break;
      } // if
      else{
        interleaver = itpp::reverse(interleaver);
      } // else
    } // for reverse
    std::cout << std::endl;   // ##
  
    *output = interleaver;
    return success;
  }

  
  template < typename kind >
  itpp::Vec< kind > Interleave(const itpp::Vec< kind >& input, const itpp::ivec& interleaver)
  {
    assert(input.size() == interleaver.size());
    itpp::Vec< kind > output(input.size());

    for (int i = 0; i < input.size(); ++i){
      output[i] = input[interleaver[i]];
    } // for i

    return output;
  }

  template < typename kind >
  itpp::Vec< kind >Deinterleave(const itpp::Vec< kind >& input, const itpp::ivec& interleaver)
  {
    assert(input.size() == interleaver.size());
    itpp::Vec< kind > output(input.size());

    for (int i = 0; i < input.size(); ++i){
      output[interleaver[i]] = input[i];
    } // for i

    return output;
  }
}

#endif
