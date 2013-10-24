#include <cassert>
#include "../include/myutl.h"
#include "../include/length_adjuster.h"

namespace mylib{
  // class AdvanceAdjuster_MLC +++++++++++++++++++++++++++++++++++++
  inline std::vector< itpp::bvec > AdvanceAdjuster_MLC::Forward(const std::vector< itpp::bvec > &input) const
  {
    assert(input.size() == originalLength_.size()); // レベルのチェック
    int levels = input.size();
    
    std::vector< itpp::bvec > output(levels, itpp::bvec(0));
    std::vector< int > difference(levels); // encodingLength_との差分
    // Level 0の差分
    difference[0] = encodingLength_[0] - input[0].size();
    assert(difference[0] >= 0);
    // Level 1以降の差分
    for (int l = 1; l < levels; ++l){
      //    output[l].set_size(encodingLength_[l]);
      difference[l] = encodingLength_[l] - input[l].size() + difference[l-1];
      assert(difference[l] >= 0);
    } // for l

    // Level 0から順にoutputを作っていく
    int startPos = 0;
    for (int l = 0; l < levels-1; ++l){
      output[l] = itpp::concat(input[l].get(startPos, -1), input[l+1].mid(0, difference[l]));
      startPos = difference[l];
    } // for l
    // 最後のレベルだけ別枠
    output[levels - 1] = itpp::concat(input[levels - 1].get(startPos, -1), itpp::bvec(difference[levels - 1])); // ## bvecの初期値がわからないが
    // 最後は0で埋まっているはず

    return output;
  }

  inline std::vector< itpp::bvec > AdvanceAdjuster_MLC::Inverse(const std::vector< itpp::bvec > &input) const
  {
    assert(input.size() == encodingLength_.size()); // レベルのチェック
    int levels = input.size();

    std::vector< itpp::bvec > output(levels, itpp::bvec(0));
    std::vector< int > includedInLevel(levels); // それぞれのレベルに正しいレベルの信号が何ビット入っているか
    // Level 0の差分
    includedInLevel[0] = originalLength_[0];
    int difference = input[0].size() - includedInLevel[0];
    assert(difference >= 0);
    // Level 1以降の差分
    for (int l = 1; l < levels; ++levels){
      includedInLevel[l] = originalLength_[l] - difference;
      difference = input[l].size() - includedInLevel[l];
      assert(difference >= 0);
    } // for l
  
    // Level 0
    output[0] = input[0].mid(0, originalLength_[0]);
    int endPos = includedInLevel[0];
    // Level 1以降
    for (int l = 1; l < levels; ++l){
      output[l] = itpp::concat(input[l-1].get(endPos, -1), input[l].mid(0, includedInLevel[l]));
      endPos = includedInLevel[l];
    } // for l

    return output;
  }
  // end of AdvanceAdjuster_MLC -----------------------------------------

  // class NaturalAdjuster_MLC +++++++++++++++++++++++++++++++++++++++++
  inline std::vector< itpp::bvec > NaturalAdjuster_MLC::Forward(const std::vector<itpp::bvec> &input) const
  {
    assert(input.size() == originalLength_.size());
    int rank             = originalLength_.size();
    int levels           = encodingLength_.size();

    std::vector< itpp::bvec > output(levels);

    // outputの長さ設定
    for (int l = 0; l < levels; ++l){
      output[l].set_size(encodingLength_[l]);
    }                             // for l
  
    int l                      = 0;
    int k                      = 0;
    for (int r = 0; r < rank; ++r){
      assert( input[r].size() == originalLength_[r]);
      for (int i = 0; i < originalLength_[r]; ++i){
        if (k >= encodingLength_[l]){
          l++;
          k                    = 0;
        }
        output[l][k] = input[r][i];
        k++;
      }                           // for i
    }                             // for r

    if (l != levels-1){
      std::cerr << "Error: Information bits are too few!! in NaturalAdjuster_MLC.Forward()." << std::endl;
      exit(1);
    }
  
    while (k < encodingLength_[levels - 1]){
      output[levels-1][k] = itpp::bin(0);
      k++;
    }                             // while k

    return output;
  }

  inline std::vector< itpp::bvec > NaturalAdjuster_MLC::Inverse(const std::vector<itpp::bvec> &input) const
  {
    assert(input.size() == encodingLength_.size());
    int rank   = originalLength_.size();
    int levels = encodingLength_.size();

    // 長さチェック
    for (int l = 0; l < levels; ++l){
      assert(input[l].size() == encodingLength_[l]);
    } // for l
    
    std::vector< itpp::bvec > output(rank);

    // outputの長さの設定
    for (int r = 0; r < rank; ++r){
      output[r].set_size(originalLength_[r]);
    } // for r

    int l = 0;
    int k = 0;
    for (int r = 0; r < rank; ++r){
      for (int i = 0; i < originalLength_[r]; ++i){
        if (k >= encodingLength_[l]){
          l++;
          k = 0;
        }
        output[r][i] = input[l][k];
        k++;
      } // for i
    } // for r

    return output;
  }

  // 最後のレベルの最後は0で埋める
  std::vector< itpp::bvec > AdjustInfoLength_forMLC(const itpp::bvec &infoBits, const std::vector< int > &infoLengths_LDPC, bool randomPadding)
  {
    int levels = infoLengths_LDPC.size();
    u_int totalInfoBits = mylib::Sum(infoLengths_LDPC);
    if (static_cast< u_int >(infoBits.size()) > totalInfoBits){
      std::cerr << "Error: Information bits is longer than that of LDPC codes in AdjustInfoLength_forMLC." << std::endl;
      exit(1);
    } 

    std::vector< itpp::bvec > output(levels);
    for (int i = 0; i < levels; ++i){
      if (randomPadding){
        output[i] = itpp::randb(infoLengths_LDPC[i]);
      } // if 
      else{
        output[i].set_size(infoLengths_LDPC[i]);
        output[i].zeros();
      } // else 
    } // for i
    
    int startPos = 0;
    for (int i = 0; i < levels; ++i){
      if (startPos + infoLengths_LDPC[i] >= infoBits.size()){
        output[i].set_subvector(0, infoBits.get(startPos, -1));
        std::cout << std::endl;
        break;
      } // if
      else {
        output[i].set_subvector(0, infoBits.mid(startPos, infoLengths_LDPC[i]));
        startPos += infoLengths_LDPC[i];
      }
    } // for i

    return output;
  }

  itpp::bvec AdjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo, unsigned int totalBits)
  {
    int levels = decodedInfo.size();

    u_int totalDecodedBits = 0;
    for (int i = 0; i < levels; ++i){
      totalDecodedBits += decodedInfo[i].size();
    } // for i
    if (totalDecodedBits < totalBits){
      std::cerr << "Error: Decoded bits are shorter than an original total bits in AdjustInfoLength_fromMSD." << std::endl;
      exit(1);
    }

    itpp::bvec concatenated(0);
    for (int i = 0; i < levels; ++i){
      concatenated = itpp::concat(concatenated, decodedInfo[i]);
    } // for i
    
    itpp::bvec output = concatenated.mid(0, totalBits);
    
    return output;
  }

  // ## 今のところ他のクラスは保留


  // // マルチメディアの1フレーム分のトータルの情報長を行重みと列重みからいい感じに分ける
  // inline std::vector< int > AdaptSourceLength_MLC(int totalSourceBits,
  //                                                 const std::vector< int > &rowWeights,
  //                                                 const std::vector< int > &colWeights)
  // {
  
  
  // }

}
