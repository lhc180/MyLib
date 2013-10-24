#ifndef LENGTH_ADJUSTER_H
#define LENGTH_ADJUSTER_H

#include <vector>
#include <itpp/itbase.h>

namespace mylib{
  // 最小公倍数(least common multiple)
  inline int lcm(int a, int b)
  {
    int l = a*b / itpp::gcd(a,b);
    return l;
  }


  // MLCの各レベルの情報ビットの長さを調整するクラス
  // 基底クラス
  class SourceLengthAdjuster_MLC
  {
  protected:
    std::vector< int > originalLength_;
    std::vector< int > encodingLength_;
  
  public:
    SourceLengthAdjuster_MLC(): originalLength_(0), encodingLength_(0)
    { }
    SourceLengthAdjuster_MLC(const std::vector< int > &originalLength,
                             const std::vector< int > &encodingLength):
      originalLength_(originalLength), encodingLength_(encodingLength)
    { }

    virtual ~SourceLengthAdjuster_MLC()
    {
      originalLength_.clear();
      encodingLength_.clear();
    }

    virtual void Set(const std::vector< int > &originalLength, const std::vector< int > &encodingLength)
    {
      originalLength_ = originalLength;
      encodingLength_ = encodingLength;
    }
  
    virtual std::vector< itpp::bvec > Forward(const std::vector< itpp::bvec > &input) const = 0;
    virtual std::vector< itpp::bvec > Inverse(const std::vector< itpp::bvec > &input) const = 0;
  };

  // 低いレベルの最後に高いレベルの最初の方をくっつける。
  class AdvanceAdjuster_MLC: public SourceLengthAdjuster_MLC
  {
  public:
    // Constructor
    AdvanceAdjuster_MLC()
    { }
  
    AdvanceAdjuster_MLC(const std::vector< int > &originalLength,
                        const std::vector< int > &encodingLength):
      SourceLengthAdjuster_MLC(originalLength, encodingLength)
    { }
    // Destructor
    virtual ~AdvanceAdjuster_MLC()
    { }

    virtual std::vector< itpp::bvec > Forward(const std::vector< itpp::bvec > &input) const;
    virtual std::vector< itpp::bvec > Inverse(const std::vector< itpp::bvec > &input) const;
  
  };

  // 固定数のビットごとを取り出して、低いレベルの最後に付加していく
  // 量子化ビットのMSBだけ取り出して付加するのに使える
  class SkipAdjuster_MLC: virtual public SourceLengthAdjuster_MLC
  {
  protected:
    int interval_;                    // スキップ番号
    
  public:
    // Constructor
    SkipAdjuster_MLC()
    { }
    SkipAdjuster_MLC(const std::vector< int > &originalLength,
                     const std::vector< int > &encodingLength,
                     int interval):
      SourceLengthAdjuster_MLC(originalLength, encodingLength), interval_(interval)
    { }
    // Destructor
    virtual ~SkipAdjuster_MLC()
    { }

    void Set(const std::vector< int > &originalLength, const std::vector< int > &encodingLength, int interval)
    {
      SourceLengthAdjuster_MLC::Set(originalLength, encodingLength);
      interval_ = interval;
    }
  
    virtual std::vector< itpp::bvec > Forward(const std::vector< itpp::bvec > &input) const;
    virtual std::vector< itpp::bvec > Inverse(const std::vector< itpp::bvec > &input) const;
  };

  // ランダムに取り出して低いレベルの最後に付加していく
  class RandomAdjuster_MLC: virtual public SourceLengthAdjuster_MLC
  {
  public:
    // Constructor
    RandomAdjuster_MLC()
    { }
    RandomAdjuster_MLC(const std::vector< int > &originalLength,
                       const std::vector< int > &encodingLength):
      SourceLengthAdjuster_MLC(originalLength, encodingLength)
    { }
    // Destructor
    virtual ~RandomAdjuster_MLC()
    { }

    virtual std::vector< itpp::bvec > Forward(const std::vector< itpp::bvec > &input) const;
    virtual std::vector< itpp::bvec > Inverse(const std::vector< itpp::bvec > &input) const;
  };

  // スキップしたものをランダムに取り出す
  class SkipRandomAdjuster_MLC: public  SkipAdjuster_MLC, RandomAdjuster_MLC
  {
  public:
    SkipRandomAdjuster_MLC(const std::vector< int > &originalLength,
                           const std::vector< int > &encodingLength,
                           int skip):
      SkipAdjuster_MLC(originalLength, encodingLength, skip)
    { }
    virtual ~SkipRandomAdjuster_MLC()
    { }

    virtual std::vector< itpp::bvec > Forward(const std::vector< itpp::bvec > &input) const;
    virtual std::vector< itpp::bvec > Inverse(const std::vector< itpp::bvec > &input) const;
  };

  // 各ビット位置それぞれで情報ビットを格納しておいたものを
  // 上のビットから順に埋めていく。
  // これまでと少し使い勝手が違う。
  // つまりinput[0]がMSBが格納されたitpp::bvecになる
  class NaturalAdjuster_MLC: public SourceLengthAdjuster_MLC
  {
  public:
    NaturalAdjuster_MLC()
    { }
    NaturalAdjuster_MLC(const std::vector< int > &originalLength,
                        const std::vector< int > &encodingLength):
      SourceLengthAdjuster_MLC(originalLength, encodingLength)
    { }
    virtual ~NaturalAdjuster_MLC();

    virtual std::vector< itpp::bvec > Forward(const std::vector< itpp::bvec > &input) const;
    virtual std::vector< itpp::bvec > Inverse(const std::vector< itpp::bvec > &input) const;
  
  };

  /************************************************************************************
   * AdjustInfoLength_forMLC -- 1次元の情報ビットをMLC用に各レベルの
   * LDPCの情報長に合うように先頭から割り当てていく
   * 
   * Arguments:
   *   infoBits -- 1次元の情報ビット
   *   infoLengths_LDPC -- 各レベルのLDPCの情報長
   *
   * Return Value:
   *   std::vector< itpp::bvec >  -- 各レベルに割り当てられた情報ビット
   ************************************************************************************/
  std::vector< itpp::bvec > AdjustInfoLength_forMLC(const itpp::bvec &infoBits, const std::vector< int > &infoLengths_LDPC, bool randomPadding = false);
  
  /************************************************************************************
   * AdjustInfoLength_fromMSD -- AdjustInfoLength_forMLCによって分割された情報ビットを1次元に並び替える
   * 
   * Arguments:
   *   decodedInfo -- MLCによって復号された多次元の情報ビット
   *   totalBits -- 一次元情報ビットの総情報長
   *
   * Return Value:
   *   itpp::bvec -- 1次元情報ビット
   ************************************************************************************/
  itpp::bvec AdjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo, unsigned int totalBits);
  
  
  // マルチメディアの1フレーム分のトータルの情報長を行重みと列重みからいい感じに分ける
  std::vector< int > AdaptSourceLength_MLC(int totalSourceBits,
                                           const std::vector< int > &rowWeights,
                                           const std::vector< int> &colWeights);
  itpp::ivec AdaptSourceLength_MLC(int totalSourceBits,
                                   const itpp::ivec &rowWeights,
                                   const itpp::ivec &colWeights);

}


#endif







