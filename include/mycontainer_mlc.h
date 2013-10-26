#ifndef MYCONTAINER_MLC_H
#define MYCONTAINER_MLC_H

#include <myutl.h>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
/************************************************************************************
 * mycontainer_mlc
 *   
 * MLC用にオリジナルのコンテナを作る
 *
 * Contents:
 *   FunctionName of ClassName
 *
 * Last Updated: <2013/10/26 21:37:59 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/


namespace mylib{

  /************************************************************************************
   * BvecStack
   * 
   * bvecにsubvectorをセットしていって、容量を超える入力が入ったらfalseを返す
   ************************************************************************************/
  class BvecStack
  {
  private:
    itpp::bvec container_;
    std::vector< int > sizesOfElements_;
    int lastIndex_;
        
  public:
    BvecStack(): container_(0), lastIndex_(0)
    { }
    
    explicit BvecStack(int num, int init = 0)
    {
      Init(num, init);
    }
    
    virtual ~BvecStack()
    { }

    // padding bits of init=0,1,-1.
    // -1 provides random bits.
    void Init(int num, int init = -1)
    {
      lastIndex_ = 0;
      sizesOfElements_.resize(0);
      
      switch (init){
      case 0:
        container_.set_size(num);
        container_.zeros();
        break;

      case 1:
        container_.set_size(num);
        container_.ones();
        break;
      
      case -1:
        container_ = itpp::randb(num);
        break;
      
      default:
        std::cerr << "Error: the value of \"init\" is unavailable in BvecStack." << std::endl;
        exit(1);
        break;
      }
    }
    
    
    bool Add(const itpp::bvec& input)
    {
      if (lastIndex_ + input.size() > container_.size()){
        return false;
      } // if
      else{
        container_.set_subvector(lastIndex_, input);
        sizesOfElements_.push_back(input.size());
        lastIndex_ += input.size();
        return true;
      } // else 
    }

    itpp::bvec Get()
    {
      return container_;
    }

    std::vector< int > Sizes()
    {
      return sizesOfElements_;
    }
    
  };

  
  /************************************************************************************
   * BvecStack_MLC 
   * 
   * MLC用にBvecStackをレベル分だけ用意して詰めていく
   ************************************************************************************/
  class BvecStack_MLC
  {
  private:
    std::vector< BvecStack > stacks_;
    std::vector< double > rates_;
    int totalSize_;
    int storedSize_;
    
  public:
    BvecStack_MLC(const std::vector< int >& sizes, int init = 0)
    {
      Init(sizes, init);
    }
      
    virtual ~BvecStack_MLC()
    { }

    void Init(const std::vector< int >& sizes, int init = 0)
    {
      stacks_.resize(sizes.size());
      rates_.resize(sizes.size());
      
      totalSize_ = mylib::Sum(sizes);
      storedSize_ = 0;
      
      for (u_int i = 0; i < sizes.size(); ++i){
        stacks_[i].Init(sizes[i], init);
        rates_[i] = static_cast< double >(sizes[i])/ static_cast< double >(totalSize_);
      } // for i
    }
    
    bool Add(const itpp::bvec& input)
    {
      if (storedSize_ + input.size() > totalSize_){
        return false;
      } // if totalSize_
      
      std::vector< itpp::bvec > subvecs(rates_.size());
      double sumRates = 0;
      int previousIndex = 0;
      for (u_int level = 0; level < rates_.size(); ++level){
        sumRates += rates_[level];
        int lastIndex = std::floor(input.size() * sumRates);
        if (level == rates_.size()-1){
          lastIndex = 0;
        } // if

        subvecs[level] = input.get(previousIndex, lastIndex-1); // 実際に入る数を考慮して-1してある
        
        previousIndex = lastIndex;
      } // for i

      for (u_int level = 0; level < rates_.size(); ++level){
        if (!(stacks_[level].Add(subvecs[level]))){
          return false;
        } // if i
      } // for i

      storedSize_ += input.size();
      
      return true;
    }

    std::vector< itpp::bvec > Get()
    {
      std::vector< itpp::bvec > output(rates_.size());
      
      for (u_int i = 0; i < rates_.size(); ++i){
        output[i] = stacks_[i].Get();
      } // for

      return output;
    }

    std::vector< std::vector< int > > Sizes()
    {
      std::vector< std::vector< int > > sizes(rates_.size(), std::vector< int >(0));

      for (u_int i = 0; i < rates_.size(); ++i){
        sizes[i] = stacks_[i].Sizes();
      } // for i
      
      return sizes;
    }
  };

  /************************************************************************************
   * DivideBvec -- inputをsizesごとに分割して返す

   * 
   * Arguments:
   *   input -- 入力bvec
   *   sizes -- 各subvectorのサイズ
   *
   * Return Value:
   *   output -- 分割されたbvecをstd::vectorに入れて返す
   ************************************************************************************/
  std::vector< itpp::bvec > DivideBvec(const itpp::bvec& input, const std::vector< int >& sizes)
  {
    int sum = mylib::Sum(sizes);
    assert(sum <= input.size());

    int elements = sizes.size();
    std::vector< itpp::bvec > output(elements, itpp::bvec(0));
    int lastIndex = 0;
    for (int i = 0; i < elements; ++i){
      output[i] = input.mid(lastIndex, sizes[i]);
      lastIndex += sizes[i];
    } // for i

    return output;
  }


  /************************************************************************************
   * DivideBvec_MSD -- MLCのフレームの前方から情報ビットを取り出していって合成したものを返す
   * 
   * Arguments:
   *   input -- MLCのフレーム
   *   sizes -- 各情報ビットが各レベルに何ビット含まれているか
   *
   * Return Value:
   *   output -- 複数の情報ビット
   ************************************************************************************/
  std::vector< itpp::bvec > DivideBvec_MSD(const std::vector< itpp::bvec >& input,
                                           const std::vector< std::vector< int > >& sizes,
                                           int effectiveLevel = -1)
  {
    assert(input.size() == sizes.size());
    assert(effectiveLevel <= static_cast< int >(input.size()));
    if (effectiveLevel < 0){
      effectiveLevel = input.size();
    } // if 
    
    std::vector< itpp::bvec > output(sizes[0].size(), itpp::bvec(0));
    for (int level = 0; level < effectiveLevel; ++level){
      std::vector< itpp::bvec > t_seqs = DivideBvec(input[level], sizes[level]);
      for (u_int i = 0 ; i < t_seqs.size(); ++i){
        output[i] = itpp::concat(output[i], t_seqs[i]);
      } // for i
    } // for level

    return output;
  }
  
}

#endif










