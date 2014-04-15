#ifndef TURBO_CODE_H
#define TURBO_CODE_H

/************************************************************************************
 * turbo_code.h
 *   
 * ターボ符号に関するクラス
 *
 * Contents:
 *   class Rsc
 *   class TurboCode
 *
 * Last Updated: <2014/04/12 15:17:46 from Okauchi.local by yoshito>
 ************************************************************************************/

#include <cassert>
#include <boost/rational.hpp>
#include <itpp/itcomm.h>

namespace mylib{
  struct encodeTable{
    itpp::bin output_;
    int nextState_;
  };

  
  class Rsc
  {
  private:
    const static boost::rational< int > codeRate_; // 符号化率は1/2で固定
    const int constraint_;
    const int memory_;
    const int stateNum_;
    const unsigned int feedforward_;    // Octal Form
    const unsigned int feedback_;       // Octal Form
    std::vector< std::vector< encodeTable > > encodeTable_;
    std::vector< std::vector < int > > revEncodeTable_;
    mutable itpp::vec lambda_;
    mutable int lastState_;
    mutable itpp::bvec tailbitTable_;
    
  protected:
    virtual double Jacobian(double x1, double x2) const;

    virtual void CalcLambda(const itpp::cvec& received, const itpp::vec &logLikelihood_in,
                            double n0, bool knowLastState) const;

    virtual void CalcLLR_out(const itpp::cvec &received, const itpp::vec &logLikelihood_in,
                             itpp::vec *logLikelihood_out, double n0) const;

    // TurboCode内で使う用
    static void CalcLLR_out(const itpp::cvec &received, const itpp::vec &lambda,
                            const itpp::vec &logLikelihood_in, itpp::vec *logLikelihood_out, double n0);

    static itpp::bvec HardDecision(const itpp::vec& lambda);
    
  public:
    // feedforwardとfeedbackは8進表示なので0..で格納しておく
    // feedforward:g_1, feedback:g_0
    // constraint=5, feedforward=021, feedback=037
    Rsc(int constraint = 3, unsigned int feedforward = 05, unsigned int feedback = 07);
    
    virtual ~Rsc()                     // 継承はしない前提
    { }
    
    // Terminate is not supported.
    void GenParity(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec GenParity(const itpp::bvec& input) const{
      itpp::bvec output; GenParity(input, &output); return output;
    }

    void Terminate(itpp::bvec *tailbits, itpp::bvec *parity) const;

    // Supporting Termination
    void GenParityWithTerm(const itpp::bvec& input, itpp::bvec* tailbits, itpp::bvec* output) const;
    
    void Encode(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec Encode(const itpp::bvec& input) const{
      itpp::bvec output; Encode(input, &output); return output;
    }

    void EncodeWithTerm(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec EncodeWithTerm(const itpp::bvec& input) const{
      itpp::bvec output;
      EncodeWithTerm(input, &output);
      return output;
    }

    // Log-Map 復号
    void Decode(const itpp::cvec& received, const itpp::vec& logLikelihood_in,
                itpp::vec* logLikelihood_out, double n0, bool knowLastState = true) const;
    
    // BPSK_cにしかまだ対応していない
    void HardDecision(itpp::bvec* outputBits) const;
    itpp::bvec HardDecision() const
    {
      itpp::bvec outputBits;
      HardDecision(&outputBits);
      return outputBits;
    }

    int Constraint() const
    { return constraint_; }
    
    static boost::rational< int > CodeRate(){
      return codeRate_;
    }

    friend class TurboCode;
  };
  
  class TurboCode
  {
  private:
    const itpp::ivec interleaver_;
    const Rsc rsc1_, rsc2_;
    static const boost::rational< int > codeRate_; // 符号化率は1/3で固定
    const bool termination_;

    virtual void doEncode(const itpp::bvec& input, itpp::bvec* output) const; // NVI
    virtual void doEncodeWithTerm(const itpp::bvec& input, itpp::bvec* output) const;
    
    virtual void doDecode(const itpp::cvec& receivedSignal, itpp::bvec* output,
                          double n0, int iteration) const;

    virtual void doDecodeWithTerm(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                  double n0, int iteration) const;
    
    // ++++ Zero Padding ++++
    virtual void doDecodeWithZeroPadding(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                         double n0, int numPads, int iteration) const;

    virtual void doDecodeWithZeroPadding_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                              double n0, int numPads, int iteration) const;

    virtual void doDecodeWithZeroPaddingHDinMAP1(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                                 double n0, int numPads, int iteration) const;

    virtual void doDecodeWithZeroPaddingHDinMAP1_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                                      double n0, int numPads, int iteration) const;

    
    // MAP1に入力される外部値のみ補正
    virtual void doDecodeWithZP1(const itpp::cvec& received, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;
    virtual void doDecodeWithZP1_term(const itpp::cvec& received, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;

    // MAP2に入力される外部値のみ補正
    virtual void doDecodeWithZP2(const itpp::cvec& received, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;
    virtual void doDecodeWithZP2_term(const itpp::cvec& received, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;

    
    // ++++ Cyclic Suffix ++++
    virtual void doDecodeWithCyclicSuffix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                          double n0, int numPads, int iteration) const;

    virtual void doDecodeWithCyclicSuffix_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                               double n0, int numPads, int iteration) const;

    // MAP1に入力される外部値のみ補正
    virtual void doDecodeWithCS1(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;
    virtual void doDecodeWithCS1_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;

    // MAP2に入力される外部値のみ補正
    virtual void doDecodeWithCS2(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;
    virtual void doDecodeWithCS2_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;
    
    // ++++ Cyclic Prefix ++++
    virtual void doDecodeWithCyclicPrefix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                          double n0, int numPads, int iteration) const;

    virtual void doDecodeWithCyclicPrefix_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                               double n0, int numPads, int iteration) const;

    // MAP1に入力される外部値のみ補正
    virtual void doDecodeWithCP1(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;
    virtual void doDecodeWithCP1_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;

    // MAP2に入力される外部値のみ補正
    virtual void doDecodeWithCP2(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;
    virtual void doDecodeWithCP2_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                 double n0, int numPads, int iteration) const;
    
    // ++++ Inversed Prefix ++++
    virtual void doDecodeWithInversedPrefix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                            double n0, int numPads, int iteration) const;

    virtual void doDecodeWithInversedPrefix_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                                 double n0, int numPads, int iteration) const;

    // ++++ Cyclic Infix ++++
    virtual void doDecodeWithCyclicInfix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                         double n0, int start, int numPads, int iteration) const;
    
    virtual void doDecodeWithCyclicInfix_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                              double n0, int start, int numPads, int iteration) const;
    
  protected:
    virtual void SeparateReceivedSignalForZeroPadding(const itpp::cvec& received,
                                                      itpp::cvec* in1, itpp::cvec* in2, int numPads) const;
    virtual void ModifySignalForZeroPadding(itpp::cvec* received, int numPads) const;
    virtual void ModifyLLRForZeroPadding(itpp::vec* llr, int numPads) const;
    
    virtual void ModifyLLRForCyclicSuffix(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForCyclicPrefix(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForInversedPrefix(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForCyclicInfix(itpp::vec* llr, int start, int numPads) const;

    virtual void SeparateReceivedSignal(const itpp::cvec& receivedSignal,
                                        itpp::cvec* in1, itpp::cvec* in2) const;
    
  public:
    explicit TurboCode(itpp::ivec interleaver, int constraint = 3, int feedforward = 05, int feedback = 07,
                       bool termination = false):
      interleaver_(interleaver), rsc1_(constraint, feedforward, feedback), rsc2_(constraint, feedforward, feedback),
      termination_(termination)
    { }
        
    virtual ~TurboCode()
    { }

    // ここから全てNVIパターンを使った実装
    void Encode(const itpp::bvec& input, itpp::bvec* output) const
    {
      if (termination_){
        doEncodeWithTerm(input, output);
      } // if termination_
      else {
        doEncode(input, output);
      }
      
    }
    itpp::bvec Encode(const itpp::bvec& input) const
    {
      itpp::bvec output;
      Encode(input, &output);
      return output;
    }

    // Log-Map Decode
    void Decode(const itpp::cvec& receivedSignal, itpp::bvec* output, double n0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithTerm(receivedSignal, output, n0, iteration);
      } // if termination_
      else {
        doDecode(receivedSignal, output, n0, iteration); 
      } 
    }
    itpp::bvec Decode(const itpp::cvec& receivedSignal, double n0, int iteration = 10) const
    {
      itpp::bvec output;
      Decode(receivedSignal, &output, n0, iteration);
      return output;
    }

    // ++++++++ Zero Padding ++++++++
    void DecodeWithZeroPadding(const itpp::cvec& receivedSignal, itpp::bvec* output,
                               double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithZeroPadding_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithZeroPadding(receivedSignal, output, n0, numPads, iteration);
      } // else 
    }
    itpp::bvec DecodeWithZeroPadding(const itpp::cvec& receivedSignal,
                                     double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZeroPadding(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    // Hard DecisionをMAP Decoder 1で行う
    void DecodeWithZeroPaddingHDinMAP1(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                       double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithZeroPaddingHDinMAP1_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithZeroPaddingHDinMAP1(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithZeroPaddingHDinMAP1(const itpp::cvec& receivedSignal,
                                             double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZeroPaddingHDinMAP1(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }
    
    // MAP1に入る外部値だけ補正する
    void DecodeWithZP1(const itpp::cvec& receivedSignal, itpp::bvec* output,
                       double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithZP1_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithZP1(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithZP1(const itpp::cvec& receivedSignal,
                             double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZP1(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    // MAP2に入る外部値だけ補正する
    void DecodeWithZP2(const itpp::cvec& receivedSignal, itpp::bvec* output,
                       double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithZP2_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithZP2(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithZP2(const itpp::cvec& receivedSignal,
                             double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZP2(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }
    
    
    // ++++++++ Cyclic Suffix ++++++++
    void DecodeWithCyclicSuffix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCyclicSuffix_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithCyclicSuffix(receivedSignal, output, n0, numPads, iteration);
      } // else
    }
    itpp::bvec DecodeWithCyclicSuffix(const itpp::cvec& receivedSignal,
                                      double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCyclicSuffix(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    // MAP1に入る外部値だけ補正する
    void DecodeWithCS1(const itpp::cvec& receivedSignal, itpp::bvec* output,
                       double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCS1_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithCS1(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithCS1(const itpp::cvec& receivedSignal,
                             double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCS1(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    // MAP2に入る外部値だけ補正する
    void DecodeWithCS2(const itpp::cvec& receivedSignal, itpp::bvec* output,
                       double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCS2_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithCS2(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithCS2(const itpp::cvec& receivedSignal,
                             double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCS2(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    
    // ++++++++ Cyclic Prefix ++++++++
    void DecodeWithCyclicPrefix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCyclicPrefix_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithCyclicPrefix(receivedSignal, output, n0, numPads, iteration);
      } // else
    }
    itpp::bvec DecodeWithCyclicPrefix(const itpp::cvec& receivedSignal,
                                      double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCyclicPrefix(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    // MAP1に入る外部値だけ補正する
    void DecodeWithCP1(const itpp::cvec& receivedSignal, itpp::bvec* output,
                       double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCP1_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithCP1(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithCP1(const itpp::cvec& receivedSignal,
                             double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCP1(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    // MAP2に入る外部値だけ補正する
    void DecodeWithCP2(const itpp::cvec& receivedSignal, itpp::bvec* output,
                       double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCP2_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithCP2(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithCP2(const itpp::cvec& receivedSignal,
                             double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCP2(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    
    void DecodeWithInversedPrefix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                  double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithInversedPrefix_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithInversedPrefix(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithInversedPrefix(const itpp::cvec& receivedSignal,
                                        double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithInversedPrefix(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }
    
    void DecodeWithCyclicInfix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                               double n0, int start = 0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCyclicInfix_term(receivedSignal, output, n0, start, numPads, iteration);
      } // if
      else{
        doDecodeWithCyclicInfix(receivedSignal, output, n0, start, numPads, iteration);
      } // else 
    }
    itpp::bvec DecodeWithCyclicInfix(const itpp::cvec& receivedSignal, double n0, int start = 0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCyclicInfix(receivedSignal, &output, n0, start, numPads, iteration);
      return output;
    }
    
    static boost::rational< int > CodeRate()
    { return codeRate_; }

    static boost::rational< int > CodeRateWithTerm(int infoLength, int constraint)
    {
      int memory = constraint - 1;
      return boost::rational< int >(infoLength,
                                    codeRate_.denominator()*(infoLength + memory) + memory);
    }

    boost::rational< int > CodeRateWithTerm() const
    {
      return CodeRateWithTerm(interleaver_.size(), rsc1_.Constraint());
    }

    friend class ExitChart_TurboCode;
    
  };

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
  
  bool S_RandomInterleaver(itpp::ivec* output, int size, int spread, int maxTrial = 500)
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
