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
 * Last Updated: <2014/04/22 14:15:57 from dr-yst-no-pc.local by yoshito>
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
    virtual void doDecodeWithZP(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                         double n0, int numPads, int iteration) const;

    virtual void doDecodeWithZP_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                              double n0, int numPads, int iteration) const;

    virtual void doDecodeWithZP_Judge(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                      double n0, int numPads, int numJudgeBits, int firstIteration,
                                      int secondIteration) const;

    virtual void doDecodeWithZP_Judge_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                           double n0, int numPads, int numJudgeBits, int firstIteration,
                                           int secondIteration) const;

    virtual void doDecodeWithZP_Judge(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                      double n0, const std::vector< int >& numPads,
                                      const std::vector< int > &numJudgeBits,
                                      int firstIteration, int secondIteration) const;

    virtual void doDecodeWithZP_Judge_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                           double n0, const std::vector< int >& numPads,
                                           const std::vector< int >& numJudgeBits, int firstIteration,
                                           int secondIteration) const;

    
    // ++++ Cyclic Suffix ++++
    virtual void doDecodeWithCS(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                          double n0, int numPads, int iteration) const;

    virtual void doDecodeWithCS_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                               double n0, int numPads, int iteration) const;

    
    // ++++ Cyclic Prefix ++++
    virtual void doDecodeWithCP(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                          double n0, int numPads, int iteration) const;

    virtual void doDecodeWithCP_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                               double n0, int numPads, int iteration) const;
    
    // ++++ Inversed Prefix ++++
    virtual void doDecodeWithIP(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                            double n0, int numPads, int iteration) const;

    virtual void doDecodeWithIP_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                                 double n0, int numPads, int iteration) const;

    // ++++ Cyclic Infix ++++
    virtual void doDecodeWithCI(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                         double n0, int start, int numPads, int iteration) const;
    
    virtual void doDecodeWithCI_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                              double n0, int start, int numPads, int iteration) const;
    
  protected:
    virtual void ModifyLLRForZP(itpp::vec* llr, int numPads) const;
    
    virtual void ModifyLLRForCS(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForCP(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForIP(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForCI(itpp::vec* llr, int start, int numPads) const;

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
    void DecodeWithZP(const itpp::cvec& receivedSignal, itpp::bvec* output,
                               double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithZP_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithZP(receivedSignal, output, n0, numPads, iteration);
      } // else 
    }
    itpp::bvec DecodeWithZP(const itpp::cvec& receivedSignal,
                                     double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZP(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }    

    
    void DecodeWithZP_Judge(const itpp::cvec& receivedSignal, itpp::bvec* output,
                            double n0, int numPads, int numJudgeBits, int firstIteration = 10,
                            int secondIteration = 10) const
    {
      if (termination_){
        doDecodeWithZP_Judge_term(receivedSignal, output, n0, numPads, numJudgeBits, firstIteration, secondIteration);
      } // if
      else{
        doDecodeWithZP_Judge(receivedSignal, output, n0, numPads, numJudgeBits, firstIteration, secondIteration);
      } // else 
    }

    itpp::bvec DecodeWithZP_Judge(const itpp::cvec& receivedSignal, double n0,
                                  int numPads, int numJudgeBits, int firstIteration = 10,
                                  int secondIteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZP_Judge(receivedSignal, &output, n0, numPads, numJudgeBits, firstIteration, secondIteration);
      return output;
    }

    void DecodeWithZP_Judge(const itpp::cvec& receivedSignal, itpp::bvec* output,
                            double n0, const std::vector< int >& numPads, std::vector< int >& numJudgeBits,
                            int firstIteration = 10, int secondIteration = 10) const
    {
      if (termination_){
        doDecodeWithZP_Judge_term(receivedSignal, output, n0, numPads, numJudgeBits,
                                  firstIteration, secondIteration);
      } // if
      else{
        doDecodeWithZP_Judge(receivedSignal, output, n0, numPads, numJudgeBits,
                             firstIteration, secondIteration);
      } // else
    }

    itpp::bvec DecodeWithZP_Judge(const itpp::cvec& receivedSignal, double n0,
                                  const std::vector< int >& numPads, std::vector< int >& numJudgeBits,
                                  int firstIteration = 10, int secondIteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZP_Judge(receivedSignal, &output, n0, numPads, numJudgeBits, firstIteration, secondIteration);
      return output;
    }
    
    // ++++++++ Cyclic Suffix ++++++++
    void DecodeWithCS(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCS_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithCS(receivedSignal, output, n0, numPads, iteration);
      } // else
    }
    itpp::bvec DecodeWithCS(const itpp::cvec& receivedSignal,
                                      double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCS(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }
    
    // ++++++++ Cyclic Prefix ++++++++
    void DecodeWithCP(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCP_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithCP(receivedSignal, output, n0, numPads, iteration);
      } // else
    }
    itpp::bvec DecodeWithCP(const itpp::cvec& receivedSignal,
                                      double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCP(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }
    
    void DecodeWithIP(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                  double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithIP_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        doDecodeWithIP(receivedSignal, output, n0, numPads, iteration);
      } // else
    }

    itpp::bvec DecodeWithIP(const itpp::cvec& receivedSignal,
                                        double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithIP(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }
    
    void DecodeWithCI(const itpp::cvec& receivedSignal, itpp::bvec* output,
                               double n0, int start = 0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        doDecodeWithCI_term(receivedSignal, output, n0, start, numPads, iteration);
      } // if
      else{
        doDecodeWithCI(receivedSignal, output, n0, start, numPads, iteration);
      } // else 
    }
    itpp::bvec DecodeWithCI(const itpp::cvec& receivedSignal, double n0, int start = 0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCI(receivedSignal, &output, n0, start, numPads, iteration);
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
