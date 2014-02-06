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
 * Last Updated: <2014/02/06 16:05:16 from Yoshitos-iMac.local by yoshito>
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
    mutable itpp::vec lambda_;
    mutable int lastState_;
    mutable itpp::bvec tailbitTable_;
    
  protected:
    virtual double Jacobian(double x1, double x2) const;

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
                      itpp::vec* logLikelihood_out, double n0) const;

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
  };
  
  class TurboCode
  {
  private:
    const itpp::ivec interleaver_;
    const Rsc rsc1_, rsc2_;
    static const boost::rational< int > codeRate_; // 符号化率は1/3で固定
    const bool termination_;

    virtual void DoEncode(const itpp::bvec& input, itpp::bvec* output) const; // NVI
    virtual void DoEncodeWithTerm(const itpp::bvec& input, itpp::bvec* output) const;
    
    virtual void DoDecode(const itpp::cvec& receivedSignal, itpp::bvec* output,
                          double n0, int iteration) const;

    virtual void DoDecodeWithTerm(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                  double n0, int iteration) const;
    
    
    virtual void DoDecodeWithZeroPadding(const itpp::cvec& receivedSignal, itpp::bvec* ouptut,
                                         double n0, int numPads, int iteration) const;

    virtual void DoDecodeWithZeroPadding_term(const itpp::cvec& receivedSignal, itpp::bvec* ouptut,
                                              double n0, int numPads, int iteration) const;

    virtual void DoDecodeWithCyclicSuffix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                          double n0, int numPads, int iteration) const;

    virtual void DoDecodeWithCyclicSuffix_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                               double n0, int numPads, int iteration) const;

    
    virtual void DoDecodeWithCyclicPrefix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                          double n0, int numPads, int iteration) const;

    virtual void DoDecodeWithCyclicPrefix_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                               double n0, int numPads, int iteration) const;


    
    
  protected:
    virtual void ModifyLLRForZeroPadding(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForCyclicSuffix(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForCyclicPrefix(itpp::vec* llr, int numPads) const;

    virtual void SeparateReceivedSignal(const itpp::cvec& receivedSignal,
                                        itpp::cvec* in1, itpp::cvec* in2) const;

    template< typename kind >
    itpp::Vec< kind > Interleave(const itpp::Vec< kind >& input) const;
    
    template< typename kind >
    itpp::Vec< kind > Deinterleave(const itpp::Vec< kind >& input) const;
    
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
        DoEncodeWithTerm(input, output);
      } // if termination_
      else {
        DoEncode(input, output);
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
        DoDecodeWithTerm(receivedSignal, output, n0, iteration);
      } // if termination_
      else {
        DoDecode(receivedSignal, output, n0, iteration); 
      } 
    }
    itpp::bvec Decode(const itpp::cvec& receivedSignal, double n0, int iteration = 10) const
    {
      itpp::bvec output;
      Decode(receivedSignal, &output, n0, iteration);
      return output;
    }

    void DecodeWithZeroPadding(const itpp::cvec& receivedSignal, itpp::bvec* output,
                               double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        DoDecodeWithZeroPadding_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        DoDecodeWithZeroPadding(receivedSignal, output, n0, numPads, iteration);
      } // else
    }
    itpp::bvec DecodeWithZeroPadding(const itpp::cvec& receivedSignal,
                                     double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZeroPadding(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    void DecodeWithCyclicSuffix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        DoDecodeWithCyclicSuffix_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        DoDecodeWithCyclicSuffix(receivedSignal, output, n0, numPads, iteration);
      } // else
    }
    itpp::bvec DecodeWithCyclicSuffix(const itpp::cvec& receivedSignal,
                                      double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCyclicSuffix(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }

    void DecodeWithCyclicPrefix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                double n0, int numPads = 0, int iteration = 10) const
    {
      if (termination_){
        DoDecodeWithCyclicPrefix_term(receivedSignal, output, n0, numPads, iteration);
      } // if
      else{
        DoDecodeWithCyclicPrefix(receivedSignal, output, n0, numPads, iteration);
      } // else
    }
    itpp::bvec DecodeWithCyclicPrefix(const itpp::cvec& receivedSignal,
                                      double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithCyclicPrefix(receivedSignal, &output, n0, numPads, iteration);
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
    
  };

  // class TurboCodeWithTerm: public TurboCode
  // {
  // private:
  //   virtual void DoEncode(const itpp::bvec& input, itpp::bvec* output) const; // NVI

  //   virtual void DoDecode(const itpp::cvec& receivedSignal, itpp::bvec* output,
  //                         double n0, int iteration) const;

  //   virtual void DoDecodeWithZeroPadding(const itpp::cvec& receivedSignal, itpp::bvec* ouptut,
  //                                        double n0, int numPads, int iteration) const;

  //   virtual void DoDecodeWithCyclicSuffix(const itpp::cvec& receivedSignal, itpp::bvec* output,
  //                                         double n0, int numPads, int iteration) const;

  //   virtual void DoDecodeWithCyclicPrefix(const itpp::cvec& receivedSignal, itpp::bvec* output,
  //                                         double n0, int numPads, int iteration) const;
    
    
  // public:
  //   explicit TurboCodeWithTerm(itpp::ivec interleaver, int constraint = 3, int feedforward = 05, int feedback = 07):
  //     TurboCode(interleaver, constraint, feedforward, feedback)
  //   { }
    
  //   virtual ~TurboCodeWithTerm()
  //   { }

  // };
  
}

#endif
