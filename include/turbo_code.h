#ifndef TURBO_CODE_H
#define TURBO_CODE_H

/************************************************************************************
 * turbo_code.h
 *   
 * ���������˴ؤ��륯�饹
 *
 * Contents:
 *   class Rsc
 *   class TurboCode
 *
 * Last Updated: <2014/01/23 18:32:19 from Yoshitos-iMac.local by yoshito>
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
    const static boost::rational< int > codeRate_; // ��沽Ψ��1/2�Ǹ���
    const int constraint_;
    const int memory_;
    const int stateNum_;
    const int feedforward_;    // Octal Form
    const int feedback_;       // Octal Form
    std::vector< std::vector< encodeTable > > encodeTable_;
    mutable itpp::vec lambda_;
    mutable int lastState_;

  protected:
    double Jacobian(double x1, double x2) const;
    
  public:
    // feedforward��feedback��8��ɽ���ʤΤ�0..�ǳ�Ǽ���Ƥ���
    // feedforward:g_1, feedback:g_0
    // constraint=5, feedforward=021, feedback=037
    Rsc(int constraint = 3, int feedforward = 05, int feedback = 07);
    
    virtual ~Rsc()                     // �Ѿ��Ϥ��ʤ�����
    { }
    
    // Terminate is not suported.
    void GenParity(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec GenParity(const itpp::bvec& input) const{
      itpp::bvec output; GenParity(input, &output); return output;
    }
    
    void Encode(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec Encode(const itpp::bvec& input) const{
      itpp::bvec output; Encode(input, &output); return output;
    }

    // complex��
    // void MapDecode(const itpp::cvec& received, const itpp::vec& logLikelihood_in,
    //                itpp::vec* logLikelihood_out, double n0);

    // Log-Map ����
    void Decode(const itpp::cvec& received, const itpp::vec& logLikelihood_in,
                      itpp::vec* logLikelihood_out, double n0) const;

    // BPSK_c�ˤ����ޤ��б����Ƥ��ʤ�
    void HardDecision(itpp::bvec* outputBits) const;
    itpp::bvec HardDecision() const
    {
      itpp::bvec outputBits;
      HardDecision(&outputBits);
      return outputBits;
    }

    static boost::rational< int > CodeRate(){
      return codeRate_;
    }
  };
  
  class TurboCode
  {
  private:
    const itpp::ivec interleaver_;
    const Rsc rsc1_, rsc2_;
    const static boost::rational< int > codeRate_; // ��沽Ψ��1/3�Ǹ���

  protected:
    virtual void ModifyLLRForZeroPadding(itpp::vec* llr, int numPads) const;
    
  public:
    explicit TurboCode(itpp::ivec interleaver, int constraint = 3, int feedforward = 05, int feedback = 07):
      interleaver_(interleaver), rsc1_(constraint, feedforward, feedback), rsc2_(constraint, feedforward, feedback)
    { }
        
    virtual ~TurboCode()
    { }

    void Encode(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec Encode(const itpp::bvec& input) const
    {
      itpp::bvec output;
      Encode(input, &output);
      return output;
    }

    // Log-Map Decode
    void Decode(const itpp::cvec& receivedSignal, itpp::bvec* output, double n0, int iteration = 10) const;
    itpp::bvec Decode(const itpp::cvec& receivedSignal, double n0, int iteration = 10) const
    {
      itpp::bvec output;
      Decode(receivedSignal, &output, n0, iteration);
      return output;
    }

    void DecodeWithZeroPadding(const itpp::cvec& receivedSignal, itpp::bvec* output,
                               double n0, int numPads = 0, int iteration = 10) const;
    itpp::bvec DecodeWithZeroPadding(const itpp::cvec& receivedSignal,
                                     double n0, int numPads = 0, int iteration = 10) const
    {
      itpp::bvec output;
      DecodeWithZeroPadding(receivedSignal, &output, n0, numPads, iteration);
      return output;
    }
    
    
    static boost::rational< int > CodeRate()
    {
      return codeRate_;
    }
    
  };

}

#endif