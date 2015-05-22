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
 * Last Updated: <2015/05/22 11:21:53 from WatanabeYoshito-no-iMac.local by yoshito>
 ************************************************************************************/

#include <cassert>
#include <map>
#include <set>
#include <boost/rational.hpp>
#include <itpp/itcomm.h>

namespace mylib{
  struct encodeTable{
    itpp::bin output_;
    int nextState_;
  };
  
  class Rsc
  {
  protected:                                       // ## private�ˤ���!!
    const static boost::rational< int > codeRate_; // ��沽Ψ��1/2�Ǹ���
    const int constraint_;
    const int memory_;
    const int stateNum_;
    const unsigned int feedforward_;    // Octal Form
    const unsigned int feedback_;       // Octal Form
    std::vector< std::vector< encodeTable > > encodeTable_;
    std::vector< std::vector < int > > revEncodeTable_;
    mutable itpp::vec lambda_;
    mutable itpp::bvec tailbitTable_;
    mutable itpp::vec jacobianTable_;
    mutable int lastState_;
    
  protected:
    virtual double Jacobian(double x1, double x2) const;

    // �ޥ������å���
    // virtual void CalcAlpha(itpp::mat *alpha, const std::vector< itpp::mat > &gamma, int nodeNum) const;
    // virtual void CalcBeta(itpp::mat *beta, const std::vector< itpp::mat > &gamma, int nodeNum) const;
    
    virtual void CalcLambda(const itpp::cvec& received, const itpp::vec &logLikelihood_in,
                            double n0, bool knowLastState) const;

    virtual void CalcLLR_out(const itpp::cvec &received, const itpp::vec &logLikelihood_in,
                             itpp::vec *logLikelihood_out, double n0) const;

    // TurboCode��ǻȤ���
    static void CalcLLR_out(const itpp::cvec &received, const itpp::vec &lambda,
                            const itpp::vec &logLikelihood_in, itpp::vec *logLikelihood_out, double n0);

    static itpp::bvec HardDecision(const itpp::vec& lambda);
    
  public:
    // feedforward��feedback��8��ɽ���ʤΤ�0..�ǳ�Ǽ���Ƥ���
    // feedforward:g_1, feedback:g_0
    // constraint=5, feedforward=021, feedback=037
    Rsc(int constraint = 3, unsigned int feedforward = 05, unsigned int feedback = 07);
    
    virtual ~Rsc()                     // �Ѿ��Ϥ��ʤ�����
    { }
    
    // Terminate is not supported.
    void GenParity(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec GenParity(const itpp::bvec& input) const{
      itpp::bvec output; GenParity(input, &output); return output;
    }

    void Terminate(itpp::bvec *tailbits, itpp::bvec *parity) const;

    // Supporting Termination
    void GenParity_term(const itpp::bvec& input, itpp::bvec* tailbits, itpp::bvec* output) const;
    
    void Encode(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec Encode(const itpp::bvec& input) const{
      itpp::bvec output; Encode(input, &output); return output;
    }

    void Encode_term(const itpp::bvec& input, itpp::bvec* output) const;
    itpp::bvec Encode_term(const itpp::bvec& input) const{
      itpp::bvec output;
      Encode_term(input, &output);
      return output;
    }

    // Log-Map ����
    void Decode(const itpp::cvec& received, const itpp::vec& logLikelihood_in,
                itpp::vec* logLikelihood_out, double n0, bool knowLastState = true) const;
    
    // BPSK_c�ˤ����ޤ��б����Ƥ��ʤ�
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
  protected:
    const itpp::ivec interleaver_;
    const Rsc rsc1_, rsc2_;
    static const boost::rational< int > codeRate_; // ��沽Ψ��1/3�Ǹ���
    int iteration_;
    const bool termination_;
    mutable itpp::vec llrToRsc1_;

  protected:
    virtual void ResetLLR_() const;
    
    virtual void doEncode(const itpp::bvec& input, itpp::bvec* output) const; // NVI
    virtual void doEncode_term(const itpp::bvec& input, itpp::bvec* output) const;
    
    virtual void doDecode(const itpp::cvec& receivedSignal, itpp::bvec* output,
                          double n0) const;

    virtual void doDecode_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                               double n0) const;

    virtual void Decoder(const itpp::cvec& in1,
                         const itpp::cvec& in2,
                         double n0, int iteration) const;
    
    virtual void Decoder_term(const itpp::cvec& in1, const itpp::cvec& in2,
                              double n0, int iteration) const;
    
    // ���������ǥ�������1��2�Ѥ�ʬ����
    virtual void SeparateReceivedSignal(const itpp::cvec& receivedSignal,
                                        itpp::cvec* in1, itpp::cvec* in2) const;
    
  public:
    explicit TurboCode(const itpp::ivec &interleaver, int constraint = 3, int feedforward = 05, int feedback = 07,
                       int iteration = 10,
                       bool termination = true):
      interleaver_(interleaver), rsc1_(constraint, feedforward, feedback),
      rsc2_(constraint, feedforward, feedback), iteration_(iteration),
      termination_(termination)
    {
      ResetLLR_();
    }

    virtual ~TurboCode()
    { }

    void SetIterations(int iteration) { iteration_ = iteration; }
    
    // ������������NVI�ѥ������Ȥä�����
    void Encode(const itpp::bvec& input, itpp::bvec* output) const
    {
      if (termination_){
        doEncode_term(input, output);
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
    void Decode(const itpp::cvec& receivedSignal, itpp::bvec* output, double n0, bool reset = true) const
    {
      if (reset){
        ResetLLR_();
      } // if 
      
      if (termination_){
        doDecode_term(receivedSignal, output, n0);
      } // if termination_
      else {
        doDecode(receivedSignal, output, n0); 
      } 
    }
    itpp::bvec Decode(const itpp::cvec& receivedSignal, double n0, bool reset = true) const
    {
      itpp::bvec output;
      Decode(receivedSignal, &output, n0, reset);
      return output;
    }
    
    static boost::rational< int > CodeRate()
    { return codeRate_; }

    static boost::rational< int > CodeRate_term(int infoLength, int constraint)
    {
      int memory = constraint - 1;
      return boost::rational< int >(infoLength,
                                    codeRate_.denominator()*(infoLength + memory) + memory);
    }

    boost::rational< int > CodeRate_term() const
    {
      return CodeRate_term(interleaver_.size(), rsc1_.Constraint());
    }
    
  };

  // Key���ѥǥ��󥰥ӥåȤο�
  // Value������
  // ɬ����������Ȥ߹�碌����Ȥ�ʤ���Ф����ʤ�
  // ������ɬ�פˤʤä���Ŭ�����䤹
  const std::map< int, int > PresetPaddingThreshold {
    {0, 0},
      {64, 57},
        {128, 100},
          {256, 179},
            {512, 329},
  };
  
  class ZeroPadding
  {
  private:
    int frameLength_;
    std::set< int > padPositions_;
    
  public:
    ZeroPadding(int frameLength = 0, const std::set< int > padPositions = std::set<int>()): frameLength_(frameLength)
    {
      SetPadPositions(padPositions);
    }
    
    virtual ~ZeroPadding() { }

    int FrameLength() const { return frameLength_; } 
    void SetFrameLength(int frameLength) { frameLength_ = frameLength; }
    
    std::set< int > PadPositions() const { return padPositions_; } // ���å�

    // �����ǥ����Ȥ����
    virtual void SetPadPositions(const std::set< int >& padPositions)
    {
      padPositions_ = padPositions;
    }
    
    // ++++ Encoder side ++++
    // input.size() + numPads_�Υ������Υǡ������֤�
    virtual itpp::bvec Pad(const itpp::bvec& input) const;

    // input.size()�Υǡ������֤�
    // �Ĥޤ�input����Υǡ�����0���֤�������
    virtual itpp::bvec Nullify(const itpp::bvec& input) const;

    // ++++ Decoder side ++++
    virtual bool JudgeZP(const itpp::bvec& input) const;

    static bool JudgeZP(const itpp::bvec& input, const std::set< int > &padPositions);

    virtual itpp::vec ModifyLLR(const itpp::vec& llr, double replacedLLR = -50) const;
  };
  
  class CZP: public ZeroPadding
  {
  protected:
    virtual void SetupPadPositions_(int numPads);
    
  public:
    CZP(int frameLength = 0, int numPads = 0): ZeroPadding(frameLength)
    {
      SetupPadPositions_(numPads);
    }
    virtual ~CZP() { }

    virtual void SetNumPads(int numPads) { SetupPadPositions_(numPads); }
    
  };
  
  class SZI: public ZeroPadding
  {
  protected:
    virtual void SetupPadPositions_(int numPads);
    
  public:
    SZI(int frameLength = 0, int numPads = 0): ZeroPadding(frameLength)
    {
      SetupPadPositions_(numPads);
    }
    virtual ~SZI() { }

    virtual void SetNumPads(int numPads) {  SetupPadPositions_(numPads); }
  };
  
  // Zero Padding
  // �¼�Turbo Decoder
  class TurboCodeWithZP: public TurboCode
  {
  private:
    ZeroPadding zeroPadding_;
    
  protected:
    // �ʲ��ϴ��쥯�饹��virtual��������Ƥ��롣
    virtual void Decoder(const itpp::cvec& in1, const itpp::cvec& in2,
                              double n0, int iteration) const;
    virtual void Decoder_term(const itpp::cvec& in1, const itpp::cvec& in2,
                                   double n0, int iteration) const;

    // Encoder�����̤�TurboCode��Ʊ����Ĥ�����פʤΤ�NVI�Ǽºݤ˸ƤФ��ؿ������Ѥ���
    // �ʲ��ϴ��쥯�饹����Ȥ�Ʊ��
    // virtual void doDecode(const itpp::cvec &receivedSignal, itpp::bvec *output, double n0) const;
    // virtual void doDecode_term(const itpp::cvec &receivedSignal, itpp::bvec *output,
    //                            double n0) const;    
    
  public:
    TurboCodeWithZP(const itpp::ivec &interleaver, int constraint = 3, int feedforward = 05, int feedback = 07,
                    int iteration = 10, const ZeroPadding &zp = ZeroPadding(0), bool termination = true):
      TurboCode(interleaver, constraint, feedforward, feedback, iteration, termination),
      zeroPadding_(zp)
    { }
    
    virtual ~TurboCodeWithZP()
    { }

    void SetZeroPadding(const ZeroPadding &zp)
    {
      zeroPadding_ = zp;      
    }
  };

  // With Decision of Zero Padding Insertion
  // Note: This is only a DECODER.
  class TurboDecodeWithZP_Judge
  {
  private:
    TurboCodeWithZP turboCodeZP_;
    itpp::ivec iterations_;
    std::vector< ZeroPadding > zeroPaddings_; // �ǽ��0�ĤΥѥǥ��󥰥ӥå�
    mutable std::vector< std::set< int > > complementPadPositions_;
    
    std::set< int > Complement(const std::set< int > &superSet, const std::set< int > &subset)
    {
      std::set< int > setComplement(superSet);
      for (auto cur_elem = subset.begin(); cur_elem != subset.end(); ++cur_elem){
        setComplement.erase(*cur_elem);
      } // for cur_elem
      return setComplement;
    }
    
  protected:
    virtual void AssertionCheck_() const
    {
      assert(iterations_.size() == 2);
    }

    virtual void buildComplementPadPositions_() 
    {
      complementPadPositions_ = std::vector< std::set< int > >(zeroPaddings_.size());
      for (int i = 0; i < static_cast< int >(zeroPaddings_.size()); ++i){
        complementPadPositions_[i] = zeroPaddings_[i].PadPositions();
        for (int j = 0; j < i; ++j){
          complementPadPositions_[i] = Complement(complementPadPositions_[i], complementPadPositions_[j]); 
        }   // j

        // ####++++
        std::cout << "complementPadPositions_[" << i << "] = " << std::endl;
        for (auto cur_pos = complementPadPositions_[i].begin(); cur_pos != complementPadPositions_[i].end();
             ++cur_pos){
          std::cout << *cur_pos << ' ';
        } // for cur_pos
        std::cout << std::endl;
        // ####----
        
      }     // i
    }
    
  public:
    TurboDecodeWithZP_Judge(const itpp::ivec& interleaver, int constraint, int feedforward, int feedback,
                          const itpp::ivec& iterations, const std::vector< ZeroPadding >& zp,
                          bool termination):
      turboCodeZP_(interleaver, constraint, feedforward, feedback, iterations[0], zp[0], termination),
      iterations_(iterations), zeroPaddings_(zp)
    {
      AssertionCheck_();
      buildComplementPadPositions_();
    }
    
    TurboDecodeWithZP_Judge(const TurboCodeWithZP& turboCode,
                            const itpp::ivec& iterations, const std::vector< ZeroPadding >& zp):
      turboCodeZP_(turboCode),
      iterations_(iterations), zeroPaddings_(zp)
    {
      AssertionCheck_();
      buildComplementPadPositions_();
    }
    virtual ~TurboDecodeWithZP_Judge()
    { }

    void Decode(const itpp::cvec& receivedSignal, itpp::bvec* output, double n0)
    {      
      turboCodeZP_.SetIterations(iterations_[0]);
      turboCodeZP_.SetZeroPadding(zeroPaddings_[0]); // 0�ΤϤ�

      turboCodeZP_.Decode(receivedSignal, output, n0);

      int padCand_i = 0;
      for (int i = 1; i < static_cast< int >(zeroPaddings_.size()); ++i){
        if (ZeroPadding::JudgeZP(*output, complementPadPositions_[i])) { 
          padCand_i = i;
        }
        else {
          break;
        }
      } // for i
      
      turboCodeZP_.SetZeroPadding(zeroPaddings_[padCand_i]);
      turboCodeZP_.SetIterations(iterations_[1]);
      turboCodeZP_.Decode(receivedSignal, output, n0, false);
    }
    
    itpp::bvec Decode(const itpp::cvec& receivedSignal, double n0)
    {
      itpp::bvec output;
      Decode(receivedSignal, &output, n0);
      return output;
    }    
  };

  // Multiple padding �ξ��Ͼ��������ΰ��֤�ɬ���礭�����ΰ��֤˴ޤޤ�롣
  class PaddingBitInserter
  {
  private:
    std::vector< ZeroPadding > zeroPaddings_; // �ǽ��0��
    itpp::vec zpProbs_;

  protected:
    void AssertionCheck_() const
    {
      assert(static_cast< int >(zeroPaddings_.size()) == zpProbs_.size());
      double totalProb = 0.0;
      for (int i = 0; i < zpProbs_.size(); ++i){
        totalProb += zpProbs_[i];
      } // for i
      it_assert(totalProb == 1.0, "Sum of probabilities is not equal to 1.");
    }

    int IndexOfPadding_() const
    {
      double randNum = itpp::randu();
      double cummProbs = 0.0;

      int index = 0;
      for (int i = 0; i < zpProbs_.size(); ++i){
        cummProbs += zpProbs_[i];
        if (randNum <= cummProbs){
          index = i;
          break;
        } // if
      } // for i
      return index;
    }
    
  public:
    PaddingBitInserter(const std::vector< ZeroPadding >& zp, const itpp::vec& probs):
      zeroPaddings_(zp), zpProbs_(probs)
    {
      AssertionCheck_();
    }
    virtual ~PaddingBitInserter() { }

    // ���Ѥ��줿zeroPaddings_�Υ���ǥ������֤�
    int Pad(const itpp::bvec &input, itpp::bvec *output) const
    {
      int index = IndexOfPadding_();
      *output = zeroPaddings_[index].Pad(input);
      return index;
    }
    
    itpp::bvec Pad(const itpp::bvec &input) const
    {
      itpp::bvec output;
      Pad(input, &output);
      return output;
    }

    // ���Ѥ��줿zeroPaddings_�Υ���ǥ������֤�
    int Nullify(const itpp::bvec &input, itpp::bvec *output) const
    {
      int index = IndexOfPadding_();
      *output = zeroPaddings_[index].Nullify(input);
      return index;
    }
    
    itpp::bvec Nullify(const itpp::bvec &input) const
    {
      itpp::bvec output;
      Nullify(input, &output);
      return output;
    }
  };

}


#endif
