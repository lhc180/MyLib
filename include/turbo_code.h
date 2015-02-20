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
 * Last Updated: <2015/02/20 20:14:28 from alcohorhythm.local by yoshito>
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
  protected:                                       // ## privateにする!!
    const static boost::rational< int > codeRate_; // 符号化率は1/2で固定
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

    // マルチスレッド用
    // virtual void CalcAlpha(itpp::mat *alpha, const std::vector< itpp::mat > &gamma, int nodeNum) const;
    // virtual void CalcBeta(itpp::mat *beta, const std::vector< itpp::mat > &gamma, int nodeNum) const;
    
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
  protected:
    const itpp::ivec interleaver_;
    const Rsc rsc1_, rsc2_;
    static const boost::rational< int > codeRate_; // 符号化率は1/3で固定
    const int iteration_;
    const bool termination_;

    virtual void doEncode(const itpp::bvec& input, itpp::bvec* output) const; // NVI
    virtual void doEncode_term(const itpp::bvec& input, itpp::bvec* output) const;
    
    virtual void doDecode(const itpp::cvec& receivedSignal, itpp::bvec* output,
                          double n0) const;

    virtual void doDecode_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                               double n0) const;
    
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
    
    virtual void ModifyAPLLR(itpp::vec* llr, double APZeroProb, int start, int num) const;
    
    // virtual void ModifyLLRForZP(itpp::vec* llr, int numPads) const;

    virtual void Decoder(itpp::vec& llrToRsc1, const itpp::cvec& in1,
                         const itpp::cvec& in2,
                         double n0, int iteration) const;
    
    virtual void Decoder_term(itpp::vec& llrToRsc1, const itpp::cvec& in1, const itpp::cvec& in2,
                              double n0, int iteration) const;

    // virtual void DecoderForZP_term(itpp::vec& llrToRsc1, const itpp::cvec& in1, const itpp::cvec& in2,
    //                                double n0, int numPads, int iteration) const;
    
    virtual void ModifyLLRForCS(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForCP(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForIP(itpp::vec* llr, int numPads) const;

    virtual void ModifyLLRForCI(itpp::vec* llr, int start, int numPads) const;

    virtual void SeparateReceivedSignal(const itpp::cvec& receivedSignal,
                                        itpp::cvec* in1, itpp::cvec* in2) const;
    
  public:
    explicit TurboCode(const itpp::ivec &interleaver, int constraint = 3, int feedforward = 05, int feedback = 07,
                       int iteration = 10,
                       bool termination = true):
      interleaver_(interleaver), rsc1_(constraint, feedforward, feedback),
      rsc2_(constraint, feedforward, feedback), iteration_(iteration),
      termination_(termination)
    { }
        
    virtual ~TurboCode()
    { }

    // ここから全てNVIパターンを使った実装
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
    void Decode(const itpp::cvec& receivedSignal, itpp::bvec* output, double n0) const
    {
      if (termination_){
        doDecode_term(receivedSignal, output, n0);
      } // if termination_
      else {
        doDecode(receivedSignal, output, n0); 
      } 
    }
    itpp::bvec Decode(const itpp::cvec& receivedSignal, double n0) const
    {
      itpp::bvec output;
      Decode(receivedSignal, &output, n0);
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

  // Zero Padding
  class TurboCodeWithZP: public TurboCode
  {
  protected:
    mutable int numPads_;

    virtual void DecoderForZP(itpp::vec& llrToRsc1, const itpp::cvec& in1, const itpp::cvec& in2,
                              double n0, int iteration) const;
    virtual void DecoderForZP_term(itpp::vec& llrToRsc1, const itpp::cvec& in1, const itpp::cvec& in2,
                                   double n0, int iteration) const;
    
    virtual void ModifyLLR(itpp::vec* llr) const;

    // Encoderは普通のTurboCodeと同じやつで大丈夫なのでNVIで実際に呼ばれる関数だけ変える
    virtual void doDecode(const itpp::cvec &receivedSignal, itpp::bvec *output, double n0) const;
    virtual void doDecode_term(const itpp::cvec &receivedSignal, itpp::bvec *output,
                               double n0) const;
    
    virtual void doDecode_ModOne(const itpp::cvec &receivedSignal, itpp::bvec *output,
                                 int MAPIndex, double n0) const;
    virtual void doDecode_term_ModOne(const itpp::cvec &receivedSignal, itpp::bvec *output, 
                                      int MAPIndex, double n0) const;
    
    
  public:
    TurboCodeWithZP(const itpp::ivec &interleaver, int constraint, int feedforward, int feedback,
                    int iteration, int numPads, bool termination):
      TurboCode(interleaver, constraint, feedforward, feedback, iteration, termination),
      numPads_(numPads)
    { }
    
    virtual ~TurboCodeWithZP()
    { }

    void Decode_ModOne(const itpp::cvec& receivedSignal, itpp::bvec* output, int MAPIndex, double n0) const
    {
      if (termination_){
        doDecode_term_ModOne(receivedSignal, output, MAPIndex, n0);
      } // if termination_
      else {
        doDecode_ModOne(receivedSignal, output, MAPIndex, n0); 
      } 
    }
    itpp::bvec Decode_ModOne(const itpp::cvec& receivedSignal, int MAPIndex, double n0) const
    {
      itpp::bvec output;
      Decode_ModOne(receivedSignal, &output, MAPIndex, n0);
      return output;
    }
  };

  // With Decision of Zero Padding Insertion
  class TurboCodeWithZP_Judge: public TurboCodeWithZP
  {
  private:
    itpp::ivec numPadsCandidates_;
    itpp::ivec judgeBits_;
    int levels_;
    int secondIteration_;

  protected:
    virtual void doDecode(const itpp::cvec &receivedSignal, itpp::bvec *output, double n0) const;
    virtual void doDecode_term(const itpp::cvec &receivedSignal, itpp::bvec *output,
                               double n0) const;

    virtual bool checkPaddingInsertion(const itpp::cvec &receivedSignal, double n0) const;
    virtual bool checkPaddingInsertion_term(const itpp::cvec &receivedSignal, double n0) const;
    
  public:
    // 1段階しか無い場合
    TurboCodeWithZP_Judge(const itpp::ivec& interleaver, int constraint, int feedforward, int feedback,
                          int firstIteration, int secondIteration, int numPads, int judgeBits, bool termination):
      TurboCodeWithZP(interleaver, constraint, feedforward, feedback, firstIteration, 0, termination),
      numPadsCandidates_(1),judgeBits_(1),levels_(1), secondIteration_(secondIteration)
    {
      numPadsCandidates_[0] = numPads;
      judgeBits_[0] = judgeBits;

      assert(numPads >= judgeBits);
    }
    // 複数段階チェックする場合
    TurboCodeWithZP_Judge(const itpp::ivec& interleaver, int constraint, int feedforward, int feedback,
                          int firstIteration, int secondIteration, 
                          const itpp::ivec &numPads, const itpp::ivec &judgeBits, bool termination):
      TurboCodeWithZP(interleaver, constraint, feedforward, feedback, firstIteration, 0, termination),
      numPadsCandidates_(numPads), judgeBits_(judgeBits), levels_(numPads.size()), secondIteration_(secondIteration)
    {
      assert(numPads.size() == judgeBits.size());
    }
    virtual ~TurboCodeWithZP_Judge()
    { }

    bool isDecidedInsertion(const itpp::cvec &receivedSignal, double n0) const
    {
      
      if (termination_){
        return checkPaddingInsertion_term(receivedSignal, n0);
      } // if
      else{
        return checkPaddingInsertion(receivedSignal, n0);
      } // else
    }
    
  };

  
  // Turbo code with sparse zero padding bits
  class TurboCodeWithSZP: public TurboCodeWithZP
  {
  protected:
    mutable itpp::ivec padsPositions_;
    virtual void ModifyLLR(itpp::vec *llr) const;
    
  public:
    // 一定間隔でnull bitを入れたとき
    TurboCodeWithSZP(const itpp::ivec &interleaver, int constraint, int feedforward, int feedback,
                     int iteration, int numPads, bool termination):
      TurboCodeWithZP(interleaver, constraint, feedforward, feedback, iteration, numPads, termination),
      padsPositions_(numPads)
    {
      if (numPads != 0){
        int interval = std::floor(static_cast< double >(interleaver_.size())/static_cast< double >(numPads));
        for (int i = 0; i < numPads; ++i){
          padsPositions_[i] = i*interval;
        } // for i
      } // if 
    }
    // 特殊なビット位置での場合はまだ対応していない
    virtual ~TurboCodeWithSZP() { }
  };
  

  class TurboCodeWithSZP_Judge: public TurboCodeWithSZP
  {
  private:
    itpp::ivec judgeBits_;
    int levels_;
    int secondIteration_;
    std::vector< itpp::ivec > multiPadsPositions_;
    
  protected:
    virtual void doDecode(const itpp::cvec &receivedSignal, itpp::bvec *output, double n0) const;
    virtual void doDecode_term(const itpp::cvec &receivedSignal, itpp::bvec *output,
                               double n0) const;

    virtual bool checkPaddingInsertion(const itpp::cvec &receivedSignal, double n0) const;
    virtual bool checkPaddingInsertion_term(const itpp::cvec &receivedSignal, double n0) const;
    
  public:
    // 1段階しか無い場合
    // しかも一定間隔のSZPにのみ対応
    TurboCodeWithSZP_Judge(const itpp::ivec& interleaver, int constraint, int feedforward, int feedback,
                           int firstIteration, int secondIteration, int numPads, int judgeBits, bool termination):
      TurboCodeWithSZP(interleaver, constraint, feedforward, feedback, firstIteration, 0, termination),
      judgeBits_(1), levels_(1), secondIteration_(secondIteration), multiPadsPositions_(1, itpp::ivec(numPads))
    {
      judgeBits_[0] = judgeBits;

      if (numPads != 0){
        int interval = std::floor(static_cast< double >(interleaver_.size())/static_cast< double >(numPads));
        for (int i = 0; i < numPads; ++i){
          multiPadsPositions_[0][i] = i*interval;
        } // for i
      } // if       
      assert(numPads >= judgeBits);
    }

    TurboCodeWithSZP_Judge(const itpp::ivec& interleaver, int constraint, int feedforward, int feedback,
                           int firstIteration, int secondIteration, 
                           const std::vector< itpp::ivec >& multiPadsPositions, const itpp::ivec &judgeBits,
                           bool termination):
      TurboCodeWithSZP(interleaver, constraint, feedforward, feedback, firstIteration, 0, termination),
      judgeBits_(judgeBits), levels_(judgeBits.size()),
      secondIteration_(secondIteration), multiPadsPositions_(multiPadsPositions)
    {
      assert(static_cast< int >(multiPadsPositions.size()) == judgeBits.size());
    }

    virtual ~TurboCodeWithSZP_Judge()
    { }

    bool isDecidedInsertion(const itpp::cvec &receivedSignal, double n0) const
    {      
      if (termination_){
        return checkPaddingInsertion_term(receivedSignal, n0);
      } // if
      else{
        return checkPaddingInsertion(receivedSignal, n0);
      } // else
    }

  };  
  

  // Encoder Side
  class ZeroPadding
  {
  private:
    const int frameLength_;
  

  protected:
    itpp::ivec padsPositions_;
    virtual void SetupPadsPositions_() = 0;
    
  public:
    ZeroPadding(int frameLength): frameLength_(frameLength)
    {
    }
    
    virtual ~ZeroPadding();

    int FrameLength() const { return frameLength_; } 

    itpp::ivec PadsPositions() { return padsPositions_; } // ゲッタ

    virtual void SetNumPads(int numPads) = 0;
    
    // input.size() + numPads_のサイズのデータを返す
    virtual itpp::bvec Pad(const itpp::bvec& input) const = 0;

    // input.size()のデータを返す
    // つまりinputの中のデータを0で置き換える
    virtual itpp::bvec Nullify(const itpp::bvec& input) const = 0;
  };
  
  class CZP: public ZeroPadding
  {
  protected:
    virtual void SetupPadsPositions_();
    
  public:
    CZP(int frameLength, int numPads_ = 0): ZeroPadding(frameLength)
    {
      SetupPadsPositions_();
    }
    virtual ~CZP();

    virtual void SetNumPads(int numPads) {
      padsPositions_ = itpp::ivec(numPads);
      SetupPadsPositions_();
    }
  };
  

  // ++++ Decoder side ++++
  class LLR_Modifier
  {
  public:
    LLR_Modifier();
    virtual ~LLR_Modifier();
    
    virtual bool JudgeZP(const itpp::bvec& input) = 0;
    
    // virtual itpp::vec ModifyLLR(const itpp::vec& llr, double replacedLLR = -50)
    // {
    //   itpp::vec outputLLR(llr);
    //   for (int i = 0; i < padsPositions_.size(); ++i){
    //     outputLLR[padsPositions_[i]] = replacedLLR;
    //   } // for i
    //   return outputLLR;
    // }
  };
}



#endif
