/************************************************************************************
 * turbo_code.cpp
 *   
 * turbo_code.hの実装部分
 *
 * Contents:
 *   class Rsc
 *   class TurboCode
 *
 * Last Updated: <2015/02/28 17:14:56 from alcohorhythm.local by yoshito>
 ************************************************************************************/
// #include <boost/thread.hpp>
#include "../include/myutl.h"
#include "../include/turbo_code.h"
#include "../include/myinterleavers.h"


namespace mylib{

  static const double LLR_THRESHOLD = 50;
  static const int JACOBIAN_TABLE_SIZE = 50000;
  static const double JACOBIAN_TABLE_SCALE = 0.005;
  
  const boost::rational< int > Rsc::codeRate_(1, 2);
  
  Rsc::Rsc(int constraint, unsigned int feedforward, unsigned int feedback):
    constraint_(constraint), memory_(constraint-1), stateNum_(static_cast< int >(itpp::pow2(memory_))),
    feedforward_(feedforward), feedback_(feedback),
    encodeTable_(static_cast< int >(itpp::pow2(memory_)), std::vector< encodeTable >(2)),
    revEncodeTable_(static_cast< int >(itpp::pow2(memory_)), std::vector< int >(2) ),
    tailbitTable_(static_cast< int >(itpp::pow2(memory_))),
    jacobianTable_(JACOBIAN_TABLE_SIZE),
    lastState_(-1)
  {
    assert(constraint_ > 0);
    assert(itpp::pow2(constraint_) > feedforward_ && itpp::pow2(constraint_) > feedback_);
    
    int mask = stateNum_ - 1;
    
    // テーブルを作る
    for (int state = 0; state < stateNum_; ++state){ // 
      for (int bit = 0; bit < 2; ++bit){
        unsigned int recursive = state & (feedback_ >> 1); // 入力を考慮していないfeedbackビット
        itpp::bin recursiveBit(0);
        for (int i = 0; i < memory_; ++i){
          recursiveBit += (recursive >> i) & 1;
        } // for i
        tailbitTable_[state] = recursiveBit;

        itpp::bin t_bit = recursiveBit + itpp::bin(((feedback_ & 1) & bit));

        // 以下いらなくなったけど一応残しておく
        //unsigned int t = ((state << 1) | bit) & feedback_; // メモリ内のフィードッバックで使うものの和
        //        itpp::bin t_bit(0);                      // tの各ビットのxor
        // for (int i = 0; i < constraint_; ++i){
        //   t_bit += (t >> i) & 1;
        // } // for i
          
        unsigned int out = ((state << 1) | static_cast< int >(t_bit)) & feedforward_;

        // outの各ビットのxor
        encodeTable_[state][bit].output_ = 0;
        for (int i = 0; i < constraint_; ++i){
          encodeTable_[state][bit].output_ += (out >> i) & 1;
        } // for i
        int nextState = ((state << 1) | static_cast< int >(t_bit)) & mask;
        encodeTable_[state][bit].nextState_ = nextState;
        revEncodeTable_[nextState][bit] = state;
      } // for bit
    } // for state

    for (int i = 0; i < JACOBIAN_TABLE_SIZE; ++i){
      jacobianTable_[i] = std::log(1.0 + std::exp(-static_cast< double >(i)*JACOBIAN_TABLE_SCALE));
    } // for i
  }
  
  void Rsc::GenParity(const itpp::bvec& input, itpp::bvec *output) const
  {
    output->set_size(input.size());

    int state = 0;
    for (int i = 0; i < input.size(); ++i){
      int bit = input[i];
      (*output)[i] = encodeTable_[state][bit].output_;
      state = encodeTable_[state][bit].nextState_;
    } // for i

    lastState_ = state;
  }

  void Rsc::Terminate(itpp::bvec *tailbits, itpp::bvec *parity) const
  {
    assert(lastState_ >= 0);
    
    tailbits->set_size(memory_);
    parity->set_size(memory_);

    int state = lastState_;
    for (int i = 0; i < memory_; ++i){
      int tailbit = tailbitTable_[state];
      (*tailbits)[i] = tailbit;
      (*parity)[i] = encodeTable_[state][tailbit].output_;
      state = encodeTable_[state][tailbit].nextState_;
    } // for i
    lastState_ = state;
    assert(lastState_ == 0);
    
  }

  void Rsc::GenParity_term(const itpp::bvec &input, itpp::bvec *tailbits, itpp::bvec *output) const
  {
    GenParity(input, output);

    itpp::bvec tailparity;
    Terminate(tailbits, &tailparity);
            
    (*output) = itpp::concat(*output, tailparity);
  }
  
  void Rsc::Encode(const itpp::bvec& input, itpp::bvec* output) const
  {
    itpp::bvec parity;
    GenParity(input, &parity);

    output->set_size(input.size()*2);
    for (int i = 0; i < input.size(); ++i){
      output[2*i] = input[i];
      output[2*i + 1] = parity[i];
    } // for i
  }

  void Rsc::Encode_term(const itpp::bvec &input, itpp::bvec *output) const
  {
    itpp::bvec parity, tailbits;
    GenParity_term(input, &tailbits, &parity);
    itpp::bvec inputWithTailbits = itpp::concat(input, tailbits);
    
    int systematicSize = inputWithTailbits.size();
    output->set_size(systematicSize*2);
    for (int i = 0; i < systematicSize; ++i){
      output[2*i] = inputWithTailbits[i];
      output[2*i + 1] = parity[i];
    } // for i
  }
  
  
  inline double Rsc::Jacobian(double x1, double x2) const
  {
    double y = std::max(x1, x2);

    // double temp = y + std::log(1.0 + std::exp(-std::abs(x2-x1)));    

    double temp = y + jacobianTable_[std::abs(itpp::round_i(
                                              itpp::SISO::threshold((x1-x2)/JACOBIAN_TABLE_SCALE,
                                                                    JACOBIAN_TABLE_SIZE)))];
    
    return temp;
  }

  // void Rsc::CalcAlpha(itpp::mat *alpha, const std::vector<itpp::mat> &gamma, int nodeNum) const
  // {
  //   for (int i = 1; i < nodeNum; ++i){
  //     for (int state = 0; state < stateNum_; ++state){
  //       // for (int bit = 0; bit < 2; ++bit){
  //       int primState0 = revEncodeTable_[state][0];
  //       int primState1 = revEncodeTable_[state][1];
  //       (*alpha)(i, state) = Jacobian((*alpha)(i-1, primState0) + gamma[i-1](primState0, 0),
  //                                     (*alpha)(i-1, primState1) + gamma[i-1](primState1, 1));
  //       // } // for bit
  //     } // for state
  //   } // for i
  // }

  // void Rsc::CalcBeta(itpp::mat *beta, const std::vector<itpp::mat> &gamma, int nodeNum) const
  // {
  //   for (int i = nodeNum - 2; i >= 0; --i){
  //     for (int state = 0; state < stateNum_; ++state){
  //       // for (int bit = 0; bit < 2; ++bit){
  //       int nextState0 = encodeTable_[state][0].nextState_;
  //       int nextState1 = encodeTable_[state][1].nextState_;
  //       (*beta)(i, state) = Jacobian((*beta)(i+1, nextState0) + gamma[i](state, 0),
  //                                    (*beta)(i+1, nextState1) + gamma[i](state, 1));
  //       // } // for bit
  //     } // for state
  //   } // for i
  // }
  
  void Rsc::CalcLambda(const itpp::cvec &received, const itpp::vec &logLikelihood_in, double n0,
                       bool knowLastState) const
  {
    itpp::BPSK_c bpsk;
    const int branchNum = received.size() * codeRate_.numerator() / codeRate_.denominator(); // レートは1/2固定
    const int nodeNum = branchNum + 1;
    
    itpp::mat alpha(nodeNum, stateNum_), beta(nodeNum, stateNum_);
    std::vector< itpp::mat > gamma(branchNum, itpp::mat(stateNum_, 2)); // [a](b,c)で時刻aの状態bに入力c

    lambda_.set_size(branchNum);          // デコードし終わったときのLambda
    
    itpp::mat logPrioriProb(branchNum, 2);
    
    // p.162の下部
    for (int i = 0; i < branchNum; ++i){
      // double t_exp = std::exp(logLikelihood_in[i]);
      logPrioriProb(i, 0) = -Jacobian(0, logLikelihood_in[i]);
      // std::log(1.0 + t_exp);
      logPrioriProb(i, 1) = logLikelihood_in[i] + logPrioriProb(i, 0);
    } // for i
    
    alpha.zeros();
    beta.zeros();
    
    // for (int i = 0; i < nodeNum; ++i){
    for (int state = 0; state < stateNum_; ++state){
      alpha(0, state) = -mylib::INFTY;
      beta(nodeNum - 1, state) = -mylib::INFTY;
    } // for state
    // } // for i
    alpha(0,0) = 0;

    if (knowLastState){                  
      beta(nodeNum - 1, lastState_) = 0;
    } // if
    else{
      for (int state = 0; state < stateNum_; ++state){
        beta(nodeNum - 1, state) = 0;
      } // for state
    } // else 
    
    
    // gamma
    {
      itpp::bvec originalCode(codeRate_.denominator()); // 元の送信信号
      for (int i = 0; i < branchNum; ++i){
        for (int state = 0; state < stateNum_; ++state){
          for (int bit = 0; bit < 2; ++bit){
            originalCode[0] = bit;
            originalCode[1] = encodeTable_[state][bit].output_;
            itpp::cvec originalTransSymbol = bpsk.modulate_bits(originalCode);

            double distance = itpp::sqr(std::real(originalTransSymbol[0] - received[2*i])) +
              itpp::sqr(std::real(originalTransSymbol[1] - received[2*i + 1]));
          
            gamma[i](state, bit) = logPrioriProb(i, bit) - distance / n0;
          } // for bit
        } // for state
      } // for i
    }   // gamma

    // boost::thread alphaThread = boost::thread(boost::bind(&Rsc::CalcAlpha, this, &alpha, gamma, nodeNum));
    // boost::thread betaThread =  boost::thread(boost::bind(&Rsc::CalcBeta, this, &beta, gamma, nodeNum));
    // alphaThread.join();
    // betaThread.join();
    
    // alpha
    for (int i = 1; i < nodeNum; ++i){
      for (int state = 0; state < stateNum_; ++state){
        // for (int bit = 0; bit < 2; ++bit){
        int primState0 = revEncodeTable_[state][0];
        int primState1 = revEncodeTable_[state][1];
        alpha(i, state) = Jacobian(alpha(i-1, primState0) + gamma[i-1](primState0, 0),
                                   alpha(i-1, primState1) + gamma[i-1](primState1, 1));
        // } // for bit
      } // for state
    } // for i
    
    // beta
    for (int i = nodeNum - 2; i >= 0; --i){
      for (int state = 0; state < stateNum_; ++state){
        // for (int bit = 0; bit < 2; ++bit){
        int nextState0 = encodeTable_[state][0].nextState_;
        int nextState1 = encodeTable_[state][1].nextState_;
        beta(i, state) = Jacobian(beta(i+1, nextState0) + gamma[i](state, 0),
                                  beta(i+1, nextState1) + gamma[i](state, 1));
          
        // } // for bit
      } // for state
    } // for i

    
    
    // lambda_
    for (int i = 1; i < nodeNum; ++i){
      itpp::vec likelihood(2);
      likelihood[0] = likelihood[1] = -mylib::INFTY;
      for (int state = 0; state < stateNum_; ++state){
        for (int bit = 0; bit < 2; ++bit){
          int nextState = encodeTable_[state][bit].nextState_;
          double t_delta = alpha(i-1, state) + gamma[i-1](state, bit)
            + beta(i, nextState);
          likelihood[bit] = Jacobian(likelihood[bit], t_delta);
        } // for bit
      } // for state
      lambda_[i - 1] = likelihood[1] - likelihood[0];
    } // for i
  }
  
  void Rsc::CalcLLR_out(const itpp::cvec &received, const itpp::vec &logLikelihood_in,
                        itpp::vec *logLikelihood_out, double n0) const
  {
    const int branchNum = received.size() * codeRate_.numerator() / codeRate_.denominator(); // レートは1/2固定
    
    logLikelihood_out->set_size(branchNum);
    for (int i = 0; i < branchNum; ++i){
      (*logLikelihood_out)[i] = lambda_[i] + 4.0 / n0 * std::real(received[2*i]) - logLikelihood_in[i];
      // if (std::abs((*logLikelihood_out)[i]) > LLR_THRESHOLD){
      //   (*logLikelihood_out)[i] = itpp::sign((*logLikelihood_out)[i])*LLR_THRESHOLD;
      // } // if 
    } // for i

  }

  // static
  void Rsc::CalcLLR_out(const itpp::cvec &received, const itpp::vec &lambda,
                        const itpp::vec &logLikelihood_in, itpp::vec *logLikelihood_out, double n0)
  {
    const int branchNum = received.size() * codeRate_.numerator() / codeRate_.denominator(); // レートは1/2固定
    
    logLikelihood_out->set_size(branchNum);
    for (int i = 0; i < branchNum; ++i){
      (*logLikelihood_out)[i] = lambda[i] + 4.0 / n0 * std::real(received[2*i]) - logLikelihood_in[i];
      // if (std::abs((*logLikelihood_out)[i]) > LLR_THRESHOLD){
      //   (*logLikelihood_out)[i] = itpp::sign((*logLikelihood_out)[i])*LLR_THRESHOLD;
      // } // if 
    } // for i

  }
  
  void Rsc::Decode(const itpp::cvec &received, const itpp::vec &logLikelihood_in,
                   itpp::vec *logLikelihood_out, double n0, bool knowLastState) const
  {
    // assert(lastState_ != -1);
    if (knowLastState && lastState_ == -1){
      lastState_ = 0;
    } // if
    CalcLambda(received, logLikelihood_in, n0, knowLastState);
    CalcLLR_out(received, logLikelihood_in, logLikelihood_out, n0);    
  }
  
  void Rsc::HardDecision(itpp::bvec *outputBits) const
  {
    int length = lambda_.size();
    
    outputBits->set_size(length);

    for (int i = 0; i < length; ++i){
      if (lambda_[i] > 0){
        (*outputBits)[i] = 1;
      } // if
      else{
        (*outputBits)[i] = 0;
      } // else 
    } // for i
    
  }

  itpp::bvec Rsc::HardDecision(const itpp::vec &lambda)
  {
    int length = lambda.size();

    itpp::bvec outputBits(length);

    for (int i = 0; i < length; ++i){
      if (lambda[i] > 0){
        outputBits[i] = 1;
      } // if
      else{
        outputBits[i] = 0;
      } // else 
    } // for i

    return outputBits;
  }

  /************************************************************************************
   * TurboCode 
   ************************************************************************************/

  const boost::rational< int > TurboCode::codeRate_(1, 3);
  
    
  void TurboCode::doEncode(const itpp::bvec &input, itpp::bvec *output) const
  {
    itpp::bvec parity1 = rsc1_.GenParity(input);

    itpp::bvec interleaved = Interleave(input, interleaver_);

    itpp::bvec parity2 = rsc2_.GenParity(interleaved);

    output->set_size(3*input.size());
    for (int i = 0; i < input.size(); ++i){
      (*output)[3*i] = input[i];
      (*output)[3*i + 1] = parity1[i];
      (*output)[3*i + 2] = parity2[i];
    } // for i
    
  }

  void TurboCode::SeparateReceivedSignal(const itpp::cvec &receivedSignal, itpp::cvec *in1, itpp::cvec *in2) const
  {
    int block = interleaver_.size();
    
    itpp::cvec r(block), parity1(block), parity2(block);

    for (int i = 0; i < block; ++i){
      r[i]       = receivedSignal[3*i];
      parity1[i] = receivedSignal[3*i + 1];
      parity2[i] = receivedSignal[3*i + 2];
    } // for i
    
    itpp::cvec interleaved_r = Interleave(r, interleaver_);
    
    in1->set_size(2*block);
    in2->set_size(2*block);
    for (int i = 0; i < block; ++i){
      (*in1)[2*i]     = r[i];
      (*in1)[2*i + 1] = parity1[i];

      (*in2)[2*i]     = interleaved_r[i];
      (*in2)[2*i + 1] = parity2[i];
    } // for i
    
  }
  
  void TurboCode::doDecode(const itpp::cvec &receivedSignal, itpp::bvec *output, double n0) const
  {
    assert(receivedSignal.size() % codeRate_.denominator() == 0);

    
    itpp::cvec in1, in2;
    SeparateReceivedSignal(receivedSignal, &in1, &in2);
    
    itpp::vec llrToRsc1(interleaver_.size());
    llrToRsc1.zeros();               // ## ここで提案法入れられるかも

    Decoder(&llrToRsc1, in1, in2, n0, iteration_);

    itpp::bvec interleaved_output = rsc2_.HardDecision();
    (*output) = Deinterleave(interleaved_output, interleaver_);
        
  }
  
  void TurboCode::doEncode_term(const itpp::bvec &input, itpp::bvec *output) const
  {
    itpp::bvec parity1, tailbits1, tailParity1;
    rsc1_.GenParity(input, &parity1);
    rsc1_.Terminate(&tailbits1, &tailParity1);
    
    itpp::bvec interleaved = Interleave(input, interleaver_);
    itpp::bvec parity2, tailbits2, tailParity2;
    rsc2_.GenParity(interleaved, &parity2);
    rsc2_.Terminate(&tailbits2, &tailParity2);
    
    itpp::bvec t_output(3*interleaver_.size());
    for (int i = 0; i < interleaver_.size(); ++i){
      t_output[3*i] = input[i];
      t_output[3*i + 1] = parity1[i];
      t_output[3*i + 2] = parity2[i];
    } // for i

    itpp::bvec codeTail1(2*tailbits1.size());
    for (int i = 0; i < tailbits1.size(); ++i){
      codeTail1[2*i] = tailbits1[i];
      codeTail1[2*i + 1] = tailParity1[i];
    } // for i
    *output = itpp::concat(t_output, codeTail1);
    
    itpp::bvec codeTail2(2*tailbits2.size());
    for (int i = 0; i < tailbits2.size(); ++i){
      codeTail2[2*i] = tailbits2[i];
      codeTail2[2*i + 1] = tailParity2[i];
    } // for i
    *output = itpp::concat(*output, codeTail2);
  }

  void TurboCode::Decoder(itpp::vec *llrToRsc1, const itpp::cvec &in1,
                          const itpp::cvec &in2, double n0, int iteration) const
  {
    for (int ite = 0; ite < iteration; ++ite){
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, *llrToRsc1, &llrFromRsc1, n0);
      
      itpp::vec llrToRsc2 = Interleave(llrFromRsc1, interleaver_);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);
      
      *llrToRsc1 = Deinterleave(llrFromRsc2, interleaver_);
    } // for ite
    
  }

  void TurboCode::Decoder_term(itpp::vec *llrToRsc1, const itpp::cvec &in1, const itpp::cvec &in2,
                                      double n0, int iteration) const
  {
    int memory = rsc1_.Constraint() - 1;
    static itpp::vec llrZeros(0);
    if (llrZeros.size() != memory){
      llrZeros.set_size(memory);
      llrZeros.zeros();
    } // if

    for (int ite = 0; ite < iteration; ++ite){
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, *llrToRsc1, &llrFromRsc1, n0);

      itpp::vec llrToRsc2 = Interleave(llrFromRsc1.left(interleaver_.size()), interleaver_);
      llrToRsc2 = itpp::concat(llrToRsc2, llrZeros);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);

      *llrToRsc1 = Deinterleave(llrFromRsc2.left(interleaver_.size()), interleaver_);
      *llrToRsc1 = itpp::concat(*llrToRsc1, llrZeros);
      
    } // for ite

  }
  
  void TurboCode::doDecode_term(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                   double n0) const
  {
    int memory = rsc1_.Constraint() - 1;
    
    itpp::cvec in1, in2;
    SeparateReceivedSignal(receivedSignal, &in1, &in2);

    itpp::cvec tail1 = receivedSignal.mid(3*interleaver_.size(), 2*memory);
    itpp::cvec tail2 = receivedSignal.right(2*memory);

    in1 = itpp::concat(in1, tail1);
    in2 = itpp::concat(in2, tail2);

    itpp::vec llrToRsc1(interleaver_.size() + memory);
    llrToRsc1.zeros();

    Decoder_term(&llrToRsc1, in1, in2, n0, iteration_);
    
    itpp::bvec interleaved_output = rsc2_.HardDecision();
    (*output) = Deinterleave(interleaved_output.left(interleaver_.size()), interleaver_);
  }
  
  /************************************************************************************
   * TurboCodeWithZP 
   * 
   * Implementation of Turbo Code with Zero Padding
   ************************************************************************************/
  void TurboCodeWithZP::Decoder(itpp::vec *llrToRsc1, const itpp::cvec &in1,
                                     const itpp::cvec &in2, double n0, int iterations) const
  {
    for (int ite = 0; ite < iterations; ++ite){
      itpp::vec llrToRsc1_mod = zeroPadding_.ModifyLLR(*llrToRsc1);
      
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, llrToRsc1_mod, &llrFromRsc1, n0);

      itpp::vec llrFromRsc1_mod = zeroPadding_.ModifyLLR(llrFromRsc1);
      
      itpp::vec llrToRsc2 = Interleave(llrFromRsc1_mod.left(interleaver_.size()), interleaver_);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);

      *llrToRsc1 = Deinterleave(llrFromRsc2.left(interleaver_.size()), interleaver_);
    } // for ite
    
  }
  
  // void TurboCodeWithZP::doDecode(const itpp::cvec& receivedSignal, itpp::bvec* output,
  //                                         double n0) const
  // {
  //   assert(receivedSignal.size() % codeRate_.denominator() == 0);
    
  //   itpp::cvec in1, in2;
  //   SeparateReceivedSignal(receivedSignal, &in1, &in2);
    
  //   itpp::vec llrToRsc1(interleaver_.size());
  //   llrToRsc1.zeros();

  //   this->Decoder(&llrToRsc1, in1, in2, n0, iteration_);
    
  //   itpp::bvec interleaved_output = rsc2_.HardDecision();
  //   (*output) = Deinterleave(interleaved_output, interleaver_);
      
  // }

    
  void TurboCodeWithZP::Decoder_term(itpp::vec *llrToRsc1, const itpp::cvec& in1, const itpp::cvec& in2,
                                          double n0, int iterations) const
  {
    // std::cout << "## zeroPadding_ = " << zeroPadding_.PadPositions() << std::endl;

    int memory = rsc1_.Constraint() - 1;
    static itpp::vec llrZeros(0);
    if (llrZeros.size() != memory){
      llrZeros.set_size(memory);
      llrZeros.zeros();
    } // if
    
    for (int ite = 0; ite < iterations; ++ite){
      //      itpp::vec llrToRsc1_mod = zeroPadding_.ModifyLLR(*llrToRsc1);
      
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, *llrToRsc1, &llrFromRsc1, n0); // ##
      
      // itpp::vec llrFromRsc1_mod = zeroPadding_.ModifyLLR(llrFromRsc1);
            
      itpp::vec llrToRsc2 = Interleave(llrFromRsc1.left(interleaver_.size()), interleaver_); // ##
      llrToRsc2 = itpp::concat(llrToRsc2, llrZeros);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);

      *llrToRsc1 = Deinterleave(llrFromRsc2.left(interleaver_.size()), interleaver_);
      *llrToRsc1 = itpp::concat(*llrToRsc1, llrZeros);
    } // for ite
    
  }
  
  // void TurboCodeWithZP::doDecode_term(const itpp::cvec &receivedSignal, itpp::bvec *output,
  //                                        double n0) const
  // {    
  //   int memory = rsc1_.Constraint() - 1;
    
  //   itpp::cvec in1, in2;
  //   SeparateReceivedSignal(receivedSignal, &in1, &in2);

  //   itpp::cvec tail1 = receivedSignal.mid(3*interleaver_.size(), 2*memory);
  //   itpp::cvec tail2 = receivedSignal.right(2*memory);

  //   in1 = itpp::concat(in1, tail1);
  //   in2 = itpp::concat(in2, tail2);

  //   itpp::vec llrToRsc1(interleaver_.size() + memory);
  //   llrToRsc1.zeros();
    
  //   this->Decoder_term(&llrToRsc1, in1, in2, n0, iteration_);
    
  //   itpp::bvec interleaved_output = rsc2_.HardDecision();
  //   (*output) = Deinterleave(interleaved_output.left(interleaver_.size()), interleaver_);

  // }

  /************************************************************************************
   * ZeroPadding 
   ************************************************************************************/
  itpp::bvec ZeroPadding::Pad(const itpp::bvec& input) const
  {
    itpp::bvec output(input);
    for (int i = 0; i < padPositions_.size(); ++i){
      output.ins(padPositions_[i], itpp::bin(0));
    } // for 

    assert(output.size() == frameLength_);
    return output;
  }
  
  itpp::bvec ZeroPadding::Nullify(const itpp::bvec& input) const
  {
    assert(input.size() == frameLength_);
    
    itpp::bvec output(input);
    for (int i = 0; i < padPositions_.size(); ++i){
      output[padPositions_[i]] = itpp::bin(0);
    } // for i
    return output;
  }

  bool ZeroPadding::JudgeZP(const itpp::bvec& input, int threshold) const
  {
    int numZeros = 0;
    
    for (int i = 0; i < padPositions_.size(); ++i){
      numZeros += static_cast< int >(!input[padPositions_[i]]);
    } // for i

    if (numZeros >= threshold){
      return true;
    } // if
    else{
      return false;
    } // else 
  }

  itpp::vec ZeroPadding::ModifyLLR(const itpp::vec &llr, double replacedLLR) const
  {
    itpp::vec output(llr);
    
    for (int i = 0; i < padPositions_.size(); ++i){
      output[padPositions_[i]] = replacedLLR;
    } // for i

    return output;
  }
  
  void CZP::SetupPadPositions_(int numPads)
  {
    itpp::ivec padPositions(numPads);
    int frameLength = ZeroPadding::FrameLength();
    
    for (int i = 0; i < numPads; ++i){
      padPositions[i] = frameLength - numPads + i;
    } // for i

    ZeroPadding::SetPadPositions(padPositions);
  }

  void SZI::SetupPadPositions_(int numPads)
  {
    itpp::ivec padPositions(numPads);
    if (numPads != 0){
      int frameLength = ZeroPadding::FrameLength();
      int padInterval = std::floor(static_cast< double >(frameLength)/static_cast< double >(numPads));
    
      for (int i = 0; i < numPads; ++i){
        padPositions[i] = i * padInterval;
      } // for i
    } // if numPads
    ZeroPadding::SetPadPositions(padPositions);
  }  
}
