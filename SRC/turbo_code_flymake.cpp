/************************************************************************************
 * turbo_code.cpp
 *   
 * turbo_code.hの実装部分
 *
 * Contents:
 *   class Rsc
 *   class TurboCode
 *
 * Last Updated: <2014/02/05 16:43:57 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/
#include "../include/myutl.h"
#include "../include/turbo_code.h"

namespace mylib{

  const boost::rational< int > Rsc::codeRate_(1, 2);
  
  Rsc::Rsc(int constraint, unsigned int feedforward, unsigned int feedback):
    constraint_(constraint), memory_(constraint-1), stateNum_(static_cast< int >(itpp::pow2(memory_))),
    feedforward_(feedforward), feedback_(feedback),
    encodeTable_(static_cast< int >(itpp::pow2(memory_)), std::vector< encodeTable >(2)),
    lastState_(-1)
  {
    assert(constraint_ > 0);
    assert(itpp::pow2(constraint_) > feedforward_ && itpp::pow2(constraint_) > feedback_);
      
    // std::cout << "## feedforward = " << feedforward_ << std::endl;
    // std::cout << "## feedback = " << feedback_ << std::endl;
      
    int mask = stateNum_ - 1;

    tailbitTable_.set_size(stateNum_);
    
    // テーブルを作る
    // ##
    // std::cout << "state\tinput\toutput" << std::endl;
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
        encodeTable_[state][bit].nextState_ = ((state << 1) | static_cast< int >(t_bit)) & mask;
        // std::cout << state << '\t' << bit << '\t' << encodeTable_[state][bit].output_ << std::endl; // ##
          
      } // for bit
    } // for state

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

  void Rsc::GenParityWithTerm(const itpp::bvec &input, itpp::bvec *tailbits, itpp::bvec *output) const
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

  void Rsc::EncodeWithTerm(const itpp::bvec &input, itpp::bvec *output) const
  {
    itpp::bvec parity, tailbits;
    GenParityWithTerm(input, &tailbits, &parity);
    itpp::bvec inputWithTailbits = itpp::concat(input, tailbits);
    
    int systematicSize = inputWithTailbits.size();
    output->set_size(systematicSize*2);
    for (int i = 0; i < systematicSize; ++i){
      output[2*i] = inputWithTailbits[i];
      output[2*i + 1] = parity[i];
    } // for i
  }
  
  
  // void Rsc::MapDecode(const itpp::cvec &received, const itpp::vec &logLikelihood_in,
  //                     itpp::vec *logLikelihood_out, double n0)
  // {
  //   assert(lastState_ != -1);

  //   itpp::BPSK bpsk;
    
  //   const int branchNum = received.size() /codesPerInfoBit_; // レートは1/2固定
  //   const int nodeNum = branchNum + 1;
    
  //   itpp::mat prioriProb(branchNum, 2);
  //   itpp::mat alpha(nodeNum, stateNum_), beta(nodeNum, stateNum_);
  //   std::vector< itpp::mat > gamma(branchNum, itpp::mat(stateNum_, 2)); // [a](b,c)で時刻aの状態bに入力c

  //   lambda_.set_size(branchNum);

  //   // p.162の下部
  //   for (int i = 0; i < branchNum; ++i){
  //     double t_exp = std::exp(logLikelihood_in[i]);
  //     prioriProb(i, 0) = 1.0 / (1.0 + t_exp);
  //     prioriProb(i, 1) = t_exp / (1.0 + t_exp); 
  //   } // for i

  //   alpha.zeros();
  //   alpha(0,0) = 1;
  //   beta.zeros();
  //   beta(nodeNum - 1, lastState_);
    
  //   for (int i = 0; i < branchNum; ++i){
  //     for (int state = 0; state < stateNum_; ++state){
  //       for (int bit = 0; bit < 2; ++bit){
  //         itpp::bvec originalCode(codesPerInfoBit_); // 元の送信信号
  //         originalCode[0] = bit;
  //         originalCode[1] = encodeTable_[state][bit].output_;
  //         itpp::vec originalTransSymbol = bpsk.modulate_bits(originalCode);

  //         double distance = itpp::sqr(originalTransSymbol[0] - received[2*i]) +
  //           itpp::sqr(originalTransSymbol[1] - received[2*i +1]);
          
  //         gamma[i](state, bit) = prioriProb(i, bit) * std::exp(-distance / n0);
  //       } // for bit
  //     } // for state
  //   } // for i

  //   // alpha
  //   for (int i = 1; i < nodeNum; ++i){
  //     double norm = 0;
  //     for (int state = 0; state < stateNum_; ++state){
  //       for (int bit = 0; bit < 2; ++bit){
  //         norm += alpha(i-1, state) * gamma[i-1](state, bit);
  //       } // for bit
  //     } // for state

  //     for (int state = 0; state < stateNum_; ++state){
  //       for (int bit = 0; bit < 2; ++bit){
  //         int nextState = encodeTable_[state][bit].nextState_;
  //         alpha(i, nextState) += alpha(i-1, state) * gamma[i-1](state, bit) / norm;
  //       } // for bit
  //     } // for state
  //   } // for i

  //   // beta
  //   for (int i = nodeNum - 2; i >= 0; --i){
  //     double norm = 0;
  //     for (int state = 0; state < stateNum_; ++state){
  //       for (int bit = 0; bit < 2; ++bit){
  //         int nextState = encodeTable_[state][bit].nextState_;
  //         norm += beta(i+1, nextState) * gamma[i](state, bit);
  //       } // for bit
  //     } // for state

  //     for (int state = 0; state < stateNum_; ++state){
  //       for (int bit = 0; bit < 2; ++bit){
  //         int nextState = encodeTable_[state][bit].nextState_;
  //         beta(i, state) += beta(i+1, nextState) * gamma[i](state, bit) / norm;
  //       } // for bit
  //     } // for state
  //   } // for i

  //   // lambda_
  //   for (int i = 1; i < nodeNum; ++i){
  //     itpp::vec likelihood(2);
  //     likelihood.zeros();
  //     for (int state = 0; state < stateNum_; ++state){
  //       for (int bit = 0; bit < 2; ++bit){
  //         int nextState = encodeTable_[state][bit].nextState_;
  //         likelihood[bit] += alpha(i-1, state) * gamma[i-1](state, bit) * beta(i, nextState);
  //        } // for bit
  //     } // for state
  //     lambda_[i - 1] = std::log(likelihood[1]) - std::log(likelihood[0]);
  //   } // for i

  //   // ## マッピングと値2の確認 ##
  //   for (int i = 0; i < branchNum; ++i){
  //     (*logLikelihood_out)[i] = lambda_[i] + 2.0 / n0 * std::real(received[2*i]) - logLikelihood_in[i];
  //   } // for i
    
  // }

  inline double Rsc::Jacobian(double x1, double x2) const
  {
    // ## max関数に置き換える
    double y = std::max(x1, x2);

    // std::cout << "## y = " << y << std::endl;
    
    double temp = y + std::log(1.0 + std::exp(-std::abs(x2-x1)));

    // std::cout << "## jac = " << temp << std::endl;
    
    return temp;
    
  }


  void Rsc::Decode(const itpp::cvec &received, const itpp::vec &logLikelihood_in,
                   itpp::vec *logLikelihood_out, double n0) const
  {
    assert(lastState_ != -1);

    itpp::BPSK_c bpsk;
    
    const int branchNum = received.size() * codeRate_.numerator() / codeRate_.denominator(); // レートは1/2固定
    const int nodeNum = branchNum + 1;
    
    itpp::mat alpha(nodeNum, stateNum_), beta(nodeNum, stateNum_);
    std::vector< itpp::mat > gamma(branchNum, itpp::mat(stateNum_, 2)); // [a](b,c)で時刻aの状態bに入力c

    lambda_.set_size(branchNum);          // デコードし終わったときのLambda

    
    itpp::mat logPrioriProb(branchNum, 2);
    // p.162の下部
    for (int i = 0; i < branchNum; ++i){
      double t_exp = std::exp(logLikelihood_in[i]);
      logPrioriProb(i, 0) = -std::log(1.0 + t_exp);
      logPrioriProb(i, 1) = logLikelihood_in[i] + logPrioriProb(i, 0);
    } // for i
    
    alpha.zeros();
    beta.zeros();

    for (int i = 0; i < nodeNum; ++i){
      for (int state = 0; state < stateNum_; ++state){
        alpha(i, state) = -mylib::INFTY;
        beta(i, state) = -mylib::INFTY;
      } // for state
    } // for i
    alpha(0,0) = 0;
    beta(nodeNum - 1, lastState_) = 0;

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
    
    // alpha
    for (int i = 1; i < nodeNum; ++i){
      for (int state = 0; state < stateNum_; ++state){
        for (int bit = 0; bit < 2; ++bit){
          int nextState = encodeTable_[state][bit].nextState_;
          alpha(i, nextState) = Jacobian(alpha(i, nextState),
                                         alpha(i-1, state) + gamma[i-1](state, bit));
        } // for bit
      } // for state
    } // for i
    
    // beta
    for (int i = nodeNum - 2; i >= 0; --i){
      for (int state = 0; state < stateNum_; ++state){
        for (int bit = 0; bit < 2; ++bit){
          int nextState = encodeTable_[state][bit].nextState_;
          beta(i, state) = Jacobian(beta(i, state),
                                    beta(i+1, nextState) + gamma[i](state, bit));
          
        } // for bit
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

    // ## マッピングと値2の確認 ##

    logLikelihood_out->set_size(branchNum);
    for (int i = 0; i < branchNum; ++i){
      (*logLikelihood_out)[i] = lambda_[i] + 2.0 / n0 * std::real(received[2*i]) - logLikelihood_in[i];
    } // for i
    
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


  // +++++++++++++++ TurboCode +++++++++++++++++++

  const boost::rational< int > TurboCode::codeRate_(1, 3);
  
  // ## 後にテンプレートクラスにする
  template< typename kind >
  itpp::Vec< kind > TurboCode::Interleave(const itpp::Vec< kind >& input) const
  {
    assert(input.size() == interleaver_.size());
    itpp::Vec< kind > output(input.size());

    for (int i = 0; i < input.size(); ++i){
      output[i] = input[interleaver_[i]];
    } // for i

    return output;
  }

  template< typename kind >
  itpp::Vec< kind > TurboCode::Deinterleave(const itpp::Vec< kind >& input) const
  {
    assert(input.size() == interleaver_.size());
    itpp::Vec< kind > output(input.size());

    for (int i = 0; i < input.size(); ++i){
      output[interleaver_[i]] = input[i];
    } // for i

    return output;
  }
  
  
  void TurboCode::DoEncode(const itpp::bvec &input, itpp::bvec *output) const
  {
    itpp::bvec parity1 = rsc1_.GenParity(input);

    itpp::bvec interleaved = Interleave(input);

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
    
    itpp::cvec interleaved_r = Interleave(r);

    in1->set_size(2*block);
    in2->set_size(2*block);
    for (int i = 0; i < block; ++i){
      in1[2*i]     = r[i];
      in1[2*i + 1] = parity1[i];

      in2[2*i]     = interleaved_r[i];
      in2[2*i + 1] = parity2[i];
    } // for i

  }
  
  void TurboCode::DoDecode(const itpp::cvec &receivedSignal, itpp::bvec *output, double n0, int iteration) const
  {
    assert(receivedSignal.size() % codeRate_.denominator() == 0);

    itpp::cvec in1, in2;
    SeparateReceivedSignal(receivedSignal, &in1, &in2);
    
    itpp::vec llrToRsc1(interleaver_.size());
    llrToRsc1.zeros();               // ## ここで提案法入れられるかも
    
    for (int ite = 0; ite < iteration; ++ite){
      
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, llrToRsc1, &llrFromRsc1, n0);
      
      itpp::vec llrToRsc2 = Interleave(llrFromRsc1);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);
      
      llrToRsc1 = Deinterleave(llrFromRsc2);
    } // for ite

    itpp::bvec interleaved_output = rsc2_.HardDecision();
    (*output) = Deinterleave(interleaved_output);
        
  }

  void TurboCode::ModifyLLRForZeroPadding(itpp::vec *llr, int numPads) const
  {
    int length = interleaver_.size();
    int numEffectiveBits = length - numPads;
    
    for (int i = numEffectiveBits; i < length; ++i){
      (*llr)[i] = -30;
    } // for i
  }
    
  
  void TurboCode::DoDecodeWithZeroPadding(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                        double n0, int numPads, int iteration) const
  {
    assert(receivedSignal.size() % codeRate_.denominator() == 0);

    itpp::cvec in1, in2;
    SeparateReceivedSignal(receivedSignal, &in1, &in2);
    
    itpp::vec llrToRsc1(interleaver_.size());
    llrToRsc1.zeros();               // ## ここで提案法入れられるかも

    ModifyLLRForZeroPadding(&llrToRsc1, numPads);
    
    for (int ite = 0; ite < iteration; ++ite){
      
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, llrToRsc1, &llrFromRsc1, n0);

      ModifyLLRForZeroPadding(&llrFromRsc1, numPads); // ## いらないかも
      
      itpp::vec llrToRsc2 = Interleave(llrFromRsc1);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);
      
      llrToRsc1 = Deinterleave(llrFromRsc2);

      ModifyLLRForZeroPadding(&llrToRsc1, numPads); // ## いらないかも
    } // for ite

    itpp::bvec interleaved_output = rsc2_.HardDecision();
    (*output) = Deinterleave(interleaved_output);
      
  }

  void TurboCode::ModifyLLRForCyclicSuffix(itpp::vec *llr, int numPads) const
  {
    int infoLength = interleaver_.size();
    int numEffectiveBits = infoLength - numPads;
    
    for (int i = 0; i < numPads; ++i){
      (*llr)[numEffectiveBits - numPads + i] += (*llr)[numEffectiveBits + i];
      (*llr)[numEffectiveBits + i] = (*llr)[numEffectiveBits - numPads + i];
    } // for i

  }

  void TurboCode::DoDecodeWithCyclicSuffix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                      double n0, int numPads, int iteration) const
  {
    assert(receivedSignal.size() % codeRate_.denominator() == 0);

    itpp::cvec in1, in2;
    SeparateReceivedSignal(receivedSignal, &in1, &in2);
    
    itpp::vec llrToRsc1(interleaver_.size());
    llrToRsc1.zeros();
        
    for (int ite = 0; ite < iteration; ++ite){
      
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, llrToRsc1, &llrFromRsc1, n0);

      ModifyLLRForCyclicSuffix(&llrFromRsc1, numPads); // ## いらないかも
      
      itpp::vec llrToRsc2 = Interleave(llrFromRsc1);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);
      
      llrToRsc1 = Deinterleave(llrFromRsc2);

      ModifyLLRForCyclicSuffix(&llrToRsc1, numPads); // ## いらないかも
    } // for ite

    itpp::bvec interleaved_output = rsc2_.HardDecision();
    (*output) = Deinterleave(interleaved_output);

  }

  void TurboCode::ModifyLLRForCyclicPrefix(itpp::vec *llr, int numPads) const
  {
    int infoLength = interleaver_.size();
    int numEffectiveBits = infoLength - numPads;
    
    for (int i = 0; i < numPads; ++i){
      (*llr)[i] += (*llr)[numEffectiveBits + i];
      (*llr)[numEffectiveBits + i] = (*llr)[i];
    } // for i

  }

  void TurboCode::DoDecodeWithCyclicPrefix(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                      double n0, int numPads, int iteration) const
  {
    assert(receivedSignal.size() % codeRate_.denominator() == 0);
    
    itpp::cvec in1, in2;
    SeparateReceivedSignal(receivedSignal, &in1, &in2);
        
    itpp::vec llrToRsc1(interleaver_.size());
    llrToRsc1.zeros();
        
    for (int ite = 0; ite < iteration; ++ite){
      
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, llrToRsc1, &llrFromRsc1, n0);

      ModifyLLRForCyclicPrefix(&llrFromRsc1, numPads); // ## いらないかも
      
      itpp::vec llrToRsc2 = Interleave(llrFromRsc1);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);
      
      llrToRsc1 = Deinterleave(llrFromRsc2);

      ModifyLLRForCyclicPrefix(&llrToRsc1, numPads); // ## いらないかも
    } // for ite

    itpp::bvec interleaved_output = rsc2_.HardDecision();
    (*output) = Deinterleave(interleaved_output);

  }

  void TurboCode::DoEncodeWithTerm(const itpp::bvec &input, itpp::bvec *output) const
  {
    itpp::bvec parity1, tailbits1, tailParity1;
    rsc1_.GenParity(input, &parity1);
    rsc1_.Terminate(&tailbits1, &tailParity1);
    
    itpp::bvec interleaved = Interleave(input);
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
      codeTail2[2*i] = tailParity2[i];
    } // for i
    *output = itpp::concat(*output, codeTail2);
  }

  void TurboCode::DoDecodeWithTerm(const itpp::cvec& receivedSignal, itpp::bvec* output,
                                   double n0, int iteration) const
  {
    itpp::cvec in1, in2;
    SeparateReceivedSignal(receivedSignal, &in1, &in2);

    itpp::cvec tail1 = receivedSignal.mid(3*interleaver_.size(), 2*(rsc1_.Constraint()-1));
    itpp::cvec tail2 = receivedSignal.right(2*(rsc2_.Constraint()-1));

    in1 = itpp::concat(in1, tail1);
    in2 = itpp::concat(in2, tail2);

    itpp::vec llrToRsc1(interleaver_.size() + rsc1_.Constraint() - 1);
    llrToRsc1.zeros();
    itpp::vec llrFromRsc1_tail(rsc1_.Constraint()-1), llrFromRsc2_tail(rsc2_.Constraint()-1);
    
    for (int ite = 0; ite < iteration; ++ite){
      itpp::vec llrFromRsc1;
      rsc1_.Decode(in1, llrToRsc1, &llrFromRsc1, n0);

      llrFromRsc1_tail = llrFromRsc1.right(rsc1_.Constraint()-1);
      itpp::vec llrToRsc2 = Interleave(llrFromRsc1.left(interleaver_.size()));
      llrToRsc2 = itpp::concat(llrToRsc2, llrFromRsc2_tail);

      itpp::vec llrFromRsc2;
      rsc2_.Decode(in2, llrToRsc2, &llrFromRsc2, n0);
      
      llrFromRsc2_tail = llrFromRsc2.right(rsc2_.Constraint()-1);
      llrToRsc1 = Deinterleave(llrFromRsc2.left(interleaver_.size()));
      llrToRsc1 = itpp::concat(llrToRsc1, llrFromRsc1_tail);
    } // for ite

    itpp::bvec interleaved_output = rsc2_.HardDecision();
    (*output) = Deinterleave(interleaved_output.left(interleaver_.size()))
    
  }

  
}

















