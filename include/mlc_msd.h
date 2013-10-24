#ifndef MLC_MSD_H
#define MLC_MSD_H

#include <vector>
#include "../include/myldpc.h"
#include "../include/myutl.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
/*************************************************/
// これはMultilevel Coding, Multistage Decodingのクラスである
// 符号化方法は今のところLDPCのみを想定している
// それぞれ別々のクラスで、実装する
//
/***********************************************/
/***************** to do list ******************/
// error concealmentに対応させる
// 
/**********************************************/


/***************** class cMLC_MSD *********************/
/************************************************************************************
 * MlcMsd 
 * 
 * ModはMSB first
 ************************************************************************************/

namespace mylib{
  class MlcMsd{
  protected:
    std::vector<Ldpc> vecLDPC_;
    itpp::Modulator_2D modulator_;
    int numLevels_;			// bitsPerSymbolと一緒
    unsigned codeLength_;
    double avrCodeRate_;		// シンボルレート
    bool setDone_;

    virtual void LevelDownModulator(std::vector< itpp::Modulator_2D >* vecMod, const itpp::bvec& decodedCode);
    
  public:
    MlcMsd() : numLevels_(0), codeLength_(0), avrCodeRate_(0), setDone_(false)
    { }
    MlcMsd(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod);

    virtual ~MlcMsd()
    {
      vecLDPC_.clear();
    }
    
    void Set(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod)
    {
      // vecLDPC.resize(vLDPC.size());
      //   for(int i = 0; i < vLDPC.size(); i++){
      //     &vecLDPC[i] = dynamic_cast< cLDPC_forMSD >(&vLDPC[i]);
      // }
      vecLDPC_ = vLDPC;

      modulator_ = Mod;
  
      assert(static_cast<int>(vecLDPC_.size()) == modulator_.bits_per_symbol());

      numLevels_ = vecLDPC_.size();
  
      codeLength_ = vecLDPC_[0].CodeLength();
      for(int level = 1; level < numLevels_; level++){
        assert(codeLength_ == vecLDPC_[level].CodeLength());
      } 

      avrCodeRate_ = 0.0;
      for(int level = 0; level < numLevels_; level++){
        avrCodeRate_ += vecLDPC_[level].CodeRate();
      }

      setDone_ = true;

    }

    void Encode(const std::vector<itpp::bvec> &vBits, itpp::cvec &symbols);

    itpp::cvec Encode(const std::vector<itpp::bvec> &vBits)
    {
      itpp::cvec symbols;
      Encode(vBits, symbols);
      return symbols;
    }

    // ## デバッグ用
    void EncodeTest(const std::vector< itpp::bvec >& vBits, itpp::cvec& symbols, std::vector< itpp::bvec >& codes);
    
    
    // 各レベルのイテレーション回数を返す
    // errConc = -1のときは隠蔽無し
    // 例えば1なら、レベル0は確実に出力、レベル1以降で誤りあればその時点で0で埋める
    // レベル以上の数字も-1と同じ
    itpp::ivec Decode(const itpp::cvec &symbols,
                      std::vector< itpp::bvec > &vDecodedBits,
                      double            N0,
                      unsigned          loopMax,
                      int              errConc = -1); // どのレベルでエラー隠蔽するか

    std::vector<itpp::bvec> Decode(const itpp::cvec &symbols,
                                   double N0,
                                   unsigned loopMax,
                                   int errConc = -1)
    {
      std::vector< itpp::bvec > vDecodedBits(modulator_.bits_per_symbol());
      Decode(symbols, vDecodedBits, N0, loopMax, errConc);
      return vDecodedBits;
    }

    // ##
    itpp::bvec DecodeTest(const itpp::cvec& symbols,
                          const std::vector< itpp::bvec >& knownBits,
                          double n0,
                          unsigned loopMax);
    
    
    std::vector< int > InfoLengths()
    {
      std::vector< int > lengths(numLevels_);
      for (int i = 0; i < numLevels_; ++i){
        lengths[i] = vecLDPC_[i].InfoLength();
      } // for i

      return lengths;
    }
    
    double SymbolRate(){
      return avrCodeRate_;
    }
  };

  /*************** end of MLC_MSD ********************/

  /************************************************************************************
   * MlcMsdForBP
   *
   * 8PSK BlockPartitioning用のMSD
   * レベル1と2をマルチスレッドで復号する
   * MlcMsdの派生クラスで作る
   ************************************************************************************/

  class MlcMsdWith8pskBp: public MlcMsd
  {
  public:
    // constructor
    MlcMsdWith8pskBp(): MlcMsd() { }
    MlcMsdWith8pskBp(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod)
    {
      this->Set(vLDPC, Mod);
    }
    
    // destructor
    virtual ~MlcMsdWith8pskBp() { }

    void Set(std::vector<Ldpc>& vLDPC, itpp::Modulator_2D &Mod)
    {
      if (Mod.bits_per_symbol() != 3){
        std::cerr << "Error: Mod.bits_per_symbol() != 3 in MlcMsdFor8pskBp." << std::endl;
        exit(8);
      } // if
      MlcMsd::Set(vLDPC, Mod);
    }
    
    std::vector< itpp::bvec > Decode(const itpp::cvec &symbols,
                                     double N0,
                                     unsigned loopMax,
                                     int errConc = -1);
    
  };
  
  // -------------------- MlcMsdForBP --------------------
  
  // MLC用にLDPCをセットする
  // 平均符号化率を返す
  double SetLdpcForMlc(std::vector< Ldpc > &vecLDPC_,
                     unsigned codeLength_,
                     itpp::ivec &vecRowWeight,
                     itpp::ivec &vecColWeight);

  // // セットパーティショニングマッピングにする
  // void SetModulatorNatural(itpp::Modulator_2D &Mod);

  /************************************************************************************
   * NaturalPSK 
   * 
   * セットパーティショニング(ナチュラル)マッピングのPSK
   ************************************************************************************/
  class SetPartitioningPSK: public itpp::PSK
  {
  protected:
    virtual void Init();
    
  public:
    explicit SetPartitioningPSK(int m): itpp::PSK(m)
    {
      Init();
    }
    virtual ~SetPartitioningPSK()
    { }
  };

  /************************************************************************************
   * ReverseSpPSK 
   * 
   * セットパーティショニングの逆マッピング
   ************************************************************************************/
  class ReverseSpPSK: public itpp::PSK
  {
  protected:
    virtual void Init();
    
  public:
    explicit ReverseSpPSK(int m): itpp::PSK(m)
    {
      Init();
    }
    virtual ~ReverseSpPSK()
    { }
  };
  
  
  // itpp::Modulatorのデフォルトのグレイマッピングがそのままblock partitioning
  class BlockPartitioningPSK: public itpp::PSK
  {
  public:
    BlockPartitioningPSK(int m): itpp::PSK(m)
    { }
    virtual ~BlockPartitioningPSK()
    { }
  };


  // Proposals for MLC and MSD below.
  // These adjust lengths of info bits such that higher level info bits move to
  // lower level.
  // std::vector< itpp::bvec > AdjustInfoLength_forMLC(const std::vector< itpp::bvec > &infoBits,
  //                                                   const std::vector< int > &infoLengths_LDPC);

  // std::vector< itpp::bvec > AdjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo,
  //                                                    const std::vector< int > &infoLengths);

  // end of proposals.

  /************************************************************************************
   * DivideIntoLevels -- MLC用に入力系列を複数のレベルに分割する。
   * 
   * Arguments:
   *   input -- 入力系列
   *   numSamples -- 各レベルの情報ビット数を格納
   *
   * Return Value:
   *   Vector_2D< kind > output -- レベル分けされた情報系列が各行に割り当てられている
   ************************************************************************************/
  template< typename kind >
  inline mylib::Vector_2D< kind > DivideIntoLevels(const std::vector< kind > &input, const std::vector< int > &numSamples)
  {
    int sum = mylib::Sum(numSamples);
    assert(input.size() == sum);

    mylib::Vector_2D< kind > output(numSamples.size());
    int inputIndex = 0;
    for(int level = 0; level < numSamples.size(); level++)
      {
        output(level).resize(numSamples[level]);
        for(int i = 0; i < numSamples[level]; i++)
          {
            output(level, i) = input[inputIndex];
            inputIndex++;
          }
      }

    return output;
  }
} // end of namespace mylib

#endif
