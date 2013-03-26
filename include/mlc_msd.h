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
namespace mylib{
  class MlcMsd{
  protected:
    std::vector<Ldpc> vecLDPC;
    itpp::Modulator_2D Modulator;
    int nLevel;			// bitsPerSymbolと一緒
    unsigned nCodeLength;
    double avrCodeRate;		// シンボルレート
    bool bSetDone;

  public:
    MlcMsd() : nLevel(0), nCodeLength(0), avrCodeRate(0), bSetDone(false)
    { }
    MlcMsd(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod);

    void Set(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod)
    {
      // vecLDPC.resize(vLDPC.size());
      //   for(int i = 0; i < vLDPC.size(); i++){
      //     &vecLDPC[i] = dynamic_cast< cLDPC_forMSD >(&vLDPC[i]);
      // }
      vecLDPC = vLDPC;

      Modulator = Mod;
  
      assert(static_cast<int>(vecLDPC.size()) == Modulator.bits_per_symbol());

      nLevel = vecLDPC.size();
  
      nCodeLength = vecLDPC[0].CodeLength();
      for(int level = 1; level < nLevel; level++){
        assert(nCodeLength == vecLDPC[level].CodeLength());
      } 

      avrCodeRate = 0.0;
      for(int level = 0; level < nLevel; level++){
        avrCodeRate += vecLDPC[level].CodeRate();
      }

      bSetDone = true;

    }

    void Encode(const std::vector<itpp::bvec> &vBits, itpp::cvec &symbols);

    itpp::cvec Encode(const std::vector<itpp::bvec> &vBits)
    {
      itpp::cvec symbols;
      Encode(vBits, symbols);
      return symbols;
    }
  
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
      std::vector< itpp::bvec > vDecodedBits(Modulator.bits_per_symbol());
      Decode(symbols, vDecodedBits, N0, loopMax, errConc);
      return vDecodedBits;
    }
  
    std::vector< int > InfoLengths()
    {
      std::vector< int > lengths(nLevel);
      for (int i = 0; i < nLevel; ++i){
        lengths[i] = vecLDPC[i].InfoLength();
      } // for i

      return lengths;
    }
    
    double SymbolRate(){
      return avrCodeRate;
    }
  };

  /*************** end of MLC_MSD ********************/

  // MLC用にLDPCをセットする
  // 平均符号化率を返す
  double SetLDPC_MLC(std::vector< Ldpc > &vecLDPC,
                     unsigned nCodeLength,
                     itpp::ivec &vecRowWeight,
                     itpp::ivec &vecColWeight);

  // セットパーティショニングマッピングにする
  void SetModulatorNatural(itpp::Modulator_2D &Mod);

  // Proposals for MLC and MSD below.
  // These adjust lengths of info bits such that higher level info bits move to
  // lower level.
  std::vector< itpp::bvec > AdjustInfoLength_forMLC(const std::vector< itpp::bvec > &infoBits,
                                                    const std::vector< int > &infoLengths_LDPC);

  std::vector< itpp::bvec > AdjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo,
                                                     const std::vector< int > &infoLengths);

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
