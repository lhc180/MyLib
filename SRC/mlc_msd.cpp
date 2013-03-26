#include <cassert>
#include "../include/mlc_msd2.h"

 
/******************* class MlcMsd *************************/
namespace mylib{
  MlcMsd::MlcMsd(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod)
  {
    Set(vLDPC, Mod);
  }

  void MlcMsd::Encode(const std::vector<itpp::bvec> &vBits,
                        itpp::cvec &vSymbol)
  {
    assert(bSetDone);

    assert(static_cast<int>(vBits.size()) == nLevel);
    vSymbol.set_size(0);

    std::vector<itpp::bvec> vecBits = vBits;

    std::vector<itpp::bvec> vecCode(nLevel);

    // 符号化していく
    for(int level = 0; level < nLevel; level++){
      // 確認用
      assert(vBits[level].size() <= static_cast<int>(vecLDPC[level].InfoLength()));
  
      // 情報ビットが足りない部分は0で補完する
      while(vecBits[level].size() < static_cast<int>(vecLDPC[level].InfoLength())){
        vecBits[level] = itpp::concat(vecBits[level], itpp::bin(0));
      }

      vecLDPC[level].Encode(vecBits[level], vecCode[level]);
    
    }
  
    for(int n = 0; n < static_cast<int>(nCodeLength); n++){
      itpp::bvec currentCodes(nLevel); 
      for(int level = 0; level < nLevel; level++){
        currentCodes[nLevel-level-1] = vecCode[level][n]; // itpp::Modulatorの都合上
                                // 低い要素番号ほど上位ビットを入れていく
                                // ここを逆順にするとLevel0がもっと復号結果よくなる
      }
      vSymbol = itpp::concat(vSymbol, Modulator.modulate_bits(currentCodes));
    }

  }
      
  // 各レベルのイテレーション回数を返す
  itpp::ivec MlcMsd::Decode(const itpp::cvec &vSymbol,
                            std::vector< itpp::bvec > &vDecodedBits,
                            const double      N0,
                            const unsigned    loopMax,
                            int        errConc)
  {
    assert(bSetDone);

    assert(vSymbol.size() == static_cast<int>(nCodeLength));
    itpp::ivec vLoop(nLevel);
    vLoop.zeros();

    if (errConc < 0){
      errConc = nLevel;
    } // if 
    
    vDecodedBits.resize(nLevel);
    for (int i = 0; i < nLevel; ++i){
      vDecodedBits[i].set_size(vecLDPC[i].InfoLength());
    } // for i
    
    // std::vector< itpp::bvec > vDecodedCodes(nLevel); // 復号された符号語
    // 情報ビットではない

    std::vector< itpp::Modulator_2D > vecMod(nCodeLength);
  
    for(int i = 0; i < static_cast<int>(nCodeLength); i++){
      vecMod[i] = Modulator;
    }
 
    itpp::bvec DecodedCodes;
    // At first, decoding level 0.
    vLoop[0]        = vecLDPC[0].Decode(vecMod, vSymbol, DecodedCodes, N0, loopMax);
    bool error = false;         // どこかのレベルでエラーが起きたらtrueにする
    // エラー隠蔽処理
    if (vLoop[0] == static_cast< int >(loopMax)){
      error = true;
      if (errConc == 0){
        for (int i = 0; i < nLevel; ++i){
          vLoop[i] = loopMax;
          vDecodedBits[i].zeros();
        } // for i
        return vLoop;
      } // if errConc
    } // if   vLoop[0]
    
    vDecodedBits[0] = DecodedCodes.left(vecLDPC[0].InfoLength());
    
    // Start decoding from level 0.
    for(int level = 1; level < nLevel; level++){
      // double tempN0 = N0 / static_cast<double>(vecLDPC.size()) * level; // ##

      // errorConc == 1のときはここで終わるかも
      if (level >= errConc && error){
        for (int k = level; k < nLevel; ++k){
          vLoop[k] = loopMax;
          vDecodedBits[k].zeros();
        } // for k
        return vLoop;
      } // if error
      
      // 全ての符号に対してModulatorをセットしていく
      for(int i = 0; i < static_cast<int>(nCodeLength); i++){
        itpp::cvec oldSymbols    = vecMod[i].get_symbols(); // ひとつ前のシンボル
        itpp::ivec oldBitmap     = vecMod[i].get_bits2symbols(); // ひとつ前のビットマップ
        itpp::cvec newSymbols(0);
        itpp::ivec newBitmap(0);
        int        code          = DecodedCodes[i]; // 前レベルの復号結果
        for(int s = 0; s < oldBitmap.size(); s++){
          if((oldBitmap[s] & 1) == code){
            newBitmap            = itpp::concat(newBitmap, oldBitmap[s] >> 1);
            newSymbols           = itpp::concat(newSymbols, oldSymbols[s]);
          }
        }
        vecMod[i].set(newSymbols, newBitmap);
      } // for i
      
      vLoop[level] = vecLDPC[level].Decode(vecMod, vSymbol, DecodedCodes, N0, loopMax);

      if (vLoop[level] == static_cast< int >(loopMax)){
        error = true;
      } // if 

      if ((level >= errConc) && error){
        for (int k = level; k < nLevel; ++k){
          vLoop[k] = loopMax;
          vDecodedBits[k].zeros();
        } // for k
        return vLoop;
      } // if error
      
      vDecodedBits[level] = DecodedCodes.left(vecLDPC[level].InfoLength());
    } // for level
       
    return vLoop;
  }
	 
  /******************* end of cMLC_MSD ************************/     

  double SetLDPC_MLC(std::vector< Ldpc > &vecLDPC,
                     unsigned nCodeLength,
                     itpp::ivec &vecRowWeight,
                     itpp::ivec &vecColWeight)
  {

    double sumCodeRate = 0.0;
    for(int i = 0; i < static_cast<int>(vecLDPC.size()); i++){
      vecLDPC[i].Set(nCodeLength, vecRowWeight[i], vecColWeight[i]);
      sumCodeRate += vecLDPC[i].CodeRate();
    }

    return sumCodeRate;
  
  }

  void SetModulatorNatural(itpp::Modulator_2D &Mod)
  {
    itpp::cvec symbols = Mod.get_symbols();
    itpp::ivec bitmap = Mod.get_bits2symbols();

    std::cout << "## default symbols = \n"
              << symbols << std::endl;

    std::cout << "## default bitmap = \n"
              << bitmap << std::endl;

    itpp::ivec newBitmap(bitmap.size());

    // set partitioning mapping
    for(int i = 0; i < newBitmap.size(); i++){
      newBitmap[i] = i;
    }

    std::cout << "## new bitmap = \n"
              << newBitmap << std::endl;

    Mod.set(symbols, newBitmap);
  }

  // 低いレイヤーの最後の部分を高いレイヤーから持ってきて埋める
  inline std::vector< itpp::bvec > AdjustInfoLength_forMLC(const std::vector< itpp::bvec > &infoBits,
                                                           const std::vector< int > &infoLengths_LDPC)
  {
    assert(infoBits.size() == infoLengths_LDPC.size());

    std::vector< itpp::bvec > output(infoBits.size());
    int startPoint = 0;
    for(int level = 0; level < static_cast<int>(infoBits.size()) - 1; level++)
      {
        output[level] = infoBits[level].get(startPoint, -1);
        startPoint = infoLengths_LDPC[level] - output[level].size();
        assert(startPoint >= 0);  // The info length must not be over LDPC info length.
        output[level] = itpp::concat(output[level], infoBits[level+1].mid(0, startPoint));
        assert(output[level].size() == infoLengths_LDPC[level]); // ## can delete this line.
      } // for level

    int level = infoBits.size() - 1;
    output[level] = infoBits[level].get(startPoint, -1);

    return output;
  }

  inline std::vector< itpp::bvec > AdjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo,
                                                            const std::vector< int > &infoLengths)
  {
    assert(decodedInfo.size() == infoLengths.size());

    std::vector< itpp::bvec > output(decodedInfo.size());
    output[0] = decodedInfo[0].mid(0, infoLengths[0]);
    int startPoint = infoLengths[0];
    for(int level = 1; level < static_cast<int>(decodedInfo.size()); level++)
      {
        output[level] = decodedInfo[level-1].get(startPoint, -1);
        startPoint = infoLengths[level] - output[level].size();
        assert(startPoint >= 0);   
        output[level] = itpp::concat(output[level], decodedInfo[level].mid(0, startPoint));
        assert(output[level].size() == infoLengths[level]); // ## can delete this line.
      
      }

    return output;
  }

}
