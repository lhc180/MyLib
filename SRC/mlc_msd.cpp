// #include <boost/thread.hpp>
// #include <boost/ref.hpp>
#include <cassert>
#include "../include/mlc_msd.h"

/******************* class MlcMsd *************************/
namespace mylib{
  MlcMsd::MlcMsd(std::vector<Ldpc> &vLDPC, itpp::Modulator_2D &Mod)
  {
    Set(vLDPC, Mod);
  }
  
  void MlcMsd::Encode(const std::vector<itpp::bvec> &vBits,
                      itpp::cvec &vSymbol)
  {
    assert(setDone_);

    assert(static_cast<int>(vBits.size()) == numLevels_);
    // vSymbol.set_size(nCodeLength);

    std::vector<itpp::bvec> vecBits = vBits;

    std::vector<itpp::bvec> vecCode(numLevels_);

    // 符号化していく
    for(int level = 0; level < numLevels_; level++){
      // 確認用
      assert(vBits[level].size() <= static_cast<int>(vecLDPC_[level].InfoLength()));
  
      // 情報ビットが足りない部分は0で補完する
      while(vecBits[level].size() < static_cast<int>(vecLDPC_[level].InfoLength())){
        vecBits[level] = itpp::concat(vecBits[level], itpp::bin(0));
      }

      vecLDPC_[level].Encode(vecBits[level], vecCode[level]);
    
    } // for level

    itpp::bvec codes(numLevels_*codeLength_); 
    for(int n = 0; n < static_cast<int>(codeLength_); n++){
      for(int level = 0; level < numLevels_; level++){
        codes[n*numLevels_ + level] = vecCode[level][n]; // itpp::Modulatorの都合上
        // 低い要素番号ほど上位ビットを入れていく
        // つまりlevel 0はMSB
      }
    } // for n
    vSymbol = modulator_.modulate_bits(codes);
  }

  void MlcMsd::EncodeTest(const std::vector<itpp::bvec> &vBits,
                          itpp::cvec &symbols,
                          std::vector<itpp::bvec> &vecCode)
  {
    assert(setDone_);

    assert(static_cast<int>(vBits.size()) == numLevels_);
    // vSymbol.set_size(nCodeLength);

    std::vector<itpp::bvec> vecBits = vBits;

    vecCode.resize(numLevels_);

    // 符号化していく
    for(int level = 0; level < numLevels_; level++){
      // 確認用
      assert(vBits[level].size() <= static_cast<int>(vecLDPC_[level].InfoLength()));
  
      // 情報ビットが足りない部分は0で補完する
      while(vecBits[level].size() < static_cast<int>(vecLDPC_[level].InfoLength())){
        vecBits[level] = itpp::concat(vecBits[level], itpp::bin(0));
      }

      vecLDPC_[level].Encode(vecBits[level], vecCode[level]);
    
    } // for level

    itpp::bvec codes(numLevels_*codeLength_); 
    for(int n = 0; n < static_cast<int>(codeLength_); n++){
      for(int level = 0; level < numLevels_; level++){
        codes[n*numLevels_ + level] = vecCode[level][n]; // itpp::Modulatorの都合上
        // 低い要素番号ほど上位ビットを入れていく
        // つまりlevel 0はMSB
      }
    } // for n
    symbols = modulator_.modulate_bits(codes);
  }
  
  void MlcMsd::LevelDownModulator(std::vector< itpp::Modulator_2D >* vecMod, const itpp::bvec& decodedCode)
  {
    // 全ての符号に対してModulatorをセットしていく
    for(int i = 0; i < static_cast<int>(codeLength_); i++){
      itpp::cvec symbols      = (*vecMod)[i].get_symbols(); // ひとつ前のシンボル
      itpp::ivec bits2symbols = (*vecMod)[i].get_bits2symbols(); // ひとつ前のビットマップ
      int        mapNum       = bits2symbols.size();
          
      itpp::cvec newSymbols(mapNum/2);
      itpp::ivec newBits2Symbols(mapNum/2);
      itpp::bin  code = decodedCode[i]; // 前レベルの復号結果
      
      int startNum = 0;
      if (code == 1){
        startNum   = mapNum/2;
      }                         // if code
          
      for (int s = 0; s < mapNum/2; ++s){ // s+startNumが前の変調器のマッピング
        newBits2Symbols[s] = s;           // sが新しい変調器の変調前の入力ビット列と仮定
        newSymbols[s]      = symbols[bits2symbols[s+startNum]];
      }                         // for s

      (*vecMod)[i].set(newSymbols, newBits2Symbols);
    }                           // for i
      
  }
  
  // 各レベルのイテレーション回数を返す
  itpp::ivec MlcMsd::Decode(const itpp::cvec &vSymbol,
                            std::vector< itpp::bvec > &vDecodedBits,
                            const double      N0,
                            const unsigned    loopMax,
                            int        errConc)
  {
    assert(setDone_);

    assert(vSymbol.size() == static_cast<int>(codeLength_));
    itpp::ivec numLoops(numLevels_);
    numLoops.zeros();

    if (errConc < 0){
      errConc = numLevels_;
    } // if 
    
    vDecodedBits.resize(numLevels_);
    for (int i = 0; i < numLevels_; ++i){
      vDecodedBits[i].set_size(vecLDPC_[i].InfoLength());
    } // for i
    
      // std::vector< itpp::bvec > vDecodedCodes(nLevel); // 復号された符号語
      // 情報ビットではない

    std::vector< itpp::Modulator_2D > vecMod(codeLength_, modulator_);
    
    itpp::bvec DecodedCodes;
    // At first, decoding level 0.
    numLoops[0]        = vecLDPC_[0].Decode(vecMod, vSymbol, DecodedCodes, N0, loopMax);
    bool error = false;         // どこかのレベルでエラーが起きたらtrueにする
    // エラー隠蔽処理
    if (numLoops[0] == static_cast< int >(loopMax)){
      error = true;
      if (errConc == 0){
        for (int i = 0; i < numLevels_; ++i){
          numLoops[i] = loopMax;
          vDecodedBits[i].zeros();
        } // for i
        return numLoops;
      } // if errConc
    } // if   vLoop[0]
    
    vDecodedBits[0] = DecodedCodes.left(vecLDPC_[0].InfoLength());
    
    // Start decoding from level 0.
    for(int level = 1; level < numLevels_; level++){
      // double tempN0 = N0 / static_cast<double>(vecLDPC.size()) * level; // ##

      // errorConc == 1のときはここで終わるかも
      if (level >= errConc && error){
        for (int k = level; k < numLevels_; ++k){
          numLoops[k] = loopMax;
          vDecodedBits[k].zeros();
        } // for k
        return numLoops;
      } // if error

      LevelDownModulator(&vecMod, DecodedCodes);
      
      numLoops[level] = vecLDPC_[level].Decode(vecMod, vSymbol, DecodedCodes, N0, loopMax);

      if (numLoops[level] == static_cast< int >(loopMax)){
        error = true;
      } // if 

      if ((level >= errConc) && error){
        for (int k = level; k < numLevels_; ++k){
          numLoops[k] = loopMax;
          vDecodedBits[k].zeros();
        } // for k
        return numLoops;
      } // if error
      
      vDecodedBits[level] = DecodedCodes.left(vecLDPC_[level].InfoLength());
    } // for level
       
    return numLoops;
  }

  // ## デバッグ用
  itpp::bvec MlcMsd::DecodeTest(const itpp::cvec& symbols,
                        const std::vector< itpp::bvec >& knownBits,
                        double n0,
                        unsigned loopMax)
  {
    assert(static_cast< int >(knownBits.size()) < numLevels_);
    assert(setDone_);
    assert(symbols.size() == static_cast<int>(codeLength_));

    std::vector< itpp::Modulator_2D > vecMod(codeLength_, modulator_);
    
    for (int i = 0; i < static_cast< int >(knownBits.size()); ++i){
      LevelDownModulator(&vecMod, knownBits[i]);
    } // for i

    itpp::bvec decodedCode;
    vecLDPC_[knownBits.size()].Decode(vecMod, symbols, decodedCode, n0, loopMax);

    // ##
    // std::cout << std::endl;
    // std::cout << "## decodedCode.size() == " << decodedCode.size() << std::endl;
    // std::cout << "## loop = " << loop << std::endl;
    // std::cout << decodedCode << std::endl;
    
    return decodedCode.left(vecLDPC_[knownBits.size()].InfoLength());
  }
  
  /******************* end of cMLC_MSD ************************/     

  // ++++++++++++++++++++ MlcMsdWith8pskBp ++++++++++++++++++++

  std::vector< itpp::bvec > MlcMsdWith8pskBp::Decode(const itpp::cvec &vSymbol,
                                                     double N0,
                                                     unsigned loopMax,
                                                     int errConc)
  {
    //std::cout << "## MlcMsdWith8pskBp::Decode was called." << std::endl;
    assert(setDone_);
    
    assert(vSymbol.size() == static_cast< int >(codeLength_));
    itpp::ivec numLoops(numLevels_);
    numLoops.zeros();

    if (errConc < 0){
      errConc = numLevels_;
    } // if
    
    std::vector< itpp::bvec > vDecodedBits(numLevels_);
    for (int i = 0; i < numLevels_; ++i){
      vDecodedBits[i].set_size(vecLDPC_[i].InfoLength());
    } // for i
    
    std::vector< itpp::Modulator_2D > vecMod(codeLength_, modulator_);

    /*
    // ここから並列処理

    DecoderParas decoded0, decoded1;
    
    
    // int loop0, loop1;           // 戻り値をもらうマルチスレッドはめんどくさそうなので，引数でloop0, loop1をもらってくる
    boost::packaged_task< DecoderParas > pt0(boost::bind(&Ldpc::DecodeAtLevel, &vecLDPC_[0], modulator_, vSymbol, N0, 0, loopMax));
    boost::unique_future< DecoderParas > uf0 = pt0.get_future();
    
    boost::packaged_task< DecoderParas > pt1(boost::bind(&Ldpc::DecodeAtLevel, &vecLDPC_[1], modulator_, vSymbol, N0, 1, loopMax));
    boost::unique_future< DecoderParas > uf1 = pt1.get_future();

    boost::thread thread0(boost::ref(pt0));
    boost::thread thread1(boost::ref(pt1));

    // thread0.join();
    // thread1.join();
    
    decoded0 = uf0.get();
    decoded1 = uf1.get();
    
    numLoops[0] = decoded0.loop;
    numLoops[1] = decoded1.loop;
    */

    itpp::bvec decoded0, decoded1;
    numLoops[0] = vecLDPC_[0].Decode(vecMod, vSymbol, decoded0, N0, loopMax);
    // std::cout << "## numLoops[0] = " << numLoops[0] << std::endl;
    // std::cout << "## numLoops[1] = " << numLoops[1] << std::endl;
    // std::cout << "## loopMax = " << static_cast< int >(loopMax) << std::endl;
    
    // std::cout << "## decoded0.decodedBits.size() = " << decoded0.decodedBits.size() << std::endl;
    // std::cout << "## decoded1.decodedBits.size() = " << decoded1.decodedBits.size() << std::endl;

    bool error = false;         // どこかのでエラーが起きたらtrueにする
    
    if (numLoops[0] == static_cast< int >(loopMax) ){
      error = true;
      if (errConc <= 1){
        vDecodedBits[0].zeros();
      } // if
      else{
        vDecodedBits[0] = decoded0.left(vecLDPC_[0].InfoLength());
      } // else 
    } // if
    else{
      vDecodedBits[0] = decoded0.left(vecLDPC_[0].InfoLength());
    } // else 

    MlcMsd::LevelDownModulator(&vecMod, decoded0);
    numLoops[1] = vecLDPC_[1].Decode(vecMod, vSymbol, decoded1, N0, loopMax);

    
    if (numLoops[1] == static_cast< int >(loopMax) ){
      error = true;
      if (errConc <= 1){
        vDecodedBits[1].zeros();
      } // if
      else{
        vDecodedBits[1] = decoded1.left(vecLDPC_[1].InfoLength());
      } // else 
    } // if
    else{
      vDecodedBits[1] = decoded1.left(vecLDPC_[1].InfoLength());
    } // else 

    if (error && (errConc <= 1)){
      vDecodedBits[2].zeros();
      numLoops[2] = loopMax;
      return vDecodedBits;
    } // if

    MlcMsd::LevelDownModulator(&vecMod, decoded1);
    
    itpp::bvec decodedCode;
    numLoops[2] = vecLDPC_[2].Decode(vecMod, vSymbol, decodedCode, N0, loopMax);

    if (numLoops[2] == static_cast< int >(loopMax)){
      error = true;
    } // if 
    
    if (error && (errConc <= 2)){
      vDecodedBits[2].zeros();
      return vDecodedBits;
    } // if


    vDecodedBits[2] = decodedCode.left(vecLDPC_[2].InfoLength());
    // std::cout << "## vDecodedBits[2] = " << vDecodedBits[2] << std::endl;
    
    // std::cout << "## vDecodedBits[1] = " << vDecodedBits[1] << std::endl;
    // std::cout << "## vDecodedBits[2] = " << vDecodedBits[2] << std::endl;
    
    
    return vDecodedBits;
  }
  
  // -------------------- MlcMsdWith8pskBp --------------------
  
  /****************** class NaturalPSK ************************/
  void SetPartitioningPSK::Init()
  {
    itpp::ivec newBitmap(M);
    
    for (int i = 0; i < M; ++i){
      // std::cout << "## i = " << i << std::endl;
      // std::cout << "## i_reversed = " << i_reversed << std::endl;      
      newBitmap[i] = itpp::reverse_int(k, i);
    } // for i

    itpp::cvec newSymbols = symbols;

    set(newSymbols, newBitmap);
  }
  /******************* end of NaturalPSK **************************/

  void ReverseSpPSK::Init()
  {
    itpp::ivec newBitmap(M);

    for (int i = 0; i < M; ++i){
      newBitmap[i] = i;
    } // for i

    itpp::cvec newSymbols = symbols;

    set(newSymbols, newBitmap);
  }
  
  
  double SetLdpcForMlc(std::vector< Ldpc > &vecLDPC_,
                     unsigned codeLength_,
                     itpp::ivec &vecRowWeight,
                     itpp::ivec &vecColWeight)
  {
    std::cout << "## codeLength = " << codeLength_ << std::endl;
    
    double sumCodeRate = 0.0;
    for(int i = 0; i < static_cast<int>(vecLDPC_.size()); i++){
      vecLDPC_[i].Set(codeLength_, vecRowWeight[i], vecColWeight[i]);
      sumCodeRate += vecLDPC_[i].CodeRate();
    }

    return sumCodeRate;
  
  }

  // 低いレイヤーの最後の部分を高いレイヤーから持ってきて埋める
  // inline std::vector< itpp::bvec > AdjustInfoLength_forMLC(const std::vector< itpp::bvec > &infoBits,
  //                                                          const std::vector< int > &infoLengths_LDPC)
  // {
  //   assert(infoBits.size() == infoLengths_LDPC.size());

  //   std::vector< itpp::bvec > output(infoBits.size());
  //   int startPoint = 0;
  //   for(int level = 0; level < static_cast<int>(infoBits.size()) - 1; level++)
  //     {
  //       output[level] = infoBits[level].get(startPoint, -1);
  //       startPoint = infoLengths_LDPC[level] - output[level].size();
  //       assert(startPoint >= 0);  // The info length must not be over LDPC info length.
  //       output[level] = itpp::concat(output[level], infoBits[level+1].mid(0, startPoint));
  //       assert(output[level].size() == infoLengths_LDPC[level]); // ## can delete this line.
  //     } // for level

  //   int level = infoBits.size() - 1;
  //   output[level] = infoBits[level].get(startPoint, -1);

  //   return output;
  // }

  // inline std::vector< itpp::bvec > AdjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo,
  //                                                           const std::vector< int > &infoLengths)
  // {
  //   assert(decodedInfo.size() == infoLengths.size());

  //   std::vector< itpp::bvec > output(decodedInfo.size());
  //   output[0] = decodedInfo[0].mid(0, infoLengths[0]);
  //   int startPoint = infoLengths[0];
  //   for(int level = 1; level < static_cast<int>(decodedInfo.size()); level++)
  //     {
  //       output[level] = decodedInfo[level-1].get(startPoint, -1);
  //       startPoint = infoLengths[level] - output[level].size();
  //       assert(startPoint >= 0);   
  //       output[level] = itpp::concat(output[level], decodedInfo[level].mid(0, startPoint));
  //       assert(output[level].size() == infoLengths[level]); // ## can delete this line.
      
  //     }

  //   return output;
  // }

}
