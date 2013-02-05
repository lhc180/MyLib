#include <cassert>
#include "../include/mlc_msd.h"

 
/******************* class cMLC_MSD *************************/

cMLC_MSD::cMLC_MSD(std::vector<cLDPC> &vLDPC, itpp::Modulator_2D &Mod)
{
  set(vLDPC, Mod);
}

void cMLC_MSD::set(std::vector<cLDPC> &vLDPC, itpp::Modulator_2D &Mod)
{
  // vecLDPC.resize(vLDPC.size());
//   for(int i = 0; i < vLDPC.size(); i++){
//     &vecLDPC[i] = dynamic_cast< cLDPC_forMSD >(&vLDPC[i]);
// }
  vecLDPC = vLDPC;

  Modulator = Mod;
  
  assert(static_cast<int>(vecLDPC.size()) == Modulator.bits_per_symbol());

  nLevel = vecLDPC.size();
  
  nCodeLength = vecLDPC[0].getCodeLength();
  for(int level = 1; level < nLevel; level++){
    assert(nCodeLength == vecLDPC[level].getCodeLength());
  } 

  avrCodeRate = 0.0;
  for(int level = 0; level < nLevel; level++){
    avrCodeRate += vecLDPC[level].getRate();
  }

  bSetDone = true;

}

void cMLC_MSD::encode(const std::vector<itpp::bvec> &vBits,
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
    assert(vBits[level].size() <= static_cast<int>(vecLDPC[level].getInputLength()));
  
    // 情報ビットが足りない部分は0で補完する
    while(vecBits[level].size() < static_cast<int>(vecLDPC[level].getInputLength())){
      vecBits[level] = itpp::concat(vecBits[level], itpp::bin(0));
    }

    vecLDPC[level].encode(vecBits[level], vecCode[level]);
    
  }
  
  for(int n = 0; n < static_cast<int>(nCodeLength); n++){
    itpp::bvec currentCodes(nLevel); 
    for(int level = 0; level < nLevel; level++){
      currentCodes[nLevel-level-1] = vecCode[level][n]; // itpp::Modulatorの都合上
							// 低い要素番号ほど上位ビットを入れていく
    }
      vSymbol = itpp::concat(vSymbol, Modulator.modulate_bits(currentCodes));
  }

}
      
// 各レベルのイテレーション回数を返す
itpp::ivec cMLC_MSD::decode(const itpp::cvec &vSymbol,
			    std::vector< itpp::bvec > &vDecodedBits,
			    const double N0,
			    const unsigned loopMax)
{
  assert(bSetDone);

  assert(vSymbol.size() == static_cast<int>(nCodeLength));
  itpp::ivec vLoop(nLevel);


  vDecodedBits.resize(nLevel);
  // std::vector< itpp::bvec > vDecodedCodes(nLevel); // 復号された符号語
  // 情報ビットではない

  std::vector< itpp::Modulator_2D > vecMod(nCodeLength);
  
  for(int i = 0; i < static_cast<int>(nCodeLength); i++){
    vecMod[i] = Modulator;
  }
 
  itpp::bvec DecodedCodes;
  // At first, decoding level 0.
  vLoop[0] = vecLDPC[0].decode(vecMod, vSymbol, DecodedCodes, N0, loopMax);
  vDecodedBits[0] = DecodedCodes.left(vecLDPC[0].getInputLength());

  // Start decoding from level 0.
  for(int level = 1; level < nLevel; level++){
    // double tempN0 = N0 / static_cast<double>(vecLDPC.size()) * level; // ##
    
    // 全ての符号に対してModulatorをセットしていく
    for(int i = 0; i < static_cast<int>(nCodeLength); i++){
      itpp::cvec oldSymbols = vecMod[i].get_symbols(); // ひとつ前のシンボル
      itpp::ivec oldBitmap = vecMod[i].get_bits2symbols(); // ひとつ前のビットマップ
      itpp::cvec newSymbols(0);
      itpp::ivec newBitmap(0);
      int code = DecodedCodes[i];
      for(int s = 0; s < oldBitmap.size(); s++){
	if((oldBitmap[s] & 1) == code){
	  newBitmap = itpp::concat(newBitmap, oldBitmap[s] >> 1);
	  newSymbols = itpp::concat(newSymbols, oldSymbols[s]);
	}
      }
      vecMod[i].set(newSymbols, newBitmap);
    }
    
    vLoop[level] = vecLDPC[level].decode(vecMod, vSymbol, DecodedCodes, N0, loopMax);
    
    vDecodedBits[level] = DecodedCodes.left(vecLDPC[level].getInputLength());
       
     
  }
       
  return vLoop;
}
	 
/******************* end of cMLC_MSD ************************/     

double setLDPC_MLC(std::vector< cLDPC > &vecLDPC,
		   unsigned nCodeLength,
		   itpp::ivec &vecRowWeight,
		   itpp::ivec &vecColWeight)
{

  double sumCodeRate = 0.0;
  for(int i = 0; i < static_cast<int>(vecLDPC.size()); i++){
    vecLDPC[i].setEncoder(nCodeLength, vecRowWeight[i], vecColWeight[i]);
    sumCodeRate += vecLDPC[i].getRate();
  }

  return sumCodeRate;
  
}

void setModulatorNatural(itpp::Modulator_2D &Mod)
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
std::vector< itpp::bvec > adjustInfoLength_forMLC(const std::vector< itpp::bvec > &infoBits,
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

std::vector< itpp::bvec > adjustInfoLength_fromMSD(const std::vector< itpp::bvec > &decodedInfo,
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

