/*********************/
// cMpegToBinaryとcBinatyToMpegの定義
//
// Padding bitは伝送しない
// なぜならpaddingは受信側でも実装できるから
/***********************/
/******** to do list ***********/
// itppのbvecに変える
/*********************************/

#include <iostream>
#include <cassert>
#include <itpp/itbase.h>
#include "../include/mpeg_binary.h"


/************* class cMpegToBinary ******************/
// dataからnumビット分取り出してvBinベクトルに付け加える
// ただし、MSBを最初にする
inline void getBits(itpp::bvec &vInBin, unsigned short data, int num)
{
  assert(num > 0);

  for(int i = 1; i <= num; i++){
    itpp::bin tBin = (data >> (num-i)) & 1;
    vInBin = itpp::concat(vInBin,tBin);
  }
}

void cMpegToBinary::set(const CMpegFrame &tMpegFrame)
{
  MpegFrame = tMpegFrame;

  // binary初期化
  vSideBin.set_size(0);		
  vMainBin.set_size(0);
  
  getBits(vSideBin, MpegFrame.m_nMainDataBegin, 9); // ここはもしかしたら
					  // 省くかもしれない
  
  // Private Bitsも一応
  if(MpegFrame.GetChannels() == 1){
    getBits(vSideBin, 0, 5);	
  }
  else{
    getBits(vSideBin, 0, 3);
  }

  for(int ch=0; ch<MpegFrame.GetChannels(); ch++){
    for(int i=0; i<4; i++){
      getBits(vSideBin, MpegFrame.m_nScfSi[ch][i],1);
    }
  }


  for(int gr = 0; gr < 2; gr++){
    for(int ch = 0; ch < MpegFrame.GetChannels(); ch++){
      getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nPart23Length,12);
      getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nBigValues,9);
      getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nGlobalGain,8);
      getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nScalefacCompress,4);
      getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nWindowSwitchingFlag,1);
      if( MpegFrame.m_GranuleInfo[ch][gr].nWindowSwitchingFlag ){
	getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nBlockType,2);
	getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nMixedBlockFlag,1);
	
	// このケースでは、Region0とRegion1のみであるため
	// ハフマンテーブルは２つで問題ない。
	for( int window=0; window<2; window++ ){
	  getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nTableSelect[window],5);
	}

	for(int window=0; window<3; window++ ){
	  getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nSubblockGain[window],3);
	}
      }
      else { //if NOT m_bWindowSwitching
	// ノーマルロングブロック（ブロックタイプ0）
	for( int window=0; window<3; window++ ){
	  getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nTableSelect[window],5);
	}

	getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nRegion0Count,4);
	getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nRegion1Count,3);
      } 

      getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nPreFlag,1);
      getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nScalefacScale,1);
      getBits(vSideBin, MpegFrame.m_GranuleInfo[ch][gr].nCount1TableSelect,1);
    } //for ch
  } //for gr

  // メインデータの処理は一括
  int nMainDataSize = MpegFrame.m_nMainDataSize;
  // paddingの場合は1バイト分減らす
  if(MpegFrame.m_nPaddingBit){
    nMainDataSize -= 1;
  }

  for(int i = 0; i < nMainDataSize; i++){
    getBits(vMainBin, MpegFrame.m_pbyteMainData[i],8);
  }
  
  bSet = true;
}

void cMpegToBinary::usage()
{
  std::cout << "\nError in cMpegToBinary.\n"
	    << "Must set a MpegFrame.\n";
}

//
// バイナリデータを返す
itpp::bvec cMpegToBinary::getSideInfoBinary()
{
  if(!bSet){
    usage();
    exit(1);
  }

  return vSideBin;
}

itpp::bvec cMpegToBinary::getMainDataBinary()
{
  if(!bSet){
    usage();
    exit(1);
  }

  return vMainBin;
}

itpp::bvec cMpegToBinary::getFrameBinary()
{
  if(!bSet){
    usage();
    exit(1);
  }

  itpp::bvec vBin = itpp::concat(vSideBin, vMainBin);

  return vBin;
}
/********** end of cMpegToBinary *****************/

/********** class cBinaryToMpeg *******************/
// vBinのfrom番目からnumBitビットの整数を返す
inline unsigned short getInt(itpp::bvec &vBin,int &from, int numBit)
{
  assert((numBit > 0) && (from >= 0) );

  unsigned short out = 0;

  for(int i = 1; i <= numBit; i++){
    int temp = vBin[from];
    out += temp << (numBit-i);
    from++;
  }

  return out;
}

// フレーム全体のバイナリデータが入力のとき
void cBinaryToMpeg::set(const CMpegFrame &tMpegFrame, 
			const itpp::bvec &vInBin)
{
  // Side info とメインデータを切り離す
  int nSideInfoLength;
  if(tMpegFrame.GetChannels() == 1){
    nSideInfoLength = 17*8;
  }
  else{
    nSideInfoLength = 32*8;
  }

  itpp::bvec vInSideBin = vInBin.get(0,nSideInfoLength-1);
  itpp::bvec vInMainBin = vInBin.get(nSideInfoLength,-1);

  set(tMpegFrame,vInSideBin,vInMainBin);

}
  
// vInSideBinもvInMainBinも要素数多すぎても大丈夫
void cBinaryToMpeg::set(const CMpegFrame &tMpegFrame,
	 const itpp::bvec &vInSideBin,
	 const itpp::bvec &vInMainBin)
{
  // itpp::bvec vInBin = itpp::concat(vInSideBin, vInMainBin);

//   set(tMpegFrame, vInBin);

  MpegFrame = tMpegFrame;

  vSideBin = vInSideBin;
  vMainBin = vInMainBin;
  
  int nBitPos = 0;

  MpegFrame.m_nMainDataBegin = getInt(vSideBin,nBitPos,9); // ここはもしかしたら
					  // 省くかもしれない
  
  // Private Bitsはスルー
  if(MpegFrame.GetChannels()==1){
    getInt(vSideBin,nBitPos,5);
  }
  else{
    getInt(vSideBin,nBitPos,3);
  }

  
  for(int ch=0; ch<MpegFrame.GetChannels(); ch++){
    for(int i=0; i<4; i++){
      MpegFrame.m_nScfSi[ch][i] = getInt(vSideBin,nBitPos,1);
    }
  }

  for(int gr = 0; gr < 2; gr++){
    for(int ch = 0; ch < MpegFrame.GetChannels(); ch++){
      MpegFrame.m_GranuleInfo[ch][gr].nPart23Length = getInt(vSideBin,nBitPos,12);
      MpegFrame.m_GranuleInfo[ch][gr].nBigValues = getInt(vSideBin,nBitPos,9);
      MpegFrame.m_GranuleInfo[ch][gr].nGlobalGain = getInt(vSideBin,nBitPos,8);
      MpegFrame.m_GranuleInfo[ch][gr].nScalefacCompress = getInt(vSideBin,nBitPos,4);
      MpegFrame.m_GranuleInfo[ch][gr].nWindowSwitchingFlag = getInt(vSideBin,nBitPos,1);
      if( MpegFrame.m_GranuleInfo[ch][gr].nWindowSwitchingFlag ){
	MpegFrame.m_GranuleInfo[ch][gr].nBlockType = getInt(vSideBin,nBitPos,2);
	MpegFrame.m_GranuleInfo[ch][gr].nMixedBlockFlag = getInt(vSideBin,nBitPos,1);
	
	// このケースでは、Region0とRegion1のみであるため
	// ハフマンテーブルは２つで問題ない。
	for( int window=0; window<2; window++ ){
	  MpegFrame.m_GranuleInfo[ch][gr].nTableSelect[window] = getInt(vSideBin,nBitPos,5);
	}

	for(int window=0; window<3; window++ ){
	  MpegFrame.m_GranuleInfo[ch][gr].nSubblockGain[window] = getInt(vSideBin,nBitPos,3);
	}
      }
      else { //if NOT m_bWindowSwitching
	// ノーマルロングブロック（ブロックタイプ0）
	for( int window=0; window<3; window++ ){
	  MpegFrame.m_GranuleInfo[ch][gr].nTableSelect[window] = getInt(vSideBin,nBitPos,5);
	}

	MpegFrame.m_GranuleInfo[ch][gr].nRegion0Count = getInt(vSideBin,nBitPos,4);
	MpegFrame.m_GranuleInfo[ch][gr].nRegion1Count = getInt(vSideBin,nBitPos,3);
      } 

      MpegFrame.m_GranuleInfo[ch][gr].nPreFlag = getInt(vSideBin,nBitPos,1);
      MpegFrame.m_GranuleInfo[ch][gr].nScalefacScale = getInt(vSideBin,nBitPos,1);
      MpegFrame.m_GranuleInfo[ch][gr].nCount1TableSelect = getInt(vSideBin,nBitPos,1);
    } //for ch
  } //for gr
  

  // メインデータの処理は一括

  int nMainDataSize = MpegFrame.m_nMainDataSize;
  // paddingの場合は1バイト分減らす
  if(MpegFrame.m_nPaddingBit){
    nMainDataSize -= 1;
  }

  // 初期化
  nBitPos = 0;

  for(int i = 0; i < nMainDataSize; i++){
    MpegFrame.m_pbyteMainData[i] = getInt(vMainBin,nBitPos,8);
  }
  
  // ## Debug
  // assert(vBin.size() == nBitPos);

  bSet = true;


}

void cBinaryToMpeg::usage()
{
  std::cout << "\nError in cBinaryToMpeg.\n"
	    << "Must set a MpegFrame.\n";
}

  // tMpegFrameには元のCMpegFrameを入れる
CMpegFrame cBinaryToMpeg::getMpegFrame()
{

  if(!bSet){
    usage();
    exit(1);
  }

  //  std::cout << "## Debug\n";
 
  return MpegFrame;
}
  
/********************** end of cBinaryToMpeg **************/

/********************** class cMpegToArrangedBinary *********/
void cMpegToArrangedBinary::set(const std::vector<CMpegFrame> &tvecMpegFrame)
{
  vSideBin.set_size(0);
  vMainBin.set_size(0);

  
  for(int i = 0; i < static_cast<int>(tvecMpegFrame.size()); i++){
    itpp::bvec tempvSideBin,tempvMainBin;

    cMpegToBinary MpegToBinary(tvecMpegFrame[i]);

    tempvSideBin = MpegToBinary.getSideInfoBinary();
    tempvMainBin = MpegToBinary.getMainDataBinary();

    vSideBin = itpp::concat(vSideBin, tempvSideBin);
    vMainBin = itpp::concat(vMainBin, tempvMainBin);

  }

  bSet = true;
}

void cMpegToArrangedBinary::usage()
{
  std::cout << "\nError in cMpegToArrangedBinary.\n"
	    << "Must set a MpegFrame.\n";
}

itpp::bvec cMpegToArrangedBinary::getSideInfoBinary()
{
  if(!bSet){
    usage();
    exit(1);
  }

  return vSideBin;
}

itpp::bvec cMpegToArrangedBinary::getMainDataBinary()
{
  if(!bSet){
    usage();
    exit(1);
  }

  return vMainBin;
}

itpp::bvec cMpegToArrangedBinary::getBlockBinary()
{
  if(!bSet){
    usage();
    exit(1);
  }

  itpp::bvec vBin = itpp::concat(vSideBin,vMainBin);

  return vBin;
}

/******************* end of cMpegToArrangedBinary ************/

/******************* class cArrangedBinaryToMpeg ****************/
void cArrangedBinaryToMpeg::set(const std::vector<CMpegFrame> &tvecMpegFrame,
				const itpp::bvec &vInBin)
{
  itpp::bvec vInSideBin;
  itpp::bvec vInMainBin;

  unsigned nMainStartPos;
  if(tvecMpegFrame[0].GetChannels() == 1){
    nMainStartPos = 17*8*tvecMpegFrame.size();
  }
  else{
    nMainStartPos = 32*8*tvecMpegFrame.size();
  }
  vInSideBin = vInBin.get(0,nMainStartPos-1);
  vInMainBin = vInBin.get(nMainStartPos, -1); 


  set(tvecMpegFrame,vInSideBin,vInMainBin);

}

void cArrangedBinaryToMpeg::set(const std::vector<CMpegFrame> &tvecMpegFrame,
				const itpp::bvec &vInSideBin,
				const itpp::bvec &vInMainBin)
{
  vecMpegFrame = tvecMpegFrame;

  
  for(int i = 0; i < static_cast<int>(vecMpegFrame.size()); i++){
    // padding抜きのメインデータサイズを格納
    unsigned nMainDataSize = vecMpegFrame[i].GetMainDataSize();
    if(vecMpegFrame[i].GetPaddingBit()){
      nMainDataSize--;
    }
    nMainDataSize *= 8;		// バイトからビットへ変換
    
    // assert(nMainDataSize == vInMainBin.size()/vecMpegFrame.size());


    unsigned nSideInfoSize;
    if(vecMpegFrame[i].GetChannels() == 1){
      nSideInfoSize = 17*8;
    }
    else{
      nSideInfoSize = 32*8;
    }
    // assert(nSideInfoSize == vInSideBin.size()/vecMpegFrame.size());
    
    itpp::bvec tvSideBin = vInSideBin.get(i*nSideInfoSize,
					  (i+1)*nSideInfoSize-1);
    itpp::bvec tvMainBin = vInMainBin.get(i*nMainDataSize,
					   (i+1)*nMainDataSize-1); 


    //    std::cout << "## tvBin.size() = " << tvBin.size() << std::endl;

    cBinaryToMpeg BinaryToMpeg(vecMpegFrame[i],tvSideBin,tvMainBin);
    vecMpegFrame[i] = BinaryToMpeg.getMpegFrame();
    
  } // for i
  

  bSet = true;
}

void cArrangedBinaryToMpeg::usage()
{
  std::cout << "\nError in cArrangedBinaryToMpeg.\n"
	    << "Must set a MpegFrame.\n";
}


std::vector<CMpegFrame> cArrangedBinaryToMpeg::getVecMpegFrame()
{
  if(!bSet){
    usage();
    exit(1);
  }

  return vecMpegFrame;
}
