/****************************************************************/
// class CMpegFile
// MP3ファイルに相当するクラス
// MP3ファイルにアクセスし、フレームを順次に読み込む
//
// ライフタイムはMP3に対するアクセスの間
/****************************************************************/
#include <cstdio>
#include <cstring>
#include <iostream>
#include "../include/mpeg.h"

CMpegFile::CMpegFile()
{
  m_pFile = NULL;
}

CMpegFile::~CMpegFile()
{
  if( IsOpen() ){
    Close();
  }
}

bool CMpegFile::Open( const char* pFile )
{
  m_pFile = fopen( pFile, "rb" );
  if( m_pFile==NULL ){
    return false;
  }

  return true;
}

void CMpegFile::Close()
{
  if( m_pFile!=NULL ){
    fclose(m_pFile);
    m_pFile = NULL;
  }
}

inline unsigned short CalcCRC( unsigned char* pCheckBuf, int nBytes )
{
  const unsigned int CRC16_POLYNOMIAL = 0x18005;

  unsigned int crc = 0xffff;
  for( int i=0; i<nBytes; i++ ){
    unsigned short byte = pCheckBuf[i];
    for( int i=0; i<8; i++ ){
      crc <<= 1;
      byte <<= 1;
      int nCRCCarry = (crc & 0x10000)!=0;
      int nDataCarry = (byte & 0x100)!=0;
      if( (!nCRCCarry && nDataCarry) || (nCRCCarry && !nDataCarry)){
	crc ^= CRC16_POLYNOMIAL;
      }

      crc &= 0xffff;
    }
  }

  return (unsigned short)crc;
}

bool CMpegFile::SeekMpegFrame( CMpegFrame* pMpegFrame )
{
  Uint8 byte = 0;
  Uint8 byte2[2];
  Uint8 byteHeader2[2];
  Uint32 dwPos;
  Uint16 word;
        
  //////////////////////////////////////////////////////////////////
  // フレームヘッダーの検索
  //////////////////////////////////////////////////////////////////
  
  for(;;){
    dwPos = ftell(m_pFile);
    if( fread( &byte2, 2, 1, m_pFile )==0 ){
      // std::cout << "## Error in serching sync words.\n";
      return false;
    }
                
    word = (byte2[0]<<8) + byte2[1];
    for( int i=0; (word & 0xfff0)!=0xfff0; i++ ){
      if( fread( &byte, 1, 1, m_pFile )==0 ){
	// std::cout << "## Error in serching sync words.\n";
	return false;
      }
                        
      word <<= 8;
      word |= byte;
                        
      // 4096バイトの中に同期ヘッダーが見つからない場合、
      // MP3ファイルではないと判断する
      if( i>4096*2 ){
	std::cerr << "## Error: Can not find a sync header.\n";
	return false;
      }
    }
                
    pMpegFrame->m_nID = (word & 0x0008) >> 3;
    pMpegFrame->m_nLayer = (word & 0x0006) >> 1;
    pMpegFrame->m_nProtectionBit = word & 0x0001; 
                
    if( fread( &byte2, 2, 1, m_pFile)==0 ){
      std::cout << "## Error: Illegal header file.\n";
      return false;
    }
    memcpy( byteHeader2, byte2, 2 );
    word = (byte2[0]<<8) + byte2[1];
                
    pMpegFrame->m_nBitrateIndex = (word & 0xf000) >> 12;
    pMpegFrame->m_nSamplingFrequency = (word & 0x0c00) >> 10;
    pMpegFrame->m_nPaddingBit = (word & 0x0200) >> 9;
    pMpegFrame->m_nMode = (word & 0x00c0) >> 6;
    pMpegFrame->m_nModeExtention = (word & 0x0030) >> 4;
    pMpegFrame->m_nCopyright = (word & 0x0008) >> 3;
    pMpegFrame->m_nOriginal = (word & 0x0004) >> 2;
    pMpegFrame->m_nEmphasis = word & 0x0003;
                
    if( pMpegFrame->IsValidHeader() ){
      break;
    }
                
    fseek(m_pFile, -3, SEEK_CUR );
  }

  //////////////////////////////////////////////////////////////////
  // CRCワードの読み込み
  //////////////////////////////////////////////////////////////////
  if( pMpegFrame->m_nProtectionBit==0 ){
    if( fread( byte2, 2, 1, m_pFile)==0 ){
      return false;
    }
    pMpegFrame->m_nCRCCheck = (byte2[0]<<8) + byte2[1];
  }
  else{
    pMpegFrame->m_nCRCCheck = 0;
  }
  
  //////////////////////////////////////////////////////////////////
  // サイド情報の読み込み
  //////////////////////////////////////////////////////////////////
  Uint8 byteSide[64];
  int nSideSize;
  if( pMpegFrame->GetChannels()==1 ){
    nSideSize = 17;
  }
  else{
    nSideSize = 32;
  }

  if( fread( byteSide, nSideSize, 1, m_pFile)==0 ){
    return false;
  }
  CBitStream bs( byteSide, nSideSize );
        
  if( !DecodeSideInfo( &bs, pMpegFrame ) ){
    return false;
  }
  
  //////////////////////////////////////////////////////////////////
  // CRCチェック
  //////////////////////////////////////////////////////////////////
  if( !pMpegFrame->m_nProtectionBit ){
    Uint8 pCheckBuf[34];
    memcpy( pCheckBuf, byteHeader2, 2 ); //ヘッダー
    memcpy( pCheckBuf+2, byteSide, nSideSize ); //サイド情報
    Uint16 crc = CalcCRC( pCheckBuf, nSideSize+2 );
    if( crc!=pMpegFrame->m_nCRCCheck ){
      //CRCエラーの場合、デコードを中止する
      return false;
      //continue;
    }
  }
 
  //////////////////////////////////////////////////////////////////
  // サポートフォーマットのチェック
  //////////////////////////////////////////////////////////////////
  if( pMpegFrame->m_nID!=1 || pMpegFrame->m_nLayer!=1 || pMpegFrame->m_nBitrateIndex==0 ){
    return false;
  }
 
  //////////////////////////////////////////////////////////////////
  // メインデータの読み込み
  //////////////////////////////////////////////////////////////////
  int nMainDataSize = pMpegFrame->GetMainDataSize();
  if( nMainDataSize<0 ){
    return false;
  }
  // std::cout << "## nMainSize = " << nMainSize << std::endl;
  
  if( fread( pMpegFrame->m_pbyteMainData, nMainDataSize, 1, m_pFile)==0 ){
    // std::cout << "Error: Can not get main data.\n";
    return false;
  }
  pMpegFrame->m_nMainDataSize = nMainDataSize;
  // std::cout << "## Serching header.\n";  
  return true;
}

bool CMpegFile::DecodeSideInfo( CBitStream* pbs, CMpegFrame* pMpegFrame )
{
  int gr, ch, i, window;
  pMpegFrame->m_nMainDataBegin = (Uint16)pbs->GetBits(9);
  if( pMpegFrame->GetChannels()==1 ){
    //Private Bits
    pbs->SkipBits(5);
  }
  else{
    //Private Bits
    pbs->SkipBits(3);
  }

  for( ch=0; ch<pMpegFrame->GetChannels(); ch++){
    for( i=0; i<4; i++){
      pMpegFrame->m_nScfSi[ch][i] = (Uint16)pbs->GetBits(1);
    }
  }

  for( gr=0; gr<2; gr++){
    for( ch=0; ch<pMpegFrame->GetChannels(); ch++){
      pMpegFrame->m_GranuleInfo[ch][gr].nPart23Length = (Uint16)pbs->GetBits(12);
      pMpegFrame->m_GranuleInfo[ch][gr].nBigValues = (Uint16)pbs->GetBits(9);
      pMpegFrame->m_GranuleInfo[ch][gr].nGlobalGain = (Uint16)pbs->GetBits(8);
      pMpegFrame->m_GranuleInfo[ch][gr].nScalefacCompress = (Uint16)pbs->GetBits(4);
      pMpegFrame->m_GranuleInfo[ch][gr].nWindowSwitchingFlag = pbs->GetBits(1);
      if( pMpegFrame->m_GranuleInfo[ch][gr].nWindowSwitchingFlag ){
	pMpegFrame->m_GranuleInfo[ch][gr].nBlockType = (Uint16)pbs->GetBits(2);
	pMpegFrame->m_GranuleInfo[ch][gr].nMixedBlockFlag = pbs->GetBits(1);
	if (pMpegFrame->m_GranuleInfo[ch][gr].nBlockType==0){
	  // このケースで、ブロックタイプ0は不正
	  return false;
	}

	// このケースでは、Region0とRegion1のみであるため
	// ハフマンテーブルは２つで問題ない。
	for( window=0; window<2; window++ ){
	  pMpegFrame->m_GranuleInfo[ch][gr].nTableSelect[window] = (Uint16)pbs->GetBits(5);
	}

	for( window=0; window<3; window++ ){
	  pMpegFrame->m_GranuleInfo[ch][gr].nSubblockGain[window] = (Uint16)pbs->GetBits(3);
	}
      }
      else { //if NOT m_bWindowSwitching
	// ノーマルロングブロック（ブロックタイプ0）
	for( int window=0; window<3; window++ ){
	  pMpegFrame->m_GranuleInfo[ch][gr].nTableSelect[window] = (Uint16)pbs->GetBits(5);
	}

	pMpegFrame->m_GranuleInfo[ch][gr].nRegion0Count = (Uint16)pbs->GetBits(4);
	pMpegFrame->m_GranuleInfo[ch][gr].nRegion1Count = (Uint16)pbs->GetBits(3);
	pMpegFrame->m_GranuleInfo[ch][gr].nMixedBlockFlag = false;
	pMpegFrame->m_GranuleInfo[ch][gr].nBlockType = 0;
      } 

      pMpegFrame->m_GranuleInfo[ch][gr].nPreFlag = (Uint16)pbs->GetBits(1);
      pMpegFrame->m_GranuleInfo[ch][gr].nScalefacScale = (Uint16)pbs->GetBits(1);
      pMpegFrame->m_GranuleInfo[ch][gr].nCount1TableSelect = (Uint16)pbs->GetBits(1);
    } //for ch
  } //for gr

  return true;
}
