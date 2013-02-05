/*************************************************************/
// class CMpegFrame
// MP3フレームに相当するクラス
// ヘッダー、サイド情報のデコードを行い、取得したデータを保持する
//
// ライフタイムは1フレーム処理の間
/**************************************************************/

#include <iostream>
#include <cstdio>
#include <cstring>
#include "../include/mpeg.h"

CMpegFrame::CMpegFrame()
{
  m_nCRCCheck = 0;
  m_nMainDataSize = 0;
  memset( m_GranuleInfo, 0, sizeof(m_GranuleInfo) );
  memset( m_nScfSi, 0, sizeof(m_nScfSi) );
}

CMpegFrame::~CMpegFrame()
{
}

BOOL CMpegFrame::IsValidHeader() const
{
  if( m_nLayer>=0 && m_nLayer<=2
      &&  m_nBitrateIndex>=1 && m_nBitrateIndex<=14
      &&  m_nSamplingFrequency>=0 && m_nSamplingFrequency<=2
      &&  m_nEmphasis!=2 )
    {
      return TRUE;
    }
  else{
    return FALSE;
  }
}

int CMpegFrame::GetFrameSize() const
{
  int nSlots = (int)(144.0 * (double)GetBitrate() / (double)GetFrequency() );
  // static UINT16 bPrePadding = 2;	// ## Paddingの値が変わったときを出力するために使う
  //  static int nEachPadding = 0;	// ## paddingごとにフレームサイズの和を出力する

  if( m_nPaddingBit ){
    nSlots++;
  }

  /************* ## ここはデバッグ用 ******************/
  
  //   if(bPrePadding != m_nPaddingBit){
  //     if(m_nPaddingBit){
  //       std::cout << "\n## Padding. Frame size = " << nSlots << "\n";
  //     }
  //     else{				// ## Debug
  //       std::cout << "\n## Not Padding. Frame size = " << nSlots << "\n";
  //     }
  //     bPrePadding = m_nPaddingBit;
  //   }
  
  // ## 別デバッグ
//   nEachPadding += nSlots;
//   if(!m_nPaddingBit){
//     std::cout << "\n## nEachPadding = " << nEachPadding << "\n";
//     nEachPadding = 0;
//   }


  /************************************/
  return nSlots;
}

int CMpegFrame::GetMainDataSize() const
{
  int nSlots = GetFrameSize();

  //Header Bytes
  nSlots -= 4;
        
  //CRC Bytes
  if( m_nProtectionBit==0 ){
    nSlots -= 2; 
  }

  //SideInfo Bytes
  if( GetChannels()==1 ){
    nSlots -= 17;
  }
  else{
    nSlots -=32;
  }

  return nSlots;
}

int CMpegFrame::GetSampleNum() const
{
  return 1152;
}

int g_nFrequency[] = { 44100, 48000, 32000, 0 };
int g_nBitrate[] =
  { 0, 32000, 40000, 48000, 56000, 64000, 80000, 96000, 112000, 128000, 160000, 192000, 224000, 256000, 320000 };


int CMpegFrame::GetFrequency() const
{
  return g_nFrequency[m_nSamplingFrequency];
}

// 少し変えた
int CMpegFrame::GetBitrate() const
{
  return g_nBitrate[m_nBitrateIndex];
}
