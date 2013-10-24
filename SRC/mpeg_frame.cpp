/*************************************************************/
// class CMpegFrame
// MP3�t���[���ɑ�������N���X
// �w�b�_�[�A�T�C�h���̃f�R�[�h���s���A�擾�����f�[�^��ێ�����
//
// ���C�t�^�C����1�t���[�������̊�
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

bool CMpegFrame::IsValidHeader() const
{
  if( m_nLayer>=0 && m_nLayer<=2
      &&  m_nBitrateIndex>=1 && m_nBitrateIndex<=14
      &&  m_nSamplingFrequency>=0 && m_nSamplingFrequency<=2
      &&  m_nEmphasis!=2 )
    {
      return true;
    }
  else{
    return false;
  }
}

int CMpegFrame::GetFrameSize() const
{
  int nSlots = (int)(144.0 * (double)GetBitrate() / (double)GetFrequency() );
  // static UINT16 bPrePadding = 2;	// ## Padding�̒l���ς�����Ƃ����o�͂��邽�߂Ɏg��
  //  static int nEachPadding = 0;	// ## padding���ƂɃt���[���T�C�Y�̘a���o�͂���

  if( m_nPaddingBit ){
    nSlots++;
  }

  /************* ## �����̓f�o�b�O�p ******************/
  
  //   if(bPrePadding != m_nPaddingBit){
  //     if(m_nPaddingBit){
  //       std::cout << "\n## Padding. Frame size = " << nSlots << "\n";
  //     }
  //     else{				// ## Debug
  //       std::cout << "\n## Not Padding. Frame size = " << nSlots << "\n";
  //     }
  //     bPrePadding = m_nPaddingBit;
  //   }
  
  // ## �ʃf�o�b�O
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

// �����ς���
int CMpegFrame::GetBitrate() const
{
  return g_nBitrate[m_nBitrateIndex];
}
