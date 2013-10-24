#include <cstring>
#include "../include/bit_stream.h"

CBitStream::CBitStream()
{
  Reset();
}

CBitStream::~CBitStream()
{
}

void CBitStream::Reset()
{
  m_nBufBitPos = 0;
  m_nBufBytePos = 0;

  m_nFramePos = 0;
  m_nSize = 0;
}

void CBitStream::SkipBits( int n )
{
  m_nBufBitPos += n;
  if( m_nBufBitPos>=8 ){
    m_nBufBytePos += m_nBufBitPos / 8;
    m_nBufBitPos %= 8;
  }
}

Uint32 CBitStream::PeekBits( int n )
{
  static bool bInit = false;
  static Uint32 wMask[32];

  if( n==0 ){
    return 0;
  }

  if( !bInit ){
    wMask[0] = 0;
    for( int i=1; i<32; i++ ){
      wMask[i] = wMask[i-1] | (1 << (32-i));
    }
    bInit = true;
  }

  Uint32 dwRet = (((Uint32)m_pbyte[m_nBufBytePos]<<24)
		  | ((Uint32)m_pbyte[m_nBufBytePos+1]<<16)
		  | ((Uint32)m_pbyte[m_nBufBytePos+2]<<8)
		  | (Uint32)m_pbyte[m_nBufBytePos+3]) & (wMask[n]>>m_nBufBitPos);

  dwRet >>= 32-m_nBufBitPos-n;

  return dwRet;
}

Uint32 CBitStream::GetBits( int n )
{
  if( n==0 ){
    return 0;
  }

  Uint32 wRet = PeekBits(n);
  SkipBits(n);

  return wRet;
}

bool CBitStream::IsOverrun()
{
  return m_nBufBytePos>=m_nSize;
}

bool CBitStream::SetData( Uint8* pData, int nSize )
{
  m_nBufBitPos = 0;
  m_nBufBytePos = 0;
  memmove(m_pbyte, pData, nSize);
  m_nSize = nSize;

  return true;
}

bool CBitStream::PutFrameMainData( const Uint8* pMainData, int nSize, int nMainDataBegin )
{
  // ビット列読み込み開始位置の設定
  m_nFramePos = m_nSize - nMainDataBegin;

  if( m_nSize>FIFO_SIZE-1000 ){
    // バッファのオーバフローが近づいたので
    // 読み取り済みの領域を削除しバッファを詰める
    memmove( m_pbyte, m_pbyte+m_nFramePos, nMainDataBegin );
    m_nFramePos = 0;
    m_nSize = nMainDataBegin;
  }

  // FIFOバッファにメインデータを追加する
  memmove( m_pbyte+m_nSize, pMainData, nSize );
  m_nSize += nSize;

  if( m_nFramePos<0 ) {
    // 初回フレームの場合、FALSEを返す
    return false;
  }

  m_nBufBytePos = m_nFramePos;
  m_nBufBitPos = 0;

  return true;
}

void CBitStream::ResetDataBitCount()
{
  m_nDataBitCount = m_nBufBytePos*8+m_nBufBitPos;
}

int CBitStream::GetDataBitCount() const
{
  return m_nBufBytePos*8+m_nBufBitPos-m_nDataBitCount;
}

bool CBitStream::SeekDataCountBitsPos( Uint32 nBits )
{
  Uint32 dwBitPos = m_nDataBitCount + nBits;
  m_nBufBytePos = dwBitPos / 8;
  m_nBufBitPos =  dwBitPos % 8;

  return true;
}
