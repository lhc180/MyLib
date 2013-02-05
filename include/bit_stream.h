#ifndef BIT_STREAM_H
#define BIT_STREAM_H

/***********************************************************/
// MP3のビットストリーム用のクラス
/***********************************************************/

// myimg.hの名前衝突と回避する
#ifndef TYPEDEF_INT32
#define TYPEDEF_INT32
typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef int INT32;
#endif

typedef char INT8;
typedef unsigned int UINT32;

#ifndef TRUE
#define TRUE    true
#endif
#ifndef FALSE
#define FALSE   false
#endif
#ifndef BOOL
#define BOOL bool
#endif

#define FIFO_SIZE 6144

class CBitStream
{
private:
  UINT8   m_pbyte[FIFO_SIZE];
  int     m_nFramePos;
  int     m_nBufBitPos;
  int     m_nBufBytePos;
  int     m_nSize;
  int     m_nDataBitCount;

public:
  CBitStream();
  CBitStream( UINT8* pData, int nSize ) { CBitStream(); SetData( pData, nSize ); }
  virtual ~CBitStream();
  void Reset();

  void SkipBits( int n );
  UINT32 PeekBits( int n );
  UINT32 GetBits( int n );
  UINT16 Get1Bit() { return (UINT16)GetBits(1); }
  void PurgeBits( int nBits ) { GetBits(nBits); }
  BOOL IsOverrun();
  BOOL SetData( UINT8* pData, int nSize );

  BOOL PutFrameMainData( const UINT8* pMainData, int nSize, int nMainDataBegin );
  void ResetDataBitCount();
  int GetDataBitCount() const;
  BOOL SeekDataCountBitsPos( UINT32 nBits );
};

#endif
