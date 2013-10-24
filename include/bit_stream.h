#ifndef BIT_STREAM_H
#define BIT_STREAM_H

/***********************************************************/
// MP3のビットストリーム用のクラス
/***********************************************************/



typedef unsigned char Uint8;
typedef unsigned short Uint16;
typedef short Int16;
typedef int Int32;             // at least 32
typedef char Int8;
typedef unsigned int Uint32;

#define FIFO_SIZE 6144

class CBitStream
{
private:
  Uint8   m_pbyte[FIFO_SIZE];
  int     m_nFramePos;
  int     m_nBufBitPos;
  int     m_nBufBytePos;
  int     m_nSize;
  int     m_nDataBitCount;

public:
  CBitStream();
  CBitStream( Uint8* pData, int nSize ) { CBitStream(); SetData( pData, nSize ); }
  virtual ~CBitStream();
  void Reset();

  void SkipBits( int n );
  Uint32 PeekBits( int n );
  Uint32 GetBits( int n );
  Uint16 Get1Bit() { return (Uint16)GetBits(1); }
  void PurgeBits( int nBits ) { GetBits(nBits); }
  bool IsOverrun();
  bool SetData( Uint8* pData, int nSize );

  bool PutFrameMainData( const Uint8* pMainData, int nSize, int nMainDataBegin );
  void ResetDataBitCount();
  int GetDataBitCount() const;
  bool SeekDataCountBitsPos( Uint32 nBits );
};

#endif
