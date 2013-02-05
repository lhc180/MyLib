#ifndef MPEG_H
#define MPEG_H

/*************************************************/
// MP3のデコード用のクラス一式
// 以下のホームページのソースコードをほぼそのまま転用
// http://web.archive.org/web/20050315005127/http://www.geocities.jp/bywnn498/mp3/
//
// C++とかいいながら、クラス以外ほとんどC言語で書かれている 
/**************************************************/
#include <cstdio>
#include "./bit_stream.h"

struct SScalefactor
{
  UINT16 nLongBlock[23];      // [Scalefactor Band#] //Max 4bits
  UINT16 nShortBlock[3][13];  // [Window][ScaleFactor Band#] //Max 3bits
};

struct SScalefactorBandIndex
{
    int nLongBlock[23];   //[Scalefactor Band# (=Critical Band#)]
    int nShortBlock[14];  //[Scalefactor Band# (=Critical Band)#]
};


struct SGranuleInfo
{
  UINT16  nPart23Length;                  // 12 Bits
  UINT16  nBigValues;                     // 9 Bits
  UINT16  nGlobalGain;                    // 8 Bits
  UINT16  nScalefacCompress;              // 4 Bits
  UINT16  nWindowSwitchingFlag;           // 1 Bit
  UINT16  nBlockType;                     // 2 Bits 
  UINT16  nMixedBlockFlag;                // 1 Bit
  UINT16  nTableSelect[3];                // 5 Bits * 3
  UINT16  nSubblockGain[3];               // 3 Bits * 3
  UINT16  nRegion0Count;                  // 4 Bits
  UINT16  nRegion1Count;                  // 3 Bits
  UINT16  nPreFlag;                       // 1 Bit 
  UINT16  nScalefacScale;                 // 1 Bit
  UINT16  nCount1TableSelect;             // 1 Bit 0:QUAD_A 1:QUAD_B

  BOOL IsMixedBlock() const { return nBlockType==2 && nMixedBlockFlag; }
  BOOL IsLongBlock() const { return nBlockType!=2; }
  BOOL IsShortBlock() const { return nBlockType==2 && !nMixedBlockFlag; }
};

struct SSampleMap
{
  UINT16 nBlockType;
  UINT16 nSubblock;
  UINT16 nOrderIndex;
  double fScale;
};

/*************************************************************/
// class CMpegFrame
// MP3フレームに相当するクラス
// ヘッダー、サイド情報のデコードを行い、取得したデータを保持する
//
// ライフタイムは1フレーム処理の間
/**************************************************************/
class CMpegFrame
{
private:
  //フレームヘッダー項目
  UINT16 m_nID;                           //1 bit
  UINT16 m_nLayer;                        //2 bits
  UINT16 m_nProtectionBit;                //1 bit
  UINT16 m_nBitrateIndex;                 //4 bits
  UINT16 m_nSamplingFrequency;            //2 bits
  UINT16 m_nPaddingBit;                   //1 bit
  UINT16 m_nMode;                         //2 bits
  UINT16 m_nModeExtention;                //2 bits
  UINT16 m_nCopyright;                    //1 bit
  UINT16 m_nOriginal;                     //1 bit
  UINT16 m_nEmphasis;                     //2 bitb

  //サイド情報項目
  UINT16  m_nCRCCheck;
  UINT16  m_nMainDataBegin;
  UINT16  m_nScfSi[2][4];                 // [Channel][] 1 Bit
  SGranuleInfo m_GranuleInfo[2][2];       // [Channel][Granule]

  //メイン情報項目
  UINT8   m_pbyteMainData[2048];
  int     m_nMainDataSize;

public:
  CMpegFrame();
  ~CMpegFrame();

  BOOL IsValidHeader() const;
  int GetFrameSize() const;
  int GetMainDataSize() const;
  int GetSampleNum() const;
  int GetFrequency() const; // Hz
  int GetBitrate() const; // Bits/Sec


  int GetID() const { return m_nID; }
  int GetLayer() const { return m_nLayer; }
  int GetProtectionBit() const { return m_nProtectionBit; }
  int GetBitrateIndex() const { return m_nBitrateIndex; }
  int GetSamplingFrequency() const { return m_nSamplingFrequency; }
  int GetPaddingBit() const { return m_nPaddingBit; }
  int GetMode() const { return m_nMode; }
  int GetModeExtention() const { return m_nModeExtention; }
  int GetCopyright() const { return m_nCopyright; }
  int GetOriginal() const { return m_nOriginal; }
  int GetEmphasis() const { return m_nEmphasis; }

  int GetVersion() const { return m_nID==0 ? 2 : 1; }
  int GetChannels() const { return m_nMode==3 ? 1 : 2; }
  int GetCRCCheck() const { return m_nCRCCheck; }

  friend class CMpegFile;
  friend class CMpegDecoder;

  // オリジナルクラス
  friend class cMpegToBinary;
  friend class cBinaryToMpeg;
};
/******************end of CMpegFrame************************/


/****************************************************************/
// class CMpegFile
// MP3ファイルに相当するクラス
// MP3ファイルにアクセスし、フレームを順次に読み込む
//
// ライフタイムはMP3に対するアクセスの間
/****************************************************************/
class CMpegFile
{
private:
  FILE* m_pFile;

  bool DecodeSideInfo( CBitStream* pbs, CMpegFrame* pmf );

public:
  CMpegFile();
  ~CMpegFile();

  bool Open( const char* pFilename );
  void Close();
  bool IsOpen() const { return m_pFile!=NULL; }
  bool SeekMpegFrame( CMpegFrame* pMpegFrame );
};
/***************end of CMpegFile*******************************/

/***********************************************************/
// class CMpegDecoder
// MP3デコーダーエンジンに相当するクラス
// メインデータのデコードを行い、PCMサンプルを生成する
//
// ライフタイムはデコードの開始から終了まで
/**********************************************************/
class CMpegDecoder
{
private:
  CBitStream m_bs;              //ビットストリームオブジェクト

  const CMpegFrame* m_pMpegFrame; //処理中のMPEGフレーム
  int               m_nGranule; //処理中のグラニュール(0 or 1)
  int               m_nChannel; //処理中のチャネル(0 or 1)

  //ImdctSynthesys用バッファ
  double m_fImdctPrevRawout[2][576]; //[Channel][Sample]

  //SubbandSynthesys用バッファ
  double m_fSubbandBuf[2][16][64]; //[ch][buf#][64samples]
  int    m_nSubbandBufIndex[2]; //[ch]

  int m_nLastError;

private:
  int  ProcessIStereo( const SScalefactor* psf, double lr[2][576] );
  void ProcessMSStereo( double x[2][576], int nISIndex );

  void DecodeScalefactors( SScalefactor sf[2][2] );
  BOOL DecodeHuffmanCode( int is[576] );
  void CreateSampleMap(  const SScalefactor* psf, SSampleMap* psm );
  void Dequantize( const SSampleMap* psm, const int is[576], double xr[576] );
  void Reorder( const SSampleMap* psm, double xr[576] );
  void JointStereoDecode( const SScalefactor* psf, double xr[2][576] );
  void Antialias( double lr[576] );
  void ImdctSynthesys( const double lr[576], double pfb[576] );
  void SubbandSynthesys( const double pfb[576], double pfbOut[576] );
  int  CreatePcm( const double pfbOut[2][576], INT16 *pcm );

public:
  CMpegDecoder();
  virtual ~CMpegDecoder();
  void Reset();

  BOOL DecodeFrame( const CMpegFrame* pMpegFrame, INT16* pcm, int* pnDataNum );
}; 
/*********************end of CMpegDecoder***************/

#endif
