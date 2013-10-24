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
  Uint16 nLongBlock[23];      // [Scalefactor Band#] //Max 4bits
  Uint16 nShortBlock[3][13];  // [Window][ScaleFactor Band#] //Max 3bits
};

struct SScalefactorBandIndex
{
    int nLongBlock[23];   //[Scalefactor Band# (=Critical Band#)]
    int nShortBlock[14];  //[Scalefactor Band# (=Critical Band)#]
};


struct SGranuleInfo
{
  Uint16  nPart23Length;                  // 12 Bits
  Uint16  nBigValues;                     // 9 Bits
  Uint16  nGlobalGain;                    // 8 Bits
  Uint16  nScalefacCompress;              // 4 Bits
  Uint16  nWindowSwitchingFlag;           // 1 Bit
  Uint16  nBlockType;                     // 2 Bits 
  Uint16  nMixedBlockFlag;                // 1 Bit
  Uint16  nTableSelect[3];                // 5 Bits * 3
  Uint16  nSubblockGain[3];               // 3 Bits * 3
  Uint16  nRegion0Count;                  // 4 Bits
  Uint16  nRegion1Count;                  // 3 Bits
  Uint16  nPreFlag;                       // 1 Bit 
  Uint16  nScalefacScale;                 // 1 Bit
  Uint16  nCount1TableSelect;             // 1 Bit 0:QUAD_A 1:QUAD_B

  bool IsMixedBlock() const { return nBlockType==2 && nMixedBlockFlag; }
  bool IsLongBlock() const { return nBlockType!=2; }
  bool IsShortBlock() const { return nBlockType==2 && !nMixedBlockFlag; }
};

struct SSampleMap
{
  Uint16 nBlockType;
  Uint16 nSubblock;
  Uint16 nOrderIndex;
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
  Uint16 m_nID;                           //1 bit
  Uint16 m_nLayer;                        //2 bits
  Uint16 m_nProtectionBit;                //1 bit
  Uint16 m_nBitrateIndex;                 //4 bits
  Uint16 m_nSamplingFrequency;            //2 bits
  Uint16 m_nPaddingBit;                   //1 bit
  Uint16 m_nMode;                         //2 bits
  Uint16 m_nModeExtention;                //2 bits
  Uint16 m_nCopyright;                    //1 bit
  Uint16 m_nOriginal;                     //1 bit
  Uint16 m_nEmphasis;                     //2 bitb

  //サイド情報項目
  Uint16  m_nCRCCheck;
  Uint16  m_nMainDataBegin;
  Uint16  m_nScfSi[2][4];                 // [Channel][] 1 Bit
  SGranuleInfo m_GranuleInfo[2][2];       // [Channel][Granule]

  //メイン情報項目
  Uint8   m_pbyteMainData[2048];
  int     m_nMainDataSize;

public:
  CMpegFrame();
  ~CMpegFrame();

  bool IsValidHeader() const;
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
  bool DecodeHuffmanCode( int is[576] );
  void CreateSampleMap(  const SScalefactor* psf, SSampleMap* psm );
  void Dequantize( const SSampleMap* psm, const int is[576], double xr[576] );
  void Reorder( const SSampleMap* psm, double xr[576] );
  void JointStereoDecode( const SScalefactor* psf, double xr[2][576] );
  void Antialias( double lr[576] );
  void ImdctSynthesys( const double lr[576], double pfb[576] );
  void SubbandSynthesys( const double pfb[576], double pfbOut[576] );
  int  CreatePcm( const double pfbOut[2][576], Int16 *pcm );

public:
  CMpegDecoder();
  virtual ~CMpegDecoder();
  void Reset();

  bool DecodeFrame( const CMpegFrame* pMpegFrame, Int16* pcm, int* pnDataNum );
}; 
/*********************end of CMpegDecoder***************/

#endif
