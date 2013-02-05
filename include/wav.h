#ifndef WAV_H
#define WAV_H

/******************************************/
// class cWavRead
// Wavファイルをより扱いやすくするためのクラス
//
// class cWavMake
// Wavファイルへの出力をより簡単なものにするクラス
//
// class cWavHandler
// 前2つのクラスを所持して、vectorやテキストファイルの入出力を行うクラス
// 
// wavからvector,vectorからwav,あとテキストファイルの入出力にも対応
// これまではデータ数分の配列を用意していたが、順次必要なデータ数を取り出していくことにより
// メモリを気にせずにwavファイルの入出力が行える
//
// リニアPCMにしか対応していない。
/*****************************************/
// --------- update info -----------
// 2012/04/01
// 仕様変更
// 特にcWavMakeのwriteWaveの仕様変更は大きい

#include <fstream>
#include <vector>
#include <itpp/itbase.h>

struct sWaveFileHeader
{
  char            hdrRiff[4];         // 'RIFF'
  unsigned int    sizeOfFile;         // ファイルサイズ - 8
  char            hdrWave[4];         // 'WAVE'
};

struct sChunk
{
  char   hdr[4];      // 'fmt ' or 'data'
  unsigned int    size;      // sizeof(PCMWAVEFORMAT) or Waveデーターサイズ
};

struct sWaveFormatPcm
{
  unsigned short  formatTag;          // WAVE_FORMAT_PCM
  unsigned short  channels;           // number of channels
  unsigned int    samplesPerSec;      // sampling rate
  unsigned int    bytesPerSec;        // samplesPerSec * channels * (bitsPerSample/8)
  unsigned short  blockAlign;         // block align
  unsigned short  bitsPerSample;      // bits per sampling
};

struct sWrSWaveFileHeader
{
  sWaveFileHeader WaveFileHeader;
  sChunk FmtChunk;
  sWaveFormatPcm WaveFormatPcm;
  sChunk DataChunk;
};

struct sLRData
{
  short Lch;
  short Rch;
};


class cWavRead
{
private:
  std::ifstream wavFile;		// wavファイル
  sWrSWaveFileHeader WrSWaveFileHeader; // ヘッダ情報
  unsigned numSamplesPerChannel;	// 波形データの1チャネル当たりの総サンプル数
  unsigned readedSamples;
  bool finishHeader;			// ヘッダの読み込みが終わっているかどうか

  void readHeader();		// ヘッダ情報の読み取り
  void readFmtChunk();		// フォーマットチャンクの読み取り
				// readHeader()から呼び出される
public:
  cWavRead(): finishHeader(false){}
  cWavRead(const char *wavFileName); // コンストラクタ
  virtual ~cWavRead();

  void open(const char *wavFileName);
  // ゲッタ
  unsigned channels() const;
  unsigned bitsPerSample() const;
  unsigned samplingFreq() const;
  unsigned fileSize() const;		// ファイルサイズ - 8
  unsigned samplesPerChannel() const;		// 総サンプル数
					// ただし1チャンネル分の数
  
  unsigned allNumSamples() const;		// 総サンプル数

  // fseekみたいに使う
  void seek(unsigned bytes, std::ios::seekdir point);

  // 波形データをsamples分読み取る
  // 読み取れたサンプル数を返す
  int getWave(std::vector<int> &output, const int samples); // チャネルに関わらずそのまま取りだす
  int getWave(itpp::ivec &output, const int samples);
  int getWave(std::vector<sLRData> &output, const int samples);

  // bool fail(); --- open関数の中に実装してあるので作る必要がない
  void close(); 		

  // friend class cWavHandler;
};

class cWavMake
{
private:
  std::ofstream wavFile;
  sWrSWaveFileHeader WrSWaveFileHeader; // ヘッダ情報
  bool finishHeader;			// ヘッダの書き込みが終了しているかどうか
  unsigned numSamplesPerChannel;
  unsigned currentSamples;	// 波形書き込み時の現在のサンプル番号を格納
  

public:
  // コンストラクタ
  cWavMake():finishHeader(false) {}
  cWavMake(const char *wavFileName); 
  virtual ~cWavMake();
  
  void open(const char *wavFileName);

  // ヘッダ情報を設定
  // samplesには1チャンネル当たりのサンプル数を入れる
  void setupParas(const unsigned samplesPerChannel = 0,
		  const int channels = 1,
		  const int samplingFreq = 44100,
		  const int bitsPerSample = 16);

  // ヘッダを書きこむ
  void writeHeader();

  // ファイルのポインタをファイルの始めに戻す
  // 最後にヘッダを書きかえるのに使える
  void seek(unsigned bytes, std::ios::seekdir point);
  // numSamplesまで達したらtrueを返す
  bool writeWave(const std::vector<int> &input);
  bool writeWave(const itpp::ivec &input);
  bool writeWave(const std::vector<sLRData> &input);
  // MP3のデコーダーから出力されたpcmを入れる場合はtrue
  // WAVの形式だとfalseにする
  bool writeWave(const short input[], const unsigned samples);
  
  void close();

  //  friend class cWavHandler;
};

/* ここからはとりあえず作る必要無いので放置
// 入力はwav,txt,vector
// 出力も同様
class cWavHandler
{
private:
  cWavRead WavRead;
  cWavMake WavMake;
  std::ifstream inputFile;
  std::ofstream outputFile;

public:
  // コンストラクタ
  cWavHandler() {};
  cWavHandler(const cWavRead wavRead, const cWavMake wavMake);		
  cWavHandler(const char *inputFileName, const char *outputFileName);

  setInput(cWavRead wavRead);
  setInput(const char *inputFileName);

  // デストラクタ
  ~cWavHandler();

}
*/

#endif
