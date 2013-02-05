/****************************************/
// wavファイルを扱うクラスの定義ファイル
/****************************************/

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include "../include/wav.h"


/***************class cWavRead*****************/
cWavRead::cWavRead(const char *wavFileName)
{
  open(wavFileName);
}

cWavRead::~cWavRead()
{
  close();
}

void cWavRead::close()
{
  wavFile.close();
}

void cWavRead::open(const char *wavFileName)
{
  wavFile.open(wavFileName, std::ios::in | std::ios::binary);
  if(wavFile.fail()){
    std::cerr << "Error: Can not open \"" << wavFileName << "\".\n";
    exit(1);
  }

  readHeader();

  finishHeader = true;

  readedSamples = 0;
}

void cWavRead::readFmtChunk()
{
  sWaveFormatPcm tempWaveFormatPcm;
  wavFile.read((char *) &tempWaveFormatPcm.formatTag, 2);
  wavFile.read((char *) &tempWaveFormatPcm.channels, 2);
  wavFile.read((char *) &tempWaveFormatPcm.samplesPerSec, 4);
  wavFile.read((char *) &tempWaveFormatPcm.bytesPerSec, 4);
  wavFile.read((char *) &tempWaveFormatPcm.blockAlign, 2);
  wavFile.read((char *) &tempWaveFormatPcm.bitsPerSample, 2);
  if(wavFile.bad()){
    std::cerr << "Error: Can not read a format chunk.\n";
    exit(1);
  }

  if(tempWaveFormatPcm.channels != 2 && tempWaveFormatPcm.channels != 1){
    std::cerr << "This number of channels is not supported.\n"
	      << "Channels = " << tempWaveFormatPcm.channels << std::endl;
    exit(1);
  }
  if(tempWaveFormatPcm.formatTag != 1){
    std::cerr << "This program supports only linear PCM.\n"
	      << "Format tag = " << tempWaveFormatPcm.formatTag << std::endl;
    exit(1);
  }
  if(tempWaveFormatPcm.bitsPerSample != 8 && tempWaveFormatPcm.bitsPerSample != 16){
    std::cerr << "This program supports only 8/16 bits sampling.\n"
	      << "bits/sec = " << tempWaveFormatPcm.bitsPerSample << std::endl;
    exit(1);
  }

  WrSWaveFileHeader.WaveFormatPcm = tempWaveFormatPcm;
}

// この中身は処理系に依存しないようにバイト単位で指定して読み取っていく
void cWavRead::readHeader()
{
  sWaveFileHeader tempWaveFileHeader;
  wavFile.read((char *) tempWaveFileHeader.hdrRiff, 4);
  wavFile.read((char *) &tempWaveFileHeader.sizeOfFile, 4);
  wavFile.read((char *) tempWaveFileHeader.hdrWave, 4);
  if(wavFile.bad()){
    std::cerr << "Error: can not read the header.\n";
    exit(1);
  }

  if(memcmp( tempWaveFileHeader.hdrRiff, "RIFF", 4) != 0){
    std::cerr << "Error: Not \"RIFF\" format.\n";
    exit(1);
  }

  // WAVE ヘッダ情報
  if (memcmp(tempWaveFileHeader.hdrWave, "WAVE", 4) != 0){
    std::cerr << "Error: Not exist \"WAVE\".\n";
    exit(1);
  }
  
  // 変数メンバに代入
  WrSWaveFileHeader.WaveFileHeader = tempWaveFileHeader;
  
  // 4Byte これ以降のバイト数 = (ファイルサイズ - 8)(Byte)
  int length = tempWaveFileHeader.sizeOfFile;

  // チャンク情報
  
  while(1){
    sChunk tempChunk;
    wavFile.read((char *) tempChunk.hdr, 4);
    wavFile.read((char *) &tempChunk.size,4);

    if(memcmp(tempChunk.hdr, "fmt ", 4) == 0){
      WrSWaveFileHeader.FmtChunk = tempChunk;
      length = tempChunk.size;
      readFmtChunk();
    }
    else if(memcmp(tempChunk.hdr, "data", 4) == 0){
      WrSWaveFileHeader.DataChunk = tempChunk;
      break;
    }
    else{
      std::cerr << "Error: Not supported part exists.\n";
      exit(1);
    }
  }

  numSamplesPerChannel =
    8*WrSWaveFileHeader.DataChunk.size/
    WrSWaveFileHeader.WaveFormatPcm.channels/
    WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
    
}  

// ゲッタ
unsigned cWavRead::channels() const
{
  return WrSWaveFileHeader.WaveFormatPcm.channels;
}

unsigned cWavRead::bitsPerSample() const
{
  return WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
}

unsigned cWavRead::samplingFreq() const
{
  return WrSWaveFileHeader.WaveFormatPcm.samplesPerSec;
}

unsigned cWavRead::fileSize() const
{
  return WrSWaveFileHeader.WaveFileHeader.sizeOfFile;
}

unsigned cWavRead::samplesPerChannel() const
{
  return numSamplesPerChannel;
}

unsigned cWavRead::allNumSamples() const
{
  return numSamplesPerChannel * WrSWaveFileHeader.WaveFormatPcm.channels;
}

void cWavRead::seek(unsigned bytes, std::ios::seekdir point)
{
  wavFile.seekg(bytes, point);
}

// samples分チャネルに関わらず取りだす
int cWavRead::getWave(std::vector<int> &output, const int samples)
{
  if(!finishHeader){
    std::cerr << "Error: Do not access a wave data before read the header info.\n";
    exit(1);
  }

  assert(samples > 0);

  output.clear();
  unsigned bits_per_sample =  WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
  unsigned i;
  if(bits_per_sample == 16){ 
    for(i = 0; i < static_cast<unsigned>(samples); i++){
      if(readedSamples == allNumSamples()){
        break;
      }
      signed short temp;
      wavFile.read((char *) &temp, bits_per_sample/8);
      output.push_back(static_cast<int>(temp));
      readedSamples++;
    }
  }
  else{
    std::cerr << "Error in cWavRead: The number of bits per sample is not available.\n";
    exit(1);
  }

  return i;
}

int cWavRead::getWave(itpp::ivec &output, const int samples)
{
  if(!finishHeader){
    std::cerr << "Error: Do not access a wave data before read the header info.\n";
    exit(1);
  }

  assert(samples > 0);

  output.clear();
  
  unsigned bits_per_sample =  WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
  unsigned i;

  if(bits_per_sample == 16){
    for(i = 0; i < static_cast<unsigned>(samples); i++){
      if(readedSamples == allNumSamples()){
        break;
      }
      signed short temp;
      wavFile.read((char *) &temp, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
      output = itpp::concat(output,static_cast<int>(temp));
      readedSamples++;
    }
  }
  else{
    std::cerr << "Error in cWavRead: The number of bits per sample is not available.\n";
    exit(1);
  }

  return i;
}


int cWavRead::getWave(std::vector<sLRData> &output, const int samples)
{
  if(!finishHeader){
    std::cerr << "Error: Do not access a wave data before read the header info.\n";
    exit(1);
  }

  output.clear();

  unsigned bits_per_sample =  WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
  int ch;
  int i;
  
  if(bits_per_sample == 16){
    if((ch = WrSWaveFileHeader.WaveFormatPcm.channels) == 1){
      for(i = 0; i < samples; i++){
        if(readedSamples == allNumSamples()){
          break;
        }
        sLRData tLRData;
        wavFile.read((char *) &tLRData.Lch, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
        tLRData.Rch = 0;
        output.push_back(tLRData);
        readedSamples++;
      }
    }
    else{
      for(i = 0; i < samples; i++){
        if(readedSamples == numSamplesPerChannel){
          break;
        }
        sLRData tLRData;
        wavFile.read((char *) &tLRData.Lch, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
        wavFile.read((char *) &tLRData.Rch, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
        output.push_back(tLRData);
        readedSamples++;
      }
    }
  }
  else{
    std::cerr << "Error in cWavRead: The number of bits per sample is not available.\n";
    exit(1);
  }

  return i;
}
    
/****************end of cWavRead****************/

/****************class cWavMake*****************/
cWavMake::cWavMake(const char *wavFileName)
{
  open(wavFileName);
}

cWavMake::~cWavMake()
{
  close();
}

void cWavMake::close()
{
  // std::cout << "Wave file was made.\n";
  wavFile.close();
}

void cWavMake::open(const char *wavFileName)
{
  wavFile.open(wavFileName, std::ios::out | std::ios::binary);
  if(wavFile.fail()){
    std::cerr << "Error: Can not open \"" << wavFileName << "\".\n";
    exit(1);
  }
  
}
  
void cWavMake::setupParas(const unsigned samples, // 全体のファイルサイズ - 8
			  const int channels,
			  const int samplingFreq,
			  const int bitsPerSample)
{
  strncpy(WrSWaveFileHeader.WaveFileHeader.hdrRiff, "RIFF", 
	  sizeof(WrSWaveFileHeader.WaveFileHeader.hdrRiff)); // "RIFF"
  
  unsigned dataSize = samples*channels*bitsPerSample/8;

  WrSWaveFileHeader.WaveFileHeader.sizeOfFile = dataSize + 36;
  
  strncpy(WrSWaveFileHeader.WaveFileHeader.hdrWave, "WAVE", 
	  sizeof(WrSWaveFileHeader.WaveFileHeader.hdrWave)); // "WAVE"
  
  strncpy(WrSWaveFileHeader.FmtChunk.hdr, "fmt ",
	  sizeof(WrSWaveFileHeader.FmtChunk.hdr)); // "fmt "

  WrSWaveFileHeader.FmtChunk.size = 16; // fmtチャンクのサイズ

  WrSWaveFileHeader.WaveFormatPcm.formatTag = 1; // リニアPCMなので1

  WrSWaveFileHeader.WaveFormatPcm.channels = channels; // チャンネル数

  WrSWaveFileHeader.WaveFormatPcm.samplesPerSec = samplingFreq; // サンプリング周波数

  int bytes = bitsPerSample/8;
  WrSWaveFileHeader.WaveFormatPcm.bytesPerSec = bytes * samplingFreq * channels; // データ速度

  WrSWaveFileHeader.WaveFormatPcm.blockAlign = bytes * channels; // ブロックサイズ

  WrSWaveFileHeader.WaveFormatPcm.bitsPerSample = bitsPerSample;

  strncpy(WrSWaveFileHeader.DataChunk.hdr, "data", 
	  sizeof(WrSWaveFileHeader.DataChunk.hdr));

  WrSWaveFileHeader.DataChunk.size = dataSize; // 波形データのサイズ
				// 36バイトはヘッダ情報のサイズ

  numSamplesPerChannel = samples;

  //   numSamples =
  //     8*WrSWaveFileHeader.DataChunk.size/
  //     WrSWaveFileHeader.WaveFormatPcm.channels/
  //     WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;

  finishHeader = false;
}

void cWavMake::writeHeader()
{
  // RIFFヘッダ
  wavFile.write((char *) WrSWaveFileHeader.WaveFileHeader.hdrRiff, 4);
  wavFile.write((char *) &WrSWaveFileHeader.WaveFileHeader.sizeOfFile, 4);
  wavFile.write((char *) WrSWaveFileHeader.WaveFileHeader.hdrWave, 4);
    
  // fmtチャンク
  wavFile.write((char *) WrSWaveFileHeader.FmtChunk.hdr, 4);
  wavFile.write((char *) &WrSWaveFileHeader.FmtChunk.size, 4);

  // fmtフォーマット
  wavFile.write((char *) &WrSWaveFileHeader.WaveFormatPcm.formatTag, 2);
  wavFile.write((char *) &WrSWaveFileHeader.WaveFormatPcm.channels, 2);
  wavFile.write((char *) &WrSWaveFileHeader.WaveFormatPcm.samplesPerSec, 4);
  wavFile.write((char *) &WrSWaveFileHeader.WaveFormatPcm.bytesPerSec, 4);
  wavFile.write((char *) &WrSWaveFileHeader.WaveFormatPcm.blockAlign, 2);
  wavFile.write((char *) &WrSWaveFileHeader.WaveFormatPcm.bitsPerSample, 2);
  
  // dataチャンク
  wavFile.write((char *) &WrSWaveFileHeader.DataChunk.hdr, 4);
  wavFile.write((char *) &WrSWaveFileHeader.DataChunk.size, 4);

  if(wavFile.bad()){
    std::cerr << "Error: can not write the header.\n";
    exit(1);
  }

  finishHeader = true;
  currentSamples = 0;
}

void cWavMake::seek(unsigned bytes, std::ios::seekdir point)
{
  wavFile.seekp(bytes, point);
}

bool cWavMake::writeWave(const std::vector<int> &input)
{
  if(!finishHeader){
    std::cerr << "Error: Can not write a wave data before write the header info.\n";
    exit(1);
  }
  
  unsigned bits_per_sample =  WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
  bool eof = false;
  int ch = WrSWaveFileHeader.WaveFormatPcm.channels;
  
  if(bits_per_sample == 16){
    for(unsigned i = 0; i < input.size(); i++){
      signed short temp = static_cast<signed short>(input[i]);
      wavFile.write((char *) &temp, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
      currentSamples++;
    }
  }
  else{
    std::cerr << "Error in cWavMake: not available bits per samples.\n";
    exit(1);
  }
  if(currentSamples == numSamplesPerChannel*ch){
    eof = true;
  }
  

  return eof;
}

bool cWavMake::writeWave(const itpp::ivec &input)
{
  if(!finishHeader){
    std::cerr << "Error: Can not write a wave data before write the header info.\n";
    exit(1);
  }
  
  unsigned bits_per_sample =  WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
  bool eof = false;
  int ch = WrSWaveFileHeader.WaveFormatPcm.channels;
  if(bits_per_sample == 16){
    for(unsigned i = 0; i < static_cast<unsigned>(input.size()); i++){
      signed short temp = static_cast<signed short>(input[i]);
      wavFile.write((char *) &temp, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
      currentSamples++;
    
    }
  }
  else{
    std::cerr << "Error in cWavMake: not available bits per samples.\n";
    exit(1);
  }

  if(currentSamples == numSamplesPerChannel*ch){
    eof = true;
  }
  

  return eof;
}


bool cWavMake::writeWave(const std::vector<sLRData> &input)
{
  if(!finishHeader){
    std::cerr << "Error: Can not write a wave data before write the header info.\n";
    exit(1);
  }
  
  unsigned bits_per_sample =  WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
  bool eof = false;
  int ch;
  
  if(bits_per_sample == 16){
    if((ch = WrSWaveFileHeader.WaveFormatPcm.channels) == 1){
      for(unsigned i = 0; i < input.size(); i++){
        wavFile.write((char *) &input[i].Lch, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
        currentSamples++;
      }
      if(currentSamples == numSamplesPerChannel){
        eof = true;
      }
    }
    else{
      for(unsigned i = 0; i < input.size(); i++){
        wavFile.write((char *) &input[i].Lch, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
        wavFile.write((char *) &input[i].Rch, WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
        currentSamples++;
      }
      if(currentSamples == numSamplesPerChannel){
        eof = true;
      }
    }
  }
  else{
    std::cerr << "Error in cWavMake: not available bits per samples.\n";
    exit(1);
  }

  return eof;
}

bool cWavMake::writeWave(const short input[], const unsigned samples)
{
  if(!finishHeader){
    std::cerr << "Error: Can not write a wave data before write the header info.\n";
    exit(1);
  }
  
  unsigned bits_per_sample =  WrSWaveFileHeader.WaveFormatPcm.bitsPerSample;
  bool eof = false;
  int ch = WrSWaveFileHeader.WaveFormatPcm.channels;
    
  //  if(bMP3){
  //     if((ch = WrSWaveFileHeader.WaveFormatPcm.channels) == 1){
  //       for(unsigned i = 0; i < samples; i++){
  // 	wavFile.write((char *) &input[i], WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
  // 	currentSamples++;
  //       }
  //       if(currentSamples == numSamplesPerChannel){
  // 	eof = true;
  //       }  
  //     }
  //     else{
  //       for(unsigned i = 0; i < samples/ch; i++){
  // 	wavFile.write((char *) &input[i], WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
  // 	wavFile.write((char *) &input[i+samples/2], WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
	
  // 	currentSamples++;
  //       }
  //       if(currentSamples == numSamplesPerChannel){
  // 	eof = true;
  //       }
  //     }
  //   } // if(bMP3)
  //else{

  if(bits_per_sample == 16){
    for(unsigned i = 0; i < samples; i++){
      wavFile.write((char *) &input[i], WrSWaveFileHeader.WaveFormatPcm.bitsPerSample/8);
      currentSamples++;
    }
  }
  else{
    std::cerr << "Error in cWavMake: not available bits per samples.\n";
    exit(1);
  }

  if(currentSamples == numSamplesPerChannel * ch){
    eof = true;
  }  
  //  }
  return eof;
}

  /*
    cWavHandler::cWavHandler(cWavRead wavFile, cWavMake wavMake)
    {
    WavRead = wavFile;
    WavMake = wavMake;
    }

    cWavHandler::cWavHandler(const char *inputFileName, const char *outputFileName)
    {
  
  */
