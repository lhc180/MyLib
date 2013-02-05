/*****************************************/
// オーディオ関連の関数
/****************************************/

#ifndef MYAUDIO_H
#define MYAUDIO_H

#include "bit_stream.h"
#include "mpeg.h"
#include "mpeg_binary.h"
#include "wav.h"
#include "snr.h"
#include "pqmf.h"
#include "quantizer.h"

// system()関数からPQevalAudioを呼び出して、ODGを返す
double calcPEAQ(const char* referenceName, const char *objectName);

// wavファイルからsnrを求める
// 曲長は短い方に合わせる
double calcSNR(const char *xFileName, const char *yFileName, const int frameLength);
// The other name
double calcSDR(const char *xFileName, const char *yFileName, const int frameLength);

// ++++++++++++++++++ 以下はwavを使わない関数++++++++++++++++++++++
// mean square error
double calcMSE(const std::vector<double> &original, const std::vector<double> &object);


#endif
