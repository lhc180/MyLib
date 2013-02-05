/*****************************************/
// �I�[�f�B�I�֘A�̊֐�
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

// system()�֐�����PQevalAudio���Ăяo���āAODG��Ԃ�
double calcPEAQ(const char* referenceName, const char *objectName);

// wav�t�@�C������snr�����߂�
// �Ȓ��͒Z�����ɍ��킹��
double calcSNR(const char *xFileName, const char *yFileName, const int frameLength);
// The other name
double calcSDR(const char *xFileName, const char *yFileName, const int frameLength);

// ++++++++++++++++++ �ȉ���wav���g��Ȃ��֐�++++++++++++++++++++++
// mean square error
double calcMSE(const std::vector<double> &original, const std::vector<double> &object);


#endif
