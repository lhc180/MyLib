#include <string>
#include <fstream>
#include <cassert>
#include <cstdlib>
// #include "../include/wav.h"
#include "../include/myaudio.h"

const std::string PQ_COMMAND("PQevalAudio");
const std::string TEMP_FILENAME("_temp.dat");


// PQevalAudio���C���X�g�[�����Ă��邱�Ƃ��}�X�g�I
double calcPEAQ(const char *reference, const char *object)
{
  std::string referenceName(reference);
  std::string objectName(object);

  std::string pqCommand(PQ_COMMAND.c_str());

  pqCommand += ' ';
  pqCommand += referenceName;
  pqCommand += ' ';
  pqCommand += objectName;
  pqCommand += " > ";
  pqCommand += TEMP_FILENAME;
  
  // �R�}���h���s
  if(system(pqCommand.c_str()) != 0){
    std::cerr << "Error at calcPEAQ.\n"
	      << "Can not execute \"" << pqCommand << "\"\n";
    exit(1);
  }

  // �����temp.dat���o���āA����PEAQ ODG�̒l������
  std::ifstream tempFile(TEMP_FILENAME.c_str());
  assert(!tempFile.fail());

  tempFile.seekg(-7,std::ios::end);
  
  double odg;
  tempFile >> odg;

  tempFile.close();

  // �쐬����temp.dat���폜
  std::string rmCommand("rm ");
  
  rmCommand += TEMP_FILENAME;

  if(system(rmCommand.c_str()) != 0){
    std::cerr << "Error at calcPEAQ.\n"
	      << "Can not execute \"" << rmCommand << "\"\n";
    exit(1);
  }

  return odg;
}
  
// WAV�t�@�C�����m��SNR�����߂�
double calcSNR(const char *xFileName, const char *yFileName, const int frameLength)
{
  cWavRead xFile(xFileName);
  cWavRead yFile(yFileName);

  assert(xFile.channels() == yFile.channels());
  assert(xFile.bitsPerSample() == yFile.bitsPerSample());
  assert(xFile.samplingFreq() == yFile.samplingFreq());

  cSNR SNR;
  
  bool bFinalFrame = false;
  while(1){
    int xFrameLength, yFrameLength;
    itpp::ivec x, y;
    xFrameLength = xFile.getWave(x, frameLength*xFile.channels());
    yFrameLength = yFile.getWave(y, frameLength*yFile.channels());
    // frameLength�����Ȃ�������Ō�̃t���[���Ƃ�������
    if(frameLength*static_cast<int>(xFile.channels()) != xFrameLength){
      bFinalFrame = true;
    }
    if(frameLength*static_cast<int>(yFile.channels()) != yFrameLength){
      bFinalFrame = true;
    }
    // �������Z�����ɍ��킹��
    if(xFrameLength > yFrameLength){
      x = x.left(yFrameLength);
    }
    else if(xFrameLength < yFrameLength){
      y = y.left(xFrameLength);
    }

    SNR.calc(x,y);

    if(bFinalFrame){
      break;
    }
  } // while 1

  return SNR.getSegmentalSNR();
}

// Just the other name
double calcSDR(const char *xFileName, const char *yFileName, const int frameLength)
{
  return calcSNR(xFileName, yFileName, frameLength);
}

// mean square error
double calcMSE(const std::vector<double> &original, const std::vector<double> &object)
{
  assert(original.size() == object.size());

  double mse = 0.0;
  for(int i = 0; i < static_cast<int>(original.size()); i++){
    mse += pow(object[i] - original[i], 2);
  }

  mse /= original.size();

  return mse;
}

