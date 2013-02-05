#ifndef PQMF_H
#define PQMF_H

/**********************************/
// PQMF
/**********************************/
/*********to do list***************/
// ��������̃f�[�^��itpp::mat�ɑΉ�������
// clear�֐�����������
/*********************************/

#include<vector>
#include<assert.h>
#include <itpp/base/vec.h>
#include "./mymatrix.h"

const unsigned FILTER_LENGTH = 512;
const unsigned SUBBAND_NUM = 32;


class cPQMF{
private:
  double x[FILTER_LENGTH],y[1024]; // ���̓t�B���^�o�b�t�@�ƍ����t�B���^�̃o�b�t�@
  double m[SUBBAND_NUM][64],n[SUBBAND_NUM][64],c[FILTER_LENGTH];
  bool start;			// �����t�B���^���X�^�[�g���ǂ���
  
public:
  cPQMF();			// �R���X�g���N�^
  
				// �f�t�H���g�f�X�g���N�^

  // �R�s�[�R���X�g���N�^�̓f�t�H���g

  void analyze(const std::vector<double> &input, std::vector<double> &output);
  void analyze(const itpp::vec &input, itpp::vec &output);
  void analyze(const std::vector<double> &input, mylib::vec_2D &output_2D);


//   void fastAnalysis(const std::vector<double> &input, std::vector<double> &output);

  // 481�T���v���̒x��������̂Œ���
  void synthesize(const std::vector<double> &input, std::vector<double> &output);
  void synthesize(const itpp::vec &input, itpp::vec &output);
  void synthesize(const mylib::vec_2D &input_2D, std::vector<double> &output);
  
//   void fastSynthesis(const std::vector<double> &input, std::vector<double> &output);
  
};



#endif
