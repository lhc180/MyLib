#ifndef MPEG_BINARY_H
#define MPEG_BINARY_H

/*********************************/
// cMpegToBinary
// MP3�̃I�[�f�B�I�f�[�^����(�T�C�h���Ɣg�`�f�[�^)
// ���o�C�i���f�[�^�ɂ���N���X
// �܂��A�o�C�i���f�[�^����MpegFrame�����
// CMpegFrame�̃t�����h�N���X�ł���
//
// cBinaryToMpeg
// �o�C�i���f�[�^����CMpegFrame�����
//
// cMpegToArrangedBinary
// cMpegToBinary�̃t�����h�N���X
// �����t���[���̃o�C�i���f�[�^���T�C�h�C���t�H���[�V�����ƃ��C���f�[�^�����ꂼ��܂Ƃ߂�
// �����1�u���b�N�ƌĂ�
//
// cArrangedBinaryToMpeg
// cBinaryToMpeg�̃t�����h�N���X
// cMpegToArrangedBinary�̋t
//
// �����̓N���X�ɂ���قǂł��������Ǎ���Ă��܂����̂�
// ���̂܂܃N���X�ɂ��Ă���
/**********************************/
/********* to do list ************/
// eigen���g���Ă݂�
/******************************/

#include <itpp/itcomm.h>
// #include <vector>
#include "./mpeg.h"

/************** class cMpegToBinay *******************/
class cMpegToBinary{
protected:
  CMpegFrame MpegFrame;		
  itpp::bvec vSideBin;		// �T�C�h�C���t�H�̃o�C�i��
  itpp::bvec vMainBin;		// ���C���f�[�^�̃o�C�i��

  bool bSet;
  
  void usage();
  
public:
  // �R���X�g���N�^
  cMpegToBinary() : bSet(false) { }  
  cMpegToBinary(const CMpegFrame &tMpegFrame)
  {
    set(tMpegFrame);
  }

  // virtual ~cMpegToBinary() { }; // --- �f�X�g���N�^
  void set(const CMpegFrame &tMpegFrame);


  itpp::bvec getSideInfoBinary();
  itpp::bvec getMainDataBinary();
  itpp::bvec getFrameBinary();

};
/****************** end of cMpegToBinary ****************/

/************** class cBinaryToMpeg *****************/
class cBinaryToMpeg{
protected:
  CMpegFrame MpegFrame;		// binary�f�[�^����MpegFrame��Ԃ��Ƃ���
				// ���̃f�[�^�����Ƀt���[�����\������
  itpp::bvec vSideBin;
  itpp::bvec vMainBin;

  bool bSet;

  void usage();

public:
  cBinaryToMpeg(): bSet(false) { } // --- �f�t�H���g�R���X�g���N�^
  cBinaryToMpeg(const CMpegFrame &tMpegFrame, const itpp::bvec &vInBin)
  {
    set(tMpegFrame, vInBin);
  }

  cBinaryToMpeg(const CMpegFrame &tMpegFrame,
		const itpp::bvec &vInSideBin,
		const itpp::bvec &vInMainBin)
  {
    set(tMpegFrame, vInSideBin, vInMainBin);
  }
  // virtual ~cBinaryToMpeg() { } // --- �f�X�g���N�^

  // �t���[���S�̂̃o�C�i���f�[�^�����͂̂Ƃ�
  void set(const CMpegFrame &tMpegFrame, const itpp::bvec &vInBin);
  // �T�C�h�C���t�H�ƃ��C���f�[�^�̃o�C�i���f�[�^��ʁX�ɓ��͂���Ƃ�
  void set(const CMpegFrame &tMpegFrame,
	   const itpp::bvec &vInSideBin,
	   const itpp::bvec &vInMainBin);

  CMpegFrame getMpegFrame();

};
/************** end of cBinaryToMpeg ****************/

/************* class cMpegToArrangedBinary **********/
class cMpegToArrangedBinary{
protected:
  itpp::bvec vSideBin;
  itpp::bvec vMainBin;

  bool bSet;

  void usage();

public:
  // �R���X�g���N�^
  cMpegToArrangedBinary():bSet(false) { }
  cMpegToArrangedBinary(const std::vector<CMpegFrame> &tvecMpegFrame)
  {
    set(tvecMpegFrame);
  }
 
  void set(const std::vector<CMpegFrame> &tvecMpegFrame);
  
  // �f�t�H���g�f�X�g���N�^

  itpp::bvec getSideInfoBinary();
  itpp::bvec getMainDataBinary();
  itpp::bvec getBlockBinary();
};

/*********** end of class cMpegToArrangedBinary **************/

/********** class cArrangedBinaryToMpeg ******************/
class cArrangedBinaryToMpeg{
protected:
  std::vector<CMpegFrame> vecMpegFrame;		// binary�f�[�^����MpegFrame��Ԃ��Ƃ���
  // ���̃f�[�^�����Ƀt���[�����\������
  
  bool bSet;

  void usage();

public:
  // �R���X�g���N�^
  cArrangedBinaryToMpeg() : bSet(false) { } 
  cArrangedBinaryToMpeg(const std::vector<CMpegFrame> &tvecMpegFrame,
			const itpp::bvec &vInBin)
  {
    set(tvecMpegFrame, vInBin);
  }
  cArrangedBinaryToMpeg(const std::vector<CMpegFrame> &tvecMpegFrame,
			const itpp::bvec &vInSideBin,
			const itpp::bvec &vInMainBin)
  {
    set(tvecMpegFrame, vInSideBin, vInMainBin);
  }

  void set(const std::vector<CMpegFrame> &tvecMpegFrame,
	   const itpp::bvec &vInBin);

  void set(const std::vector<CMpegFrame> &tvecMpegFrame,
	   const itpp::bvec &vInSideBin,
	   const itpp::bvec &vInMainBin);
			
  // virtual ~cArrangedBinaryToMpeg() { } // --- �f�X�g���N�^

  // tvecMpegFrame�̓w�b�_�Ȃǂ�p���邾��
  std::vector<CMpegFrame> getVecMpegFrame();

  
};

#endif
