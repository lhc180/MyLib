#ifndef MPEG_BINARY_H
#define MPEG_BINARY_H

/*********************************/
// cMpegToBinary
// MP3のオーディオデータ部分(サイド情報と波形データ)
// をバイナリデータにするクラス
// また、バイナリデータからMpegFrameを作る
// CMpegFrameのフレンドクラスである
//
// cBinaryToMpeg
// バイナリデータからCMpegFrameを作る
//
// cMpegToArrangedBinary
// cMpegToBinaryのフレンドクラス
// 複数フレームのバイナリデータをサイドインフォメーションとメインデータをそれぞれまとめる
// それを1ブロックと呼ぶ
//
// cArrangedBinaryToMpeg
// cBinaryToMpegのフレンドクラス
// cMpegToArrangedBinaryの逆
//
// これらはクラスにするほどでも無いけど作ってしまったので
// そのままクラスにしておく
/**********************************/
/********* to do list ************/
// eigenを使ってみる
/******************************/

#include <itpp/itcomm.h>
// #include <vector>
#include "./mpeg.h"

/************** class cMpegToBinay *******************/
class cMpegToBinary{
protected:
  CMpegFrame MpegFrame;		
  itpp::bvec vSideBin;		// サイドインフォのバイナリ
  itpp::bvec vMainBin;		// メインデータのバイナリ

  bool bSet;
  
  void usage();
  
public:
  // コンストラクタ
  cMpegToBinary() : bSet(false) { }  
  cMpegToBinary(const CMpegFrame &tMpegFrame)
  {
    set(tMpegFrame);
  }

  // virtual ~cMpegToBinary() { }; // --- デストラクタ
  void set(const CMpegFrame &tMpegFrame);


  itpp::bvec getSideInfoBinary();
  itpp::bvec getMainDataBinary();
  itpp::bvec getFrameBinary();

};
/****************** end of cMpegToBinary ****************/

/************** class cBinaryToMpeg *****************/
class cBinaryToMpeg{
protected:
  CMpegFrame MpegFrame;		// binaryデータからMpegFrameを返すときに
				// このデータを元にフレームを構成する
  itpp::bvec vSideBin;
  itpp::bvec vMainBin;

  bool bSet;

  void usage();

public:
  cBinaryToMpeg(): bSet(false) { } // --- デフォルトコンストラクタ
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
  // virtual ~cBinaryToMpeg() { } // --- デストラクタ

  // フレーム全体のバイナリデータが入力のとき
  void set(const CMpegFrame &tMpegFrame, const itpp::bvec &vInBin);
  // サイドインフォとメインデータのバイナリデータを別々に入力するとき
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
  // コンストラクタ
  cMpegToArrangedBinary():bSet(false) { }
  cMpegToArrangedBinary(const std::vector<CMpegFrame> &tvecMpegFrame)
  {
    set(tvecMpegFrame);
  }
 
  void set(const std::vector<CMpegFrame> &tvecMpegFrame);
  
  // デフォルトデストラクタ

  itpp::bvec getSideInfoBinary();
  itpp::bvec getMainDataBinary();
  itpp::bvec getBlockBinary();
};

/*********** end of class cMpegToArrangedBinary **************/

/********** class cArrangedBinaryToMpeg ******************/
class cArrangedBinaryToMpeg{
protected:
  std::vector<CMpegFrame> vecMpegFrame;		// binaryデータからMpegFrameを返すときに
  // このデータを元にフレームを構成する
  
  bool bSet;

  void usage();

public:
  // コンストラクタ
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
			
  // virtual ~cArrangedBinaryToMpeg() { } // --- デストラクタ

  // tvecMpegFrameはヘッダなどを用いるだけ
  std::vector<CMpegFrame> getVecMpegFrame();

  
};

#endif
