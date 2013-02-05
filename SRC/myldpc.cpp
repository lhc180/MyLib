// LDPCの本体ファイル

#include "../include/mycomm.h"
#include <mymatrix_utl.h>
#include <ctime>
#include <cmath>
#include <cassert>


void cLDPC::setEncoder(unsigned n, unsigned rowWeight, unsigned colWeight)
{
  HRowSize = n/rowWeight*colWeight;
  HColSize = n;
  HRowWeight= rowWeight;
  HColWeight = colWeight;

  // 列数が行重みの倍数でなかったら失敗
  if(HColSize%HRowWeight != 0 || HColSize < HRowSize){
    std::cout << "You must choose correct code length.\n";
    exit(1);
  }

  execute_setEncoder();
}
  

void cLDPC::execute_setEncoder()
{
  itpp::GF2mat tHmat,tGmat;

  makeParityCheckMatrix(tHmat);

  makeGeneratorMatrix(tHmat,tGmat);

  checkParityCondition(tHmat,tGmat);

  infoLength = tGmat.rows();

  Hmat = mylib::compressGF2matToVector_2D(tHmat);
  Hmat_trans = mylib::compressGF2matToVector_2D(tHmat.transpose());
  
  Gmat = mylib::compressGF2matToVector_2D(tGmat);
  Gmat_trans = mylib::compressGF2matToVector_2D(tGmat.transpose());
  

  bSetDone = true;
}



// パリティ検査行列を作る
void cLDPC::makeParityCheckMatrix(itpp::GF2mat &tHmat)
{
  tHmat.set_size(HRowSize, HColSize); // 検査行列のサイズを設定
  
  // 最初のsubmatrixを生成
  for(int i = 0; i < static_cast<int>(HColSize/HRowWeight); i++){
    for(int j = 0; j < static_cast<int>(HColSize); j++){
      if(j >= i*static_cast<int>(HRowWeight) && j < (i+1)*static_cast<int>(HRowWeight)){
	tHmat.set(i,j,1);
      }
      else{
	tHmat.set(i,j,0);
      }
    }
  }
  
  srand(time(NULL));
  

  // 列重み分のsubmatrixを作る
  for(int submatrixNum = 1; submatrixNum < static_cast<int>(HColWeight); submatrixNum++){
    
    // インターリーバを作る
    itpp::ivec interleaver(HColSize);
    
    mylib::makeInterleaver(interleaver);

    for(int i = 0; i < static_cast<int>(HColSize/HRowWeight); i++){
      int rowNum = i + submatrixNum*HColSize/HRowWeight;
      for(int j = 0; j < static_cast<int>(HColSize); j++){
	tHmat.set(rowNum,j,tHmat(i,interleaver(j)));
      }
    }
  } 
  
  
}


// パリティ検査行列を出力する
// void cLDPC::printParityCheckMatrix()
// {
//   std::cout << Hmat;
//   std::cout << "rank of Hmat is " << Hmat.row_rank() << std::endl;
// }


// ガウスの消去法によってPmatを生成
void cLDPC::makeGeneratorMatrix(itpp::GF2mat &tHmat, itpp::GF2mat &tGmat)
{
  
  itpp::GF2mat tempHmat = tHmat;	// ガウスの消去法後のHmat
  perm = itpp::zeros_i(tHmat.cols());	// ガウスの消去法の際の列の入れ替え操作を保持しておく

  for(int i = 0; i < perm.size(); i++){
    perm(i) = i;
  }

  // /*********ここからは自作のガウスの消去法************/
    
  // rowは行番号でもあり右側行列の対角番号でもある
  // 最右下要素から始める
  for(int index = 0; index < tempHmat.rows(); index++){
    int row = tempHmat.rows() - 1 - index;
    int col = tempHmat.cols() - 1 - index;
    // 単位行列にしたい右側行列の対角要素(row,col)が1になるように
    // 行と列の入れ替えを行なっていく
  
    int i,j;
    for(j = col; j >= 0; j--){
      for(i = row; i >= 0; i--){
	if(tempHmat(i,j) == 1){
	  goto found;
	}
      }
    }
    i = row;
    j = col;
    
  found:
    tempHmat.swap_rows(row,i);
    tempHmat.swap_cols(col,j);

    int temp = perm(j);
    perm(j) = perm(col);
    perm(col) = temp;

    // ここで(row,col)が1になっている
    // (row-1,col)から(0,col)までを0にする

    for(int ref = row-1, m = 1; ref >= 0; ref--, m++){
      // int refCol = Hmat.cols() - col + ref;
      if(tempHmat(ref,col) == 1){
	tempHmat.add_rows(ref,row);
      }
    }
    // (row+1,col)から(Hmat.rows()-1,col)を0にする      
    for(int ref = row + 1; ref < tempHmat.rows(); ref++){
      if(tempHmat(ref,col) == 1){
	tempHmat.add_rows(ref,row);
      }
    }
  }

//   /**********ここまでがガウスの消去法*****************/

  /************ここから余分な行を省く*******************/
  int startRow;
  for(startRow = 0; startRow < tempHmat.rows() - 1; startRow++){
    for(int i = 0; i < tempHmat.cols(); i++){
      if(tempHmat(startRow, i) == 1){
	goto cutting;
      }
    }
  }
  
 cutting:
  tempHmat = tempHmat.get_submatrix(startRow,0,tempHmat.rows()-1,tempHmat.cols()-1);

  // std::cout << "## Matrix H.\n" << Hmat;
  //   std::cout << "## Matrix P.\n" << Pmat;
 
  // 生成行列を作る
  itpp::GF2mat Pmat;

  Pmat = tempHmat.get_submatrix(0,0,tempHmat.rows()-1,tempHmat.cols()-tempHmat.rows() - 1);

  itpp::GF2mat Ptmat = Pmat.transpose();
  
  itpp::GF2mat Imat = itpp::gf2dense_eye(Ptmat.rows());
  
  tGmat = Imat.concatenate_horizontal(Ptmat);
  // std::cout << "## Matrix G.\n" << Gmat;

  tHmat.permute_cols(perm,0);

}

// G*Ht=0かどうかチェックする
void cLDPC::checkParityCondition(itpp::GF2mat &tHmat, itpp::GF2mat &tGmat)
{
  itpp::GF2mat MulMat;
  
  MulMat = tGmat * tHmat.transpose();
  
  // std::cout << "Matrix MulMat.\n" << MulMat;
  if(!(MulMat.is_zero())){
    std::cerr << "Error: Generator matrix can not be made.\n";
    exit(1);
  }

}

unsigned cLDPC::getCodeLength()
{
  return HColSize;
}

unsigned cLDPC::getCheckSymbolLength()
{
  return HRowSize;
}

unsigned cLDPC::getInputLength()
{
  return infoLength;
}

// 符号化
void cLDPC::encode(const itpp::bvec &input, itpp::bvec &coded)
{
  assert(bSetDone);

  coded = mylib::times(input, Gmat_trans); // input * Gmat
  
}

// 符号化率
double cLDPC::getRate()
{
  return static_cast<double>(infoLength)/static_cast<double>(HColSize);
}

// 引数の符号を返す
int cLDPC::signOfNumber(double input)
{
  int sign;
  if(input >= 0){
    sign = 1;
  }
  else{
    sign = -1;
  }

  return sign;
}

// Gallagerのf関数
double cLDPC::fFunction(double x)
{

  if(x == 0){
    x = 0.0000000001;
  }
  
  double denom = exp(x)+1;
  double nom = exp(x)-1;

  
  return log(denom/nom);
}

// 要素番号だけ格納されたelemMatに対応するように値が格納されているmatを転置する
inline mylib::vec_2D transposeCompressedVector_2D(const mylib::vec_2D &mat,
							     const mylib::ivec_2D &elemMat,
							     const mylib::ivec_2D &elemMat_trans)
{
  mylib::vec_2D mat_trans(elemMat_trans.size_rows());
  
  for(int row = 0; row < elemMat.size_rows(); row++){
    for(int col = 0; col < elemMat.size_cols(row); col++){
      int nElem = elemMat(row,col);
      double value = mat(row,col);
      mat_trans.add_cols(nElem, value);
    }
  }

  return mat_trans;
}


// 行処理
void cLDPC::rowsProcessing(mylib::vec_2D &alpha, mylib::vec_2D &beta, 
			   const itpp::vec &llrVec)
{
  for(int m = 0; m < Hmat.size_rows(); m++){
    //    std::cout << "## m = " << m << "\n";
    for(int n = 0; n < Hmat.size_cols(m); n++){
      //      std::cout << "  ## n = " << n << "\n";
      int product = 1;
      double sum = 0.0;
      
      for(int i = 0; i < Hmat.size_cols(m); i++){
	// std::cout << "    ## i = " << i << "\n";
	if(n != i){
	  product *= signOfNumber(llrVec[Hmat(m,i)] + beta(m,i));
	  sum += fFunction(fabs(llrVec[Hmat(m,i)] + beta(m,i)));
	}
      }
      
      // std::cout << "## sum = " << sum << "\n";
      alpha(m,n) = product * fFunction(sum);
      
    }
  }
}


// 列処理
// 一旦
void cLDPC::colsProcessing(mylib::vec_2D &alpha_trans,
			   mylib::vec_2D &beta)
{
  mylib::vec_2D beta_trans(Hmat_trans.size_rows(), Hmat_trans.size_cols());
  
  for(int n = 0; n < Hmat_trans.size_rows(); n++){
    for(int m = 0; m < Hmat_trans.size_cols(n); m++){
      double sum = 0.0;
      for(int i = 0; i < Hmat_trans.size_cols(n); i++){
	if(i != m){
	  sum += alpha_trans(n,i);
	}
      }
      beta_trans(n,m) = sum;
      
    }
  }

  beta = transposeCompressedVector_2D(beta_trans, Hmat_trans, Hmat);
}

// 一時推定語を求める
void cLDPC::estimateCode(const mylib::vec_2D &alpha_trans, 
			 itpp::bvec &decoded,
			 const itpp::vec &llrVec)
{
  std::vector<double> sum(decoded.size());

  for(int n = 0; n < decoded.size(); n++){
    double sum = 0.0;
    for(int m = 0; m < Hmat_trans.size_cols(n); m++){
      sum += alpha_trans(n,m);
    } // for m
    
    if(signOfNumber(llrVec[n] + sum) == 1){
      decoded[n] = 0;
    }
    else{
      decoded[n] = 1;
    }
  } // for n
}

  // パリティ検査
bool cLDPC::checkParity(const itpp::bvec &decoded)
{
  /**********************************************************/
  // ここでやることは
  // mylib::times(decoded, Hmat) --- decoded * Hmat.transpose()
  // の結果がすべて0だったらtrue,そうでなければfalseを返す。
  // これの高速化版
  /**********************************************************/
  bool result = true;
  for(int row = 0; row < Hmat.size_rows(); row++){
    itpp::bin sum = 0;
    for(int col = 0; col < Hmat.size_cols(row); col++){
      int nElem = Hmat(row,col); // 1である要素番号
      assert(nElem < decoded.size()); // ##
      sum += decoded[nElem] * itpp::bin(1);
    } // for col
    
      // sumが1となった時点で終了
    if(sum != 0){
      result = false;
      break;
    }
    
  } // for row

  return result;
}

  // 復号
unsigned cLDPC::decode(const itpp::Modulator_2D &mod,
		       const itpp::cvec &symbol,
		       itpp::bvec &decodedBits, 
		       const double N0,
		       const unsigned loopMax)
{
  assert(bSetDone);

  itpp::bvec estimatedCodes = mod.demodulate_bits(symbol);		// sum-product復号法によって得られる推定語

  assert(Hmat.is_rectangular() && Hmat_trans.is_rectangular());
  
  mylib::vec_2D beta(Hmat.size_rows(), Hmat.size_cols());

  beta.zeros();			// betaを全て0にする

  // std::cout << "rowsIndex.indexVec.size = " << rowsIndex[0].indexVec.size() << "\n";
  // std::cout << "colsIndex.indexVec.size = " << colsIndex[0].indexVec.size() << "\n";
  

  mylib::vec_2D alpha(Hmat.size_rows(), Hmat.size_cols());
  mylib::vec_2D alpha_trans(Hmat_trans.size_rows(), Hmat_trans.size_cols()); // alphaだけ転置行列も用意しとく
  itpp::vec llrVec;

  unsigned loop;
  for(loop = 0; loop < loopMax; loop++){
    // std::cout << "## loop = " << loop << "\r" << std::flush;

    // alpha.zeros();		// 一応初期化

    llrVec = calcLLR(mod, symbol, estimatedCodes, N0);    

    rowsProcessing(alpha, beta, llrVec);
    //    std::cout << "## rows process done.\n";
    alpha_trans = transposeCompressedVector_2D(alpha, Hmat, Hmat_trans);
    
    colsProcessing(alpha_trans, beta);
    //    std::cout << "## cols process done.\n";


    estimateCode(alpha_trans, estimatedCodes, llrVec);

    if(checkParity(estimatedCodes)){
      // std::cout << "\n## Success.\n";
      break;
    }

  }
  
  decodedBits = estimatedCodes.left(getInputLength()); // 復号語の最初のベクトルが元情報

  return loop;
}

  // lambdaを求める
itpp::vec cLDPC::calcLLR(const itpp::Modulator_2D &mod,
			 const itpp::cvec &receivedVec,
			 const itpp::bvec &estimatedVec,
			 const double N0)
{
  int bitsPerSymbol = mod.bits_per_symbol();
  // itpp::Vec< type > symbols = mod.get_symbols();
  // itpp::ivec bits2sym = mod.get_bits2symbols();

  //  std::cout << "\n";
  //   std::cout << symbols << "\n";
  //   std::cout << bits2sym << "\n";

  if(receivedVec.size()*bitsPerSymbol != estimatedVec.size()){
    std::cerr << "Error in ldpc::calcLLR.\n"
	      << "estimatedVec.size() is not correct.\n";
    exit(1);
  }

  itpp::vec llrVec(estimatedVec.size());

  for(int i = 0; i < receivedVec.size(); i++){ 
    // i*bitsPerSymbol目からbitsPerSymbol分取り出す
    itpp::bvec tempBits = estimatedVec.get(i*bitsPerSymbol, i*bitsPerSymbol + bitsPerSymbol -1);
        
    for(int k = 0; k < bitsPerSymbol; k++){
      // tempBitsのkビット目が0と1のものを作る
      itpp::bvec bitsAtk1 = tempBits, bitsAtk0 = tempBits;
      bitsAtk0[k] = 0;
      bitsAtk1[k] = 1;
      // これらをシンボルにする
      // kビット目が0のときと1のときのシンボル(Vecだけど要素数1)
      itpp::cvec tempSymAtk1 = mod.modulate_bits(bitsAtk1);
      itpp::cvec tempSymAtk0 = mod.modulate_bits(bitsAtk0);
      // complex型にする
      std::complex<double> symAtk1 = tempSymAtk1[0];
      std::complex<double> symAtk0 = tempSymAtk0[0];
      
      double nom = exp(-pow(abs(receivedVec[i]-symAtk0),2)/N0);
      double denom = exp(-pow(abs(receivedVec[i]-symAtk1),2)/N0);
            
      llrVec[i*bitsPerSymbol + k] = log(nom/denom);
    }
  }

  return llrVec;
}

  // this is for MLC and MSD.
itpp::vec cLDPC::calcLLR(const std::vector< itpp::Modulator_2D> &vecMod,
			 const itpp::cvec &receivedVec,
			 const double N0)
{
  
  assert(receivedVec.size() == static_cast<int>(vecMod.size()));

  itpp::vec llrVec(receivedVec.size());

  for(int i = 0; i < static_cast<int>(receivedVec.size()); i++){

    itpp::cvec symbols = vecMod[i].get_symbols();
    itpp::ivec bitmap = vecMod[i].get_bits2symbols();
    int nSymbols = symbols.size();

    double nom = 0.0;		// 分子
    double denom = 0.0;		// 分母
    // vecMod[i]の各symbolsの最下位ビットが0か1かで分母と分子の計算を分ける
    for(int s = 0; s < nSymbols; s++){
      if( (bitmap[s] & 1) == 1){
	denom += exp(-pow(abs(symbols[s] - receivedVec[i]),2)/N0);
      }
      else{
	nom += exp(-pow(abs(symbols[s] - receivedVec[i]),2)/N0);
      }
    } // for s
    
    llrVec[i] = log(nom/denom);

  } // for i
    
  return llrVec;
}

  // this is for MLC and MSD.
unsigned cLDPC::decode(const std::vector< itpp::Modulator_2D > &vecMod,
		       const itpp::cvec &symbol,
		       itpp::bvec &decodedCodes, 
		       const double N0,
		       const unsigned loopMax)
{
  assert(bSetDone);
  
  
  itpp::bvec estimatedCodes(symbol.size()); // sum-product復号法によって得られる推定語
  
  assert(Hmat.is_rectangular() && Hmat_trans.is_rectangular());

  mylib::vec_2D beta(Hmat.size_rows(), Hmat.size_cols());

  beta.zeros();			// betaを全て0にする  

  mylib::vec_2D alpha(Hmat.size_rows(), Hmat.size_cols());
  mylib::vec_2D alpha_trans(Hmat_trans.size_rows(), Hmat_trans.size_cols()); // alphaだけ転置行列も用意しとく
  itpp::vec llrVec = calcLLR(vecMod, symbol, N0);    
  
  unsigned loop;
  for(loop = 0; loop < loopMax; loop++){
    // std::cout << "## loop = " << loop << "\r" << std::flush;

    // alpha.zeros();		// 一応初期化

    

    rowsProcessing(alpha, beta, llrVec);
    //    std::cout << "## rows process done.\n";
    alpha_trans = transposeCompressedVector_2D(alpha, Hmat, Hmat_trans);

    
    colsProcessing(alpha_trans, beta);
    //    std::cout << "## cols process done.\n";


    estimateCode(alpha_trans, estimatedCodes, llrVec);

    if(checkParity(estimatedCodes)){
      // std::cout << "\n## Success.\n";
      break;
    }

  }
  // std::cout << std::endl;	// ##デバッグ用

  // delete[] rowsIndex;
  //   rowsIndex = 0;
  //   delete[] colsIndex;
  //   colsIndex = 0;


  //  decodedMat.permute_cols(perm,0);
  decodedCodes = estimatedCodes;

  return loop;
}  
