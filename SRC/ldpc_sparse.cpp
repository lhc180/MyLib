// LDPCの本体ファイル

#include <../include/ldpc_sparse.h>
#include <ctime>
#include <cmath>


void ldpc::setEncoder()
{
  // itpp::GF2mat dummyHmat;
  
  makeParityCheckMatrix();

  makeGeneratorMatrix();

  checkParityCondition();
}


// インターリーバを作る関数
static void makeInterleaver(itpp::ivec &interleaver)
{
  itpp::ivec temp(interleaver.size());

  for(int i = 0; i < temp.size(); i++){
    temp[i] = i;
  }
  
  
  for(int i = 0; i < interleaver.size(); i++){
    int num = rand()%temp.size();
    interleaver[i] = temp[num];
    temp.del(num);
  }
}

// パリティ検査行列を作る
void ldpc::makeParityCheckMatrix()
{
  Hmat.set_size(HRowSize, HColSize, HRowWeight); // 検査行列のサイズを設定
  
  // 最初のsubmatrixを生成
  for(int i = 0; i < static_cast<int>(HColSize/HRowWeight); i++){
    for(int j = 0; j < static_cast<int>(HColSize); j++){
      if(j >= i*static_cast<int>(HRowWeight) && j < (i+1)*static_cast<int>(HRowWeight)){
	Hmat.set(i,j,1);
      }
      else{
	Hmat.set(i,j,0);
      }
    }
  }
  
  srand(time(NULL));
  

  // 列重み分のsubmatrixを作る
  for(int submatrixNum = 1; submatrixNum < static_cast<int>(HColWeight); submatrixNum++){
    
    // インターリーバを作る
    itpp::ivec interleaver(HColSize);
    
    makeInterleaver(interleaver);

    for(int i = 0; i < static_cast<int>(HColSize/HRowWeight); i++){
      int rowNum = i + submatrixNum*HColSize/HRowWeight;
      for(int j = 0; j < static_cast<int>(HColSize); j++){
	Hmat.set(rowNum,j,Hmat(i,interleaver(j)));
      }
    }
  } 
  
  
}


// パリティ検査行列を出力する
void ldpc::printParityCheckMatrix()
{
  std::cout << Hmat;
  std::cout << "rank of Hmat is " << Hmat.row_rank() << std::endl;
}


// ガウスの消去法によってPmatを生成
void ldpc::makeGeneratorMatrix()
{
  
  itpp::GF2mat tempHmat = Hmat;	// ガウスの消去法後のHmat
  perm = itpp::zeros_i(Hmat.cols());	// ガウスの消去法の際の列の入れ替え操作を保持しておく

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
  
  Gmat = Imat.concatenate_horizontal(Ptmat);
  // std::cout << "## Matrix G.\n" << Gmat;

  Hmat.permute_cols(perm,0);

  Hmat.compact();
  Gmat.compact();

}

// G*Ht=0かどうかチェックする
void ldpc::checkParityCondition()
{
  itpp::GF2mat MulMat;
  
  MulMat = Gmat * Hmat.transpose();
  
  // std::cout << "Matrix MulMat.\n" << MulMat;
  if(!(MulMat.is_zero())){
    std::cerr << "Error: Generator matrix can not be made.\n";
    exit(1);
  }

}

unsigned ldpc::getCodeLength()
{
  return Hmat.cols();
}

unsigned ldpc::getCheckSymbolLength()
{
  return Hmat.rows();
}

unsigned ldpc::getInputLength()
{
  return Gmat.rows();
}

// 符号化
void ldpc::encoding(const itpp::bvec &input, itpp::bvec &coded)
{
  itpp::GF2mat tempInput(1,Gmat.rows());
  tempInput.set_row(0,input);
  
  itpp::GF2mat tempCoded = tempInput * Gmat;
  coded = tempCoded.bvecify();
  
}

// 符号化率
double ldpc::getRate()
{
  return static_cast<double>(Gmat.rows())/static_cast<double>(Gmat.cols());
}

// 引数の符号を返す
int ldpc::signOfNumber(double input)
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
double ldpc::fFunction(double x)
{

  if(x == 0){
    x = 0.0000000001;
  }
  
  double denom = exp(x)+1;
  double nom = exp(x)-1;

  
  return log(denom/nom);
}


// 行処理
void ldpc::rowsProcessing(itpp::sparse_mat &alpha, itpp::sparse_mat &beta, 
			  const itpp::vec &llrVec,
			  const std::vector<index1> &rowsIndex)
{
  for(int m = 0; m < alpha.rows(); m++){
    //    std::cout << "## m = " << m << "\n";
    for(int n = 0; n < static_cast<int>(rowsIndex[m].indexVec.size()); n++){
      //      std::cout << "  ## n = " << n << "\n";
      int product = 1;
      double sum = 0.0;
      
      for(int i = 0; i < static_cast<int>(rowsIndex[m].indexVec.size()); i++){
	// std::cout << "    ## i = " << i << "\n";
	if(n != i){
	  product *= signOfNumber(llrVec[rowsIndex[m].indexVec[i]] + beta(m,rowsIndex[m].indexVec[i]));
	  sum += fFunction(fabs(llrVec[rowsIndex[m].indexVec[i]] + beta(m,rowsIndex[m].indexVec[i])));
	}
      }
      
      // std::cout << "## sum = " << sum << "\n";
      alpha.set(m,rowsIndex[m].indexVec[n],
		product * fFunction(sum));

    }
  }
}


// 列処理
void ldpc::colsProcessing(itpp::sparse_mat &alpha, itpp::sparse_mat &beta, 
			  const std::vector<index1> &colsIndex)
{
  for(int n = 0; n < beta.cols(); n++){
    for(int m = 0; m < static_cast<int>(colsIndex[n].indexVec.size()); m++){
      double sum = 0.0;
      for(int i = 0; i < static_cast<int>(colsIndex[n].indexVec.size()); i++){
	if(i != m){
	  sum += alpha(colsIndex[n].indexVec[i],n);
	}
      }
      beta.set(colsIndex[n].indexVec[m],n, sum);
      
    }
  }
}

// 一時推定語を求める
void ldpc::estimateCode(const itpp::sprase_mat &alpha, 
			itpp::bvec &decoded,
			const itpp::vec &llrVec,
			const std::vector<index1> &colsIndex)
{
  for(int n = 0; n < decoded.size(); n++){
    double sum = 0.0;
    for(int m = 0; m < static_cast<int>(colsIndex[n].indexVec.size()); m++){
      sum += alpha(colsIndex[n].indexVec[m],n);
    }
    
    if(signOfNumber(llrVec[n] + sum) == 1){
      decoded[n] = 0;
    }
    else{
      decoded[n] = 1;
    }
  }
}

// パリティ検査
bool ldpc::checkParity(itpp::bvec &decoded)
{
  itpp::GF2mat decodedGF2(decoded, false); // GF2matを使って行ベクトルを生成
  itpp::GF2mat mulMat;
  
  mulMat = decodedGF2 * Hmat.transpose();

  bool result; 

  if(mulMat.is_zero()){
    result = true;
  }
  else{
    result = false;
  }

  return result;
}

// 復号
template<typename type>
unsigned ldpc::decoding(const itpp::Modulator< type > &mod,
			const itpp::Vec< type > &symbol,
			itpp::bvec &decodedBits, 
		        const double N0,
		        const unsigned loopMax)
{
  itpp::bvec estimatedCodes = mod.demodulate_bits(symbol);		// sum-product復号法によって得られる推定語

  itpp::sparse_mat beta(Hmat.rows(), Hmat.cols());

  beta.zeros();			// betaを全て0にする

  std::vector<index1> rowsIndex(Hmat.rows()), colsIndex(Hmat.cols()); // 各行、各列の要素が1であるインデックスを
  // 格納する

  //   rowsIndex = new index1[Hmat.rows()];
  //   colsIndex = new index1[Hmat.cols()];

  // 各行、各列の1である要素の列番号を格納していく
  for(int m = 0; m < Hmat.rows(); m++){
    for(int n = 0; n < Hmat.cols(); n++){
      if(Hmat(m,n) == 1){
	rowsIndex[m].indexVec.push_back(n);
	colsIndex[n].indexVec.push_back(m);
      }
    }
  }

  // std::cout << "rowsIndex.indexVec.size = " << rowsIndex[0].indexVec.size() << "\n";
  // std::cout << "colsIndex.indexVec.size = " << colsIndex[0].indexVec.size() << "\n";
  

  itpp::sparse_mat alpha(Hmat.rows(), Hmat.cols());
  itpp::vec llrVec;

  unsigned loop;
  for(loop = 0; loop < loopMax; loop++){
    std::cout << "## loop = " << loop << "\r" << std::flush;

    // alpha.zeros();		// 一応初期化

    llrVec = calcLLR(mod, symbol, estimatedCodes, N0);    

    rowsProcessing(alpha, beta, llrVec, rowsIndex);
    if(loop == 0){
      alpha.compact();
    }
    //    std::cout << "## rows process done.\n";

    colsProcessing(alpha, beta, colsIndex);
    if(loop == 0){
      beta.compact();
    }
//    std::cout << "## cols process done.\n";


    estimateCode(alpha, estimatedCodes, llrVec, colsIndex);

    if(checkParity(estimatedCodes)){
      std::cout << "\n## Success.\n";
      break;
    }
    
    
  }
  std::cout << std::endl;	// ##デバッグ用

  // delete[] rowsIndex;
  //   rowsIndex = 0;
  //   delete[] colsIndex;
  //   colsIndex = 0;

  /*
  itpp::GF2mat decodedMat(1,estimatedCodes.size());
  decodedMat.set_row(0,estimatedCodes);

  //  decodedMat.permute_cols(perm,0);
  itpp::bvec tempDecoded = decodedMat.bvecify();
  */

  decodedBits = estimatedCodes.left(getInputLength()); // 復号語の最初のベクトルが元情報

  return loop;
}

// lambdaを求める
template< typename type >
itpp::vec ldpc::calcLLR(const itpp::Modulator< type > &mod,
			const itpp::Vec< type > &receivedVec,
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
      itpp::Vec< type > tempSymAtk1 = mod.modulate_bits(bitsAtk1);
      itpp::Vec< type > tempSymAtk0 = mod.modulate_bits(bitsAtk0);
      // type型にする
      type symAtk1 = tempSymAtk1[0];
      type symAtk0 = tempSymAtk0[0];
      
      double denom = exp(-pow(abs(receivedVec[i]-symAtk0),2)/N0);
      double nom = exp(-pow(abs(receivedVec[i]-symAtk1),2)/N0);
            
      llrVec[i*bitsPerSymbol + k] = log(denom/nom);
    }
  }

  return llrVec;
}

