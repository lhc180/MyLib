// LDPCの本体ファイル

#include "../include/mycomm.h"
// #include <mymatrix_utl.h>
#include <ctime>
#include <cmath>
#include <cassert>

namespace mylib{
  void Ldpc::Set(unsigned n, unsigned rowWeight, unsigned colWeight)
  {
    hRowSize_ = n/rowWeight*colWeight;
    hColSize_ = n;
    hRowWeight_= rowWeight;
    hColWeight_ = colWeight;

    // 列数が行重みの倍数でなかったら失敗
    if(hColSize_%hRowWeight_ != 0 || hColSize_ < hRowSize_){
      std::cerr << "You must choose correct code length.\n";
      std::cerr << "Code length = " << n << std::endl;
      std::cerr << "Row weight = " << rowWeight << std::endl;
      exit(1);
    }

    SetupMatrix();
    setDone_ = true;
  }
  
  inline std::vector<std::vector < int > > CompressGF2matToVectorVector(const itpp::GF2mat& input, int capa)
  {
    std::vector< std::vector < int > > mat(input.rows());
    
    for(int row = 0; row < input.rows(); ++row){
      mat[row].reserve(capa);
      for(int col = 0; col < input.cols(); ++col){
	if(input(row,col) == 1){
	  mat[row].push_back(col);
	}
      }
    }

    return mat;
  }
  
  inline void Ldpc::SetupMatrix()
  {
    itpp::GF2mat t_hMat, t_gMat;

    MakeParityCheckMatrix(t_hMat);

    MakeGeneratorMatrix(t_hMat,t_gMat);

    CheckParityCondition(t_hMat,t_gMat);
  
    infoLength_ = t_gMat.rows();

    hMat_ = CompressGF2matToVectorVector(t_hMat, hRowWeight_);
    hMatTrans_ = CompressGF2matToVectorVector(t_hMat.transpose(), hColWeight_);
  
    gMat_ = CompressGF2matToVectorVector(t_gMat, t_gMat.cols()/2);
    gMatTrans_ = CompressGF2matToVectorVector(t_gMat.transpose(), t_gMat.rows()/2);

    setDone_ = true;
  }

  // パリティ検査行列を作る
  inline void Ldpc::MakeParityCheckMatrix(itpp::GF2mat &t_hMat)
  {
    t_hMat.set_size(hRowSize_, hColSize_); // 検査行列のサイズを設定

    // 最初のsubmatrixを生成
    for(int i = 0; i < static_cast<int>(hColSize_/hRowWeight_); i++){
      for(int j = i*hRowWeight_; j < (i+1)*hRowWeight_; j++){
	t_hMat.set(i,j,1);
      } // for j
    }   // for i

    // 列重み分のsubmatrixを作る
    // インターリーバを作る
    itpp::Sequence_Interleaver<int> sequence_interleaver(hColSize_); // intはダミー
    for(int submatrixNum = 1; submatrixNum < hColWeight_; submatrixNum++){
      sequence_interleaver.randomize_interleaver_sequence();
      itpp::ivec interleaver = sequence_interleaver.get_interleaver_sequence();
            
      for(int i = 0; i < static_cast<int>(hColSize_/hRowWeight_); i++){
        int rowNum = i + submatrixNum*hColSize_/hRowWeight_;
        for(int j = 0; j < static_cast<int>(hColSize_); j++){
          t_hMat.set(rowNum,j,t_hMat(i,interleaver[j]));
        } // for j
      }   // for i
    }     // for submatrixNum
  }


  // パリティ検査行列を出力する
  // void cLDPC::printParityCheckMatrix()
  // {
  //   std::cout << Hmat;
  //   std::cout << "rank of Hmat is " << Hmat.row_rank() << std::endl;
  // }


  // ガウスの消去法によってPmatを生成
  inline void Ldpc::MakeGeneratorMatrix(itpp::GF2mat &t_hMat, itpp::GF2mat &t_gMat)
  {
  
    itpp::GF2mat obj_hMat = t_hMat;	// ガウスの消去法後のHmat
    perm_.set_size(t_hMat.cols());
    for(int i = 0; i < perm_.size(); ++i){
      perm_[i] = i;
    }

    // /*********ここからは自作のガウスの消去法************/
    
    // rowは行番号でもあり右側行列の対角番号でもある
    // 最右下要素から始める
    for(int index = 0; index < obj_hMat.rows(); index++){
      int row = obj_hMat.rows() - 1 - index;
      int col = obj_hMat.cols() - 1 - index;
      // 単位行列にしたい右側行列の対角要素(row,col)が1になるように
      // 行と列の入れ替えを行なっていく
  
      int i,j;
      for(j = col; j >= 0; j--){
        for(i = row; i >= 0; i--){
          if(obj_hMat(i,j) == 1){
            goto found;
          }
        }
      }
      i = row;
      j = col;
    
    found:
      obj_hMat.swap_rows(row,i);
      obj_hMat.swap_cols(col,j);

      std::swap(perm_[j], perm_[col]);
    
    
      // ここで(row,col)が1になっている
      // (row-1,col)から(0,col)までを0にする

      for(int ref = row-1, m = 1; ref >= 0; ref--, m++){
        // int refCol = Hmat.cols() - col + ref;
        if(obj_hMat(ref,col) == 1){
          obj_hMat.add_rows(ref,row);
        }
      }
      // (row+1,col)から(Hmat.rows()-1,col)を0にする      
      for(int ref = row + 1; ref < obj_hMat.rows(); ref++){
        if(obj_hMat(ref,col) == 1){
          obj_hMat.add_rows(ref,row);
        }
      }
    }

    //   /**********ここまでがガウスの消去法*****************/

    /************ここから余分な行を省く*******************/
    int startRow;
    for(startRow = 0; startRow < obj_hMat.rows() - 1; startRow++){
      for(int i = 0; i < obj_hMat.cols(); i++){
        if(obj_hMat(startRow, i) == 1){
          goto cutting;
        }
      }
    }
  
  cutting:
    obj_hMat = obj_hMat.get_submatrix(startRow,0,obj_hMat.rows()-1,obj_hMat.cols()-1);

    // std::cout << "## Matrix H.\n" << Hmat;
    //   std::cout << "## Matrix P.\n" << Pmat;
 
    // 生成行列を作る
    itpp::GF2mat pMat;

    pMat = obj_hMat.get_submatrix(0,0,obj_hMat.rows()-1,obj_hMat.cols()-obj_hMat.rows() - 1);

    itpp::GF2mat pMatTrans = pMat.transpose();
  
    itpp::GF2mat iMat = itpp::gf2dense_eye(pMatTrans.rows());
  
    t_gMat = iMat.concatenate_horizontal(pMatTrans);
    // std::cout << "## Matrix G.\n" << Gmat;

    t_hMat.permute_cols(perm_,0);
  }

  // G*Ht=0かどうかチェックする
  inline void Ldpc::CheckParityCondition(itpp::GF2mat &t_hMat, itpp::GF2mat &t_gMat)
  {
    itpp::GF2mat MulMat;
  
    MulMat = t_gMat * t_hMat.transpose();
  
    // std::cout << "Matrix MulMat.\n" << MulMat;
    if(!(MulMat.is_zero())){
      std::cerr << "Error: Generator matrix can not be made.\n";
      exit(1);
    }

  }

  // itpp::bvecと1である要素番号が格納されたstd::vector<std::vector< int > >の掛け算
  // ただし、vector<vector <int > >の方は転置した形で格納されていなければならない
  // つまり、bvec * Mat の場合は、bvec * Mat_transposeを入れる
  // 同様に、bvec * Mat_transposeは bvec * Matとする
  inline itpp::bvec Times(const itpp::bvec &in_bvec, const std::vector<std::vector< int > >& in_mat)
  {
    itpp::bvec out(in_mat.size());
    out.zeros();

    int row = 0;
    for(std::vector< std::vector< int > >::const_iterator itRow = in_mat.begin();
        itRow != in_mat.end(); ++itRow, ++row){
      itpp::bin sum = 0;
      for(std::vector< int >::const_iterator itCol = (*itRow).begin(); itCol != (*itRow).end(); ++itCol){
	int nElem = *itCol; // 1である要素番号
	// assert(nElem < in_bvec.size()); // ##
	sum += in_bvec[nElem] * itpp::bin(1);
      }	// for col
      out[row] = sum;
    } // for row

    return out;
  }

  
  // 符号化
  // inputの長さはinfolengthと同じでなければならない
  void Ldpc::Encode(const itpp::bvec &input, itpp::bvec &coded)
  {
    assert(setDone_);

    coded = Times(input, gMatTrans_); // input * Gmat
  
  }
  
  // Gallagerのf関数
  inline double Ldpc::Ffunction(double x)
  {

    if(x == 0){
      x = 0.0000000001;
     }
  
    double nume = exp(x)+1;
    double denom = exp(x)-1;

  
    return log(nume/denom);
  }

  // 要素番号だけ格納されたelemMatに対応するように値が格納されているmatを転置する
  // inline mylib::vec_2D TransposeCompressedVector_2D(const mylib::vec_2D &mat,
  //                                                   const mylib::ivec_2D &elemMat,
  //                                                   const mylib::ivec_2D &elemMat_trans)
  // {
  //   mylib::vec_2D mat_trans(elemMat_trans.size_rows());
  
  //   for(int row = 0; row < elemMat.size_rows(); row++){
  //     for(int col = 0; col < elemMat.size_cols(row); col++){
  //       int nElem = elemMat(row,col);
  //       double value = mat(row,col);
  //       mat_trans.add_cols(nElem, value);
  //     }
  //   }

  //   return mat_trans;
  // }


  // 行処理
  inline void Ldpc::RowsProcessing(itpp::mat* alphaTrans, const itpp::mat &beta, 
                            const itpp::vec &llrVec)
  {
    std::vector< int > colsIndex(alphaTrans->rows(), 0);
    
    for(int m = 0; m < static_cast< int >(hMat_.size()); ++m){
      //    std::cout << "## m = " << m << "\n";
      for(int n = 0; n < static_cast< int >(hMat_[m].size()); ++n){
        //      std::cout << "  ## n = " << n << "\n";
        int product = 1;
        double sum = 0.0;

        for (int i = 0; i < n; ++i){
          product *= itpp::sign(llrVec[hMat_[m][i]] + beta(m,i));
          sum += Ffunction(fabs(llrVec[hMat_[m][i]] + beta(m,i)));
        } // for i
        
        for(int i = n+1; i < static_cast< int >(hMat_[m].size()); i++){
          // std::cout << "    ## i = " << i << "\n";
          product *= itpp::sign(llrVec[hMat_[m][i]] + beta(m,i));
          sum += Ffunction(fabs(llrVec[hMat_[m][i]] + beta(m,i)));
        }
      
        // std::cout << "## sum = " << sum << "\n";
        (*alphaTrans)(hMat_[m][n],colsIndex[hMat_[m][n]]) = product * Ffunction(sum); // hMat_[m][n]は列番号
        ++colsIndex[hMat_[m][n]];
      } // for n
    }   // for m
  }


  // 列処理
  // 一旦
  inline void Ldpc::ColsProcessing(const itpp::mat &alphaTrans,
                                   itpp::mat* beta)
  {
    std::vector< int > rowsIndex(beta->rows(), 0);
    
    for(int n = 0; n < static_cast< int >(hMatTrans_.size()); ++n){
      for(int m = 0; m < static_cast< int >(hMatTrans_[n].size()); ++m){
        double sum = 0.0;
        for (int i = 0; i < m; ++i){
          sum += alphaTrans(n,i);
        } // for i
        
        for(int i = m+1; i < static_cast< int >(hMatTrans_[n].size()); ++i){
          sum += alphaTrans(n,i);
        } // for i
        (*beta)(hMatTrans_[n][m], rowsIndex[hMatTrans_[n][m]]) = sum; // hMatTrans_(n,m)は行番号
        ++rowsIndex[hMatTrans_[n][m]];
      } // for m
    }   // for n
  }

  // 一時推定語を求める
  inline void Ldpc::EstimateCode(const itpp::mat &alphaTrans, 
                          itpp::bvec* decoded,
                          const itpp::vec &llrVec)
  {
    std::vector<double> sum(decoded->size());

    for(int n = 0; n < decoded->size(); ++n){
      double sum = 0.0;
      for(int m = 0; m < static_cast< int >(hMatTrans_[n].size()); ++m){
        sum += alphaTrans(n,m);
      } // for m
    
      if(itpp::sign(llrVec[n] + sum) == 1){
        (*decoded)[n] = 0;
      }
      else{
        (*decoded)[n] = 1;
      }
    } // for n
  }

  // パリティ検査
  inline bool Ldpc::CheckParity(const itpp::bvec &decoded)
  {
    /**********************************************************/
    // ここでやることは
    // mylib::times(decoded, Hmat) --- decoded * Hmat.transpose()
    // の結果がすべて0だったらtrue,そうでなければfalseを返す。
    // これの高速化版
    /**********************************************************/
    bool result = true;
    for(int row = 0; row < static_cast< int >(hMat_.size()); ++row){
      itpp::bin sum = 0;
      for(int col = 0; col < static_cast< int >(hMat_[row].size()); ++col){
        int nElem = hMat_[row][col]; // 1である要素番号
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
  int Ldpc::Decode(const itpp::Modulator_2D &mod,
                   const itpp::cvec &symbol,
                   itpp::bvec &decodedBits, 
                   const double N0,
                   const int loopMax)
  {
    assert(setDone_);

    itpp::bvec estimatedCodes = mod.demodulate_bits(symbol);		// sum-product復号法によって得られる推定語
  
    itpp::mat beta(hMat_.size(), hRowWeight_);

    beta.zeros();			// betaを全て0にする

    // std::cout << "rowsIndex.indexVec.size = " << rowsIndex[0].indexVec.size() << "\n";
    // std::cout << "colsIndex.indexVec.size = " << colsIndex[0].indexVec.size() << "\";
  
    itpp::mat alphaTrans(hMatTrans_.size(), hColWeight_); // alphaだけ転置行列も用意しとく

    int loop = 0;
    if (mod.bits_per_symbol() == 1){ // BPSK専用
      itpp::vec llrVec = CalcLLR(mod, symbol, estimatedCodes, N0);
    
      for(loop = 0; loop < loopMax; loop++){
    
        RowsProcessing(&alphaTrans, beta, llrVec);
        
        ColsProcessing(alphaTrans, &beta);
    
        EstimateCode(alphaTrans, &estimatedCodes, llrVec);

        if(CheckParity(estimatedCodes)){
          break;
        }
      }
    }
    else{
      itpp::vec llrVec;

      for(loop = 0; loop < loopMax; loop++){
    
        llrVec = CalcLLR(mod, symbol, estimatedCodes, N0);    

        RowsProcessing(&alphaTrans, beta, llrVec);
        
        ColsProcessing(alphaTrans, &beta);
    
        EstimateCode(alphaTrans, &estimatedCodes, llrVec);

        if(CheckParity(estimatedCodes)){
          break;
        }

      }
    }
  
    decodedBits = estimatedCodes.left(infoLength_); // 復号語の最初のベクトルが元情報
  
    return loop;
  }

  // 受信器側でpadding bitsの数が分かっているとき
  int Ldpc::DecodeWithPadding0(const itpp::Modulator_2D& mod,
                         const itpp::cvec& symbol,
                         itpp::bvec& decodedBits,
                         double N0,
                         int numPads,
                         int loopMax)
  {
    assert(setDone_);

    itpp::bvec estimatedCodes = mod.demodulate_bits(symbol);		// sum-product復号法によって得られる推定語

    for (int i = 0; i < numPads; ++i){
      estimatedCodes[infoLength_-1-i] = 0;
    } // for i

    itpp::mat beta(hMat_.size(), hRowWeight_);
    
    beta.zeros();			// betaを全て0にする

    // std::cout << "rowsIndex.indexVec.size = " << rowsIndex[0].indexVec.size() << "\n";
    // std::cout << "colsIndex.indexVec.size = " << colsIndex[0].indexVec.size() << "\";
    itpp::mat alphaTrans(hMatTrans_.size(), hColWeight_); // alphaだけ転置行列も用意しとく
    
    int loop = 0;
    if (mod.bits_per_symbol() == 1){ // BPSK専用
      itpp::vec llrVec = CalcLLRWithPads(mod, symbol, estimatedCodes, numPads, N0);
    
      for(loop = 0; loop < loopMax; loop++){
        // std::cout << "## loop = " << loop << "\n" << std::flush;
        RowsProcessing(&alphaTrans, beta, llrVec);
        
        ColsProcessing(alphaTrans, &beta);
    
        EstimateCode(alphaTrans, &estimatedCodes, llrVec);

        if(CheckParity(estimatedCodes)){
          break;
        }
      }
    }
    else{
      itpp::vec llrVec;

      for(loop = 0; loop < loopMax; loop++){
    
        llrVec = CalcLLRWithPads(mod, symbol, estimatedCodes, numPads, N0);    

        RowsProcessing(&alphaTrans, beta, llrVec);
        
        ColsProcessing(alphaTrans, &beta);
    
        EstimateCode(alphaTrans, &estimatedCodes, llrVec);

        if(CheckParity(estimatedCodes)){
          break;
        }

      }
    }
  
    decodedBits = estimatedCodes.left(infoLength_); // 復号語の最初のベクトルが元情報
  
    return loop;
  }



  // lambdaを求める
  itpp::vec Ldpc::CalcLLR(const itpp::Modulator_2D &mod,
                          const itpp::cvec &receivedVec,
                          const itpp::bvec &estimatedVec,
                          double N0)
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
        itpp::cvec symAtk0 = mod.modulate_bits(bitsAtk0);
        itpp::cvec symAtk1 = mod.modulate_bits(bitsAtk1);
        // complex型にする
      
        double nume = exp(-pow(abs(receivedVec[i]-symAtk0[0]),2)/N0);
        double denom = exp(-pow(abs(receivedVec[i]-symAtk1[0]),2)/N0);
            
        llrVec[i*bitsPerSymbol + k] = log(nume/denom);
      }
    }

    return llrVec;
  }

  // numPadsは最後のパディングビットの数
  itpp::vec Ldpc::CalcLLRWithPads(const itpp::Modulator_2D &mod, const itpp::cvec &receivedVec, const itpp::bvec &estimatedVec, int numPads, double N0)
  {
    int bitsPerSymbol = mod.bits_per_symbol();

    int numEffectiveBits = infoLength_ - numPads;
    
    if(receivedVec.size()*bitsPerSymbol != estimatedVec.size()){
      std::cerr << "Error in ldpc::calcLLR.\n"
                << "estimatedVec.size() is not correct.\n";
      exit(1);
    }

    itpp::vec llrVec(estimatedVec.size());
    
    int n = 0;
    for(int i = 0; i < receivedVec.size(); ++i){ 
      // i*bitsPerSymbol目からbitsPerSymbol分取り出す
      itpp::bvec tempBits = estimatedVec.get(i*bitsPerSymbol, i*bitsPerSymbol + bitsPerSymbol -1);
        
      for(int k = 0; k < bitsPerSymbol; ++k, ++n){
        
        if(n >= numEffectiveBits && n < infoLength_){
          llrVec[i*bitsPerSymbol + k] = 100;
        }
        else{
          // tempBitsのkビット目が0と1のものを作る
          itpp::bvec bitsAtk1 = tempBits, bitsAtk0 = tempBits;
          bitsAtk0[k] = 0;
          bitsAtk1[k] = 1;
          // これらをシンボルにする
          // kビット目が0のときと1のときのシンボル(Vecだけど要素数1)
          itpp::cvec symAtk0 = mod.modulate_bits(bitsAtk0);
          itpp::cvec symAtk1 = mod.modulate_bits(bitsAtk1);
          // complex型にする
      
          double nume = exp(-pow(abs(receivedVec[i]-symAtk0[0]),2)/N0);
          double denom = exp(-pow(abs(receivedVec[i]-symAtk1[0]),2)/N0);
            
          llrVec[i*bitsPerSymbol + k] = log(nume/denom);
        }
        
      }
    }

    return llrVec;

  }
  
  // this is for MLC and MSD.
  itpp::vec LdpcForMlcMsd::CalcLLR(const std::vector< itpp::Modulator_2D> &vecMod,
                          const itpp::cvec &receivedVec,
                          double N0)
  {
  
    assert(receivedVec.size() == static_cast<int>(vecMod.size()));

    itpp::vec llrVec(receivedVec.size());

    for(int i = 0; i < static_cast<int>(receivedVec.size()); i++){

      itpp::cvec symbols = vecMod[i].get_symbols();
      itpp::ivec bits2symbols = vecMod[i].get_bits2symbols();
      int nSymbols = symbols.size();

      double nume = 0.0;		// 分子
      double denom = 0.0;		// 分母
      // vecMod[i]の各symbolsの最下位ビットが0か1かで分母と分子の計算を分ける
      for (int s = 0; s < nSymbols/2; ++s){ // MSBが0
        nume += exp(-pow(abs(symbols[bits2symbols[s]] - receivedVec[i]),2)/N0);
      } // for s
      for (int s = nSymbols/2; s < nSymbols; ++s){ // MSBが1
        denom += exp(-pow(abs(symbols[bits2symbols[s]] - receivedVec[i]),2)/N0);
      } // for s
    
      llrVec[i] = log(nume/denom);

    } // for i
    
    return llrVec;
  }

  // This is for multithread MSD.
  itpp::vec LdpcForMlcMsd::CalcLLRAtLevel(const itpp::Modulator_2D &mod,
                          const itpp::cvec &receivedVec,
                          double N0,
                          int level)
  {
    int bitsPerSymbol = mod.bits_per_symbol();
    assert(level >= 0 && level < mod.bits_per_symbol());

    itpp::cvec symbols = mod.get_symbols();
    itpp::ivec bits2symbols = mod.get_bits2symbols();
    int nSymbols = symbols.size();

    // 最初にlevelの位置のビットが0か1かの違いで数字を分けておく
    std::vector< int > figure0(nSymbols/2), figure1(nSymbols/2);
    for (int i = 0, j0 = 0, j1 = 0; i < nSymbols; ++i){
      itpp::bvec bin = itpp::dec2bin(bitsPerSymbol, i);
      if (bin[level] == 0){
        figure0[j0] = i;
        ++j0;
      } // if
      else{
        figure1[j1] = i;
        ++j1;
      } // else
    } // for i
    
    itpp::vec llrVec(receivedVec.size());

    for(int i = 0; i < static_cast<int>(receivedVec.size()); i++){

      double nume = 0.0;		// 分子
      double denom = 0.0;		// 分母
      // vecMod[i]の各symbolsの最下位ビットが0か1かで分母と分子の計算を分ける
      for (int s = 0; s < nSymbols/2; ++s){ 
        nume += exp(-pow(abs(symbols[bits2symbols[figure0[s]]] - receivedVec[i]),2)/N0);
        denom += exp(-pow(abs(symbols[bits2symbols[figure1[s]]] - receivedVec[i]),2)/N0);
      }
      llrVec[i] = log(nume/denom);

    } // for i
    
    return llrVec;
  }

  
  // this is for MLC and MSD.
  int LdpcForMlcMsd::Decode(const std::vector< itpp::Modulator_2D > &vecMod,
                   const itpp::cvec &symbol,
                   itpp::bvec &decodedCodes, 
                   const double N0,
                   const int loopMax)
  {
    assert(setDone_);
    
    itpp::bvec estimatedCodes(symbol.size()); // sum-product復号法によって得られる推定語

    itpp::mat beta(hMat_.size(), hRowWeight_);

    beta.zeros();			// betaを全て0にする  

    itpp::mat alphaTrans(hMatTrans_.size(), hColWeight_); // alphaだけ転置行列も用意しとく
    
    itpp::vec llrVec = CalcLLR(vecMod, symbol, N0);    
  
    int loop;
    for(loop = 0; loop < loopMax; loop++){
      // std::cout << "## loop = " << loop << "\r" << std::flush;

      // alpha.zeros();		// 一応初期化
      RowsProcessing(&alphaTrans, beta, llrVec);
      //    std::cout << "## rows process done.\n";
    
      ColsProcessing(alphaTrans, &beta);
      //    std::cout << "## cols process done.\n";


      EstimateCode(alphaTrans, &estimatedCodes, llrVec);

      if(CheckParity(estimatedCodes)){
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

  // This is for multithread MSD.
  DecoderParas LdpcForMlcMsd::DecodeAtLevel(const itpp::Modulator_2D& mod,
                                   const itpp::cvec &symbol,
                                   // itpp::bvec &decodedCodes,
                                   double N0,
                                   int level,
                                   // int *loop,
                                   int loopMax)
  {
    assert(setDone_);

    itpp::bvec decodedCodes(symbol.size());
    
    // assert(hMat_.is_rectangular() && hMatTrans_.is_rectangular());

    itpp::mat beta(hMat_.size(), hRowWeight_);

    beta.zeros();			// betaを全て0にする
    itpp::mat alphaTrans(hMatTrans_.size(), hColWeight_); // alphaだけ転置行列も用意しとく
    
    
    itpp::vec llrVec = CalcLLRAtLevel(mod, symbol, N0, level);

    // std::cout << "## DecodeAtLevel " << level << " is called." << std::endl;
    int loop;
    for(loop = 0; loop < loopMax; ++loop){
      // std::cout << "## level, loop = " << level << " " << *loop << std::endl;
      RowsProcessing(&alphaTrans, beta, llrVec);
      ColsProcessing(alphaTrans, &beta);
      EstimateCode(alphaTrans, &decodedCodes, llrVec);

      if(CheckParity(decodedCodes)){
        break;
      }
    }

    DecoderParas decoderParas;
    decoderParas.decodedBits = decodedCodes;
    decoderParas.loop = loop;

    return decoderParas;
    // std::cout << "## decodedCodes.size() == " << decodedCodes.size() << std::endl;
    
  }

}
