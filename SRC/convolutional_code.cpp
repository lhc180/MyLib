#include "../include/convolutional_code.h"
#include "../include/myutl.h"
#include <cassert>


void cConvolutionalCode::checkInitialization()
{
  if(init != true){
    std::cout << "Error: must initialize the convolutional code cless.\n";
    exit(1);
  }
}

cConvolutionalCode::cConvolutionalCode(const int constraint, const int invRate = 2)
  : constraintLength(constraint),
    inverseRate(invRate),
    registerSize(constraint - 1),
    stateNum(static_cast<int>(pow(2,constraint - 1)))
{
  initialize();
}

void cConvolutionalCode::setCode(const int constraint, const int invRate = 2)
{
  inverseRate = invRate;
  constraintLength = constraint;
  registerSize = constraint - 1;
  stateNum = static_cast<int>(pow(2,constraint - 1));

  initialize();
    
}


void cConvolutionalCode::initialize()
{
  generatorPolynomial.resize(inverseRate);

  
  // 生成多項式を8進表示で決定する
  if(inverseRate == 2){
    if(registerSize == 2){
      generatorPolynomial[0] = 05;
      generatorPolynomial[1] = 07;
    }
    else if(registerSize == 3){
      generatorPolynomial[0] = 015;
      generatorPolynomial[1] = 017;
    }
    else if(registerSize == 4){
      generatorPolynomial[0] = 023;
      generatorPolynomial[1] = 035;
    }
    else if(registerSize == 5){
      generatorPolynomial[0] = 053;
      generatorPolynomial[1] = 075;
    }
    else if(registerSize == 6){
      generatorPolynomial[0] = 0133;
      generatorPolynomial[1] = 0171;
    }
    else if(registerSize == 7){
      generatorPolynomial[0] = 0247;
      generatorPolynomial[1] = 0371;
    }
    else if(registerSize == 8){
      generatorPolynomial[0] = 0561;
      generatorPolynomial[1] = 0753;
    }
    else if(registerSize == 9){
      generatorPolynomial[0] = 01167;
      generatorPolynomial[1] = 01545;
    }
    else if(registerSize == 10){
      generatorPolynomial[0] = 02335;
      generatorPolynomial[1] = 03661;
    }
    else if(registerSize == 11){
      generatorPolynomial[0] = 04335;
      generatorPolynomial[1] = 05723;
    }
    else if(registerSize == 12){
      generatorPolynomial[0] = 010533;
      generatorPolynomial[1] = 017661;
    }
    else if(registerSize == 13){
      generatorPolynomial[0] = 021675;
      generatorPolynomial[1] = 027123;
    }
    else{
      std::cout << "Error: must choose available value for constraints length.\n";
      exit(1);
    }
  }
  else{
    std::cout << "Error: must choose avilable value for coding rate.\n";
    exit(1);
  }

  en_shiftRegister = 0;

  termination = false;

  init = true;
}

void cConvolutionalCode::encode(const std::vector<int> &input, std::vector<int> &output)
{
  checkInitialization();	// 初期化してあるかチェック

  unsigned shiftRegister = 0; // シフトレジスタの中身はビット列で表す  

  std::vector<int>::iterator it;
  for(int n = 0; n < static_cast<int>(input.size()); n++){
    if(!((input[n] == 0) || (input[n] == 1))){
      std::cout << "Error: illegal input is used in cConvolutionalCode::encode.\n";
      std::cout << "n = " << n << std::endl;
      exit(1);
    }
    else{
      // ◯-□-□-...-□ (◯は入力ビット、□がレジスタ)という状態
      unsigned tempRegister = (input[n] << registerSize) | shiftRegister;
      for(int i = 0; i < inverseRate; i++){
	unsigned tempOutput = 0;
	// レジスタと生成多項式のorが出力の計算に使われる
	unsigned tempCalc = tempRegister & generatorPolynomial[i];
	
	for(int j = 0; j < constraintLength; j++){
	  // 出力はtempCalcの排他的論理和
	  tempOutput ^= (tempCalc >> j) & 1;
	}
	output.push_back(static_cast<int>(tempOutput));
      }
      shiftRegister = tempRegister >> 1;
    }
  }
}

void cConvolutionalCode::encode_term(const std::vector<int> &input, std::vector<int> &output)
{
  std::vector<int> copyInput(input.size());

  for(int i = 0; i < static_cast<int>(input.size()); i++){
    copyInput[i] = input[i];
  }
  
  // 終端ビット
  for(int i = 0; i < registerSize; i++){
    copyInput.push_back(0);
  }

  termination = true;

  encode(copyInput, output);
}

void cConvolutionalCode::encode(const itpp::bvec &input, itpp::bvec &output)
{
  checkInitialization();	// 初期化してあるかチェック

  output.set_size(0);
  
  for(int n = 0; n < input.size(); n++){

    // ■-□-□-...-□ (■は入力ビット、□がレジスタ)という状態
    unsigned tempRegister = (input[n] << registerSize) | en_shiftRegister;
    for(int i = 0; i < inverseRate; i++){
      unsigned tempOutput = 0;
      // レジスタと生成多項式のandが出力の計算に使われる
      unsigned tempCalc = tempRegister & generatorPolynomial[i];
	
      for(int j = 0; j < constraintLength; j++){
        // 出力はtempCalcの排他的論理和
        tempOutput ^= (tempCalc >> j) & 1;
      }
      output = itpp::concat(output, itpp::bin(tempOutput));
    }
    en_shiftRegister = tempRegister >> 1;

  }
}

void cConvolutionalCode::encode_term(const itpp::bvec &input, itpp::bvec &output)
{
  itpp::bvec copyInput = input;

  // 終端ビット
  for(int i = 0; i < registerSize; i++){
    copyInput = itpp::concat(copyInput, itpp::bin(0));
  }

  termination = true;

  encode(copyInput, output);
}
	
int cConvolutionalCode::getNextState(int lastState, int bit)
{
  unsigned tempRegister = (bit << registerSize) | lastState;

  return static_cast<int>(tempRegister >> 1);
}
  
inline int calcHammingDistance(const std::vector<int> &sequence1, const std::vector<int> &sequence2)
{
  int hamming = 0;
  for(int i = 0; i < static_cast<int>(sequence1.size()); i++){
    // パンクチャリングされていたとき
    if(sequence1[i] < 0){
      continue;
    }
    if(sequence1[i] != sequence2[i]){
      hamming++;
    }
  }

  return hamming;
}
  
inline double bpsk_mod(int input)
{
  if(input == 0)
    return 1;
  else
    return -1;
}


inline double calcEuclideanDistance(const std::vector<double> &sequence1, const std::vector<int> &sequence2)
{
  double euclidean = 0.0;

  for(int i = 0; i < static_cast<int>(sequence1.size()); i++){
    euclidean += sequence1[i] * bpsk_mod(sequence2[i]);
  }

  return euclidean;
}
    

// この中身はencodeとほとんど一緒なので、将来的にくっつけてもいい
void cConvolutionalCode::getEncodedSequence(int currentState, int input, std::vector<int> &sequence)
{
  unsigned tempRegister = (input << registerSize) | currentState;
  for(int i = 0; i < inverseRate; i++){
    unsigned tempOutput = 0;
    // レジスタと生成多項式のorが出力の計算に使われる
    unsigned tempCalc = tempRegister & generatorPolynomial[i];
    
    for(int j = 0; j < constraintLength; j++){
      // 出力はtempCalcの排他的論理和
      tempOutput ^= (tempCalc >> j) & 1;
    }
    sequence.push_back(static_cast<int>(tempOutput));
  }
}
  
int cConvolutionalCode::get1bit(int currentState)
{
  int out;
  
  out = currentState >> (registerSize - 1);

  return out;
}
  
void cConvolutionalCode::decode_hd(const std::vector<int> &input, std::vector<int> &output)
{
  // pathMetric[n*stateNum + j]はn番目のトレリスのj状態を表す
  std::vector<int> pathMetric(stateNum * (input.size() / inverseRate + 1));
  std::vector<int> stateFlow(stateNum * (input.size() / inverseRate + 1));

  output.resize(input.size()/inverseRate);

  // パスメトリックを初期化
  for(int i = 0; i < static_cast<int>(pathMetric.size()); i++){
    pathMetric[i] = input.size();
    stateFlow[i] = -1;
  }

  pathMetric[0] = 0;

  // パスメトリックの計算
  for(int n = 1; n <= static_cast<int>(input.size()/inverseRate); n++){
    for(int lastState = 0; lastState < stateNum; lastState++){
      for(int bit = 0; bit <= 1; bit++){
	int currentState = getNextState(lastState, bit); 
	if(currentState >= stateNum){
	  std::cout << "Error: illegal state.\n";
	}
  
	std::vector<int> tempInput;
	for(int i = 0; i < inverseRate; i++){
	  tempInput.push_back(input[inverseRate*(n-1) + i]);
	}
	
	std::vector<int> tempSequence;
	getEncodedSequence(lastState, bit, tempSequence);
	
	int tempPathMetric = pathMetric[(n-1)*stateNum + lastState] + calcHammingDistance(tempInput, tempSequence);
	
	// 時刻nの状態currentStateのpathMetricを値を小さい方を採用
	if(tempPathMetric < pathMetric[n*stateNum + currentState]){
	  pathMetric[n*stateNum + currentState] = tempPathMetric;
	  stateFlow[n*stateNum + currentState] = lastState;
	}
      }
    }
  }

  // トレースバック
  int currentState = 0;
  // 終端処理してある場合は0のステートからスタート
  if(termination){
    currentState = 0;
  }
  else{
    int metric = input.size();
    for(int i = 0; i < stateNum; i++){
      if(pathMetric[input.size()/inverseRate + i] < metric){
	currentState = i;
      }
    }
  }
      
  for(int n = input.size()/inverseRate; n > 0; n--){
    int lastState = stateFlow[n*stateNum + currentState];
    if(lastState < 0 || lastState >= stateNum){
      std::cout << "Error: state is illegal in trace back.\n";
      exit(1);
    }
    output[n-1] = get1bit(currentState);
	  
    currentState = lastState;
  }
}

void cConvolutionalCode::decode_sd(const std::vector<double> &input,
					   std::vector<int> &output)
{
  // pathMetric[n*stateNum + j]はn番目のトレリスのj状態を表す
  std::vector<double> pathMetric(stateNum * (input.size() / inverseRate + 1));
  std::vector<int> stateFlow(stateNum * (input.size() / inverseRate + 1));

  output.resize(input.size()/inverseRate);

  // パスメトリックを初期化
  for(int i = 0; i < static_cast<int>(pathMetric.size()); i++){
    pathMetric[i] = -mylib::INFTY;
    stateFlow[i] = -1;
  }

  pathMetric[0] = 0;

  // パスメトリックの計算
  for(int n = 1; n <= static_cast<int>(input.size()/inverseRate); n++){
    for(int lastState = 0; lastState < stateNum; lastState++){
      for(int bit = 0; bit <= 1; bit++){
	int currentState = getNextState(lastState, bit); 
	if(currentState >= stateNum){
	  std::cout << "Error: illegal state.\n";
	}
  
	std::vector<double> tempInput;
	for(int i = 0; i < inverseRate; i++){
	  tempInput.push_back(input[inverseRate*(n-1) + i]);
	}
	
	std::vector<int> tempSequence;
	getEncodedSequence(lastState, bit, tempSequence);
	
	double tempPathMetric = pathMetric[(n-1)*stateNum + lastState] + calcEuclideanDistance(tempInput, tempSequence);
	
	// 時刻nの状態currentStateのpathMetricを値を大きい方を採用
	if(tempPathMetric >= pathMetric[n*stateNum + currentState]){
	  pathMetric[n*stateNum + currentState] = tempPathMetric;
	  stateFlow[n*stateNum + currentState] = lastState;
	}
      }
    }
  }

  // トレースバック
  int currentState = 0;
  // 終端処理してある場合は0のステートからスタート
  if(termination){
    currentState = 0;
  }
  else{
    int metric = input.size();
    for(int i = 0; i < stateNum; i++){
      if(pathMetric[input.size()/inverseRate + i] < metric){
	currentState = i;
      }
    }
  }
      
  for(int n = input.size()/inverseRate; n > 0; n--){
    int lastState = stateFlow[n*stateNum + currentState];
    if(lastState < 0 || lastState >= stateNum){
      std::cout << "Error: state is illegal in trace back.\n";
      exit(1);
    }
    output[n-1] = get1bit(currentState);
	  
    currentState = lastState;
  }
}

void cConvolutionalCode::decode_sd(const itpp::vec &input,
					   itpp::bvec &output)
{
  // pathMetric[n*stateNum + j]はn番目のトレリスのj状態を表す
  std::vector<double> pathMetric(stateNum * (input.size() / inverseRate + 1));
  
  std::vector<int> stateFlow(stateNum * (input.size() / inverseRate + 1));
   
  output.set_size(input.size()/inverseRate);

  // パスメトリックを初期化
  for(int i = 0; i < static_cast<int>(pathMetric.size()); i++){
    pathMetric[i] = -mylib::INFTY;
    stateFlow[i] = -1;
  }

  pathMetric[0] = 0;

  // パスメトリックの計算
  for(int n = 1; n <= static_cast<int>(input.size()/inverseRate); n++){
    for(int lastState = 0; lastState < stateNum; lastState++){
      for(int bit = 0; bit <= 1; bit++){
	int currentState = getNextState(lastState, bit); 
	if(currentState >= stateNum){
	  std::cout << "Error: illegal state.\n";
	}
  
	std::vector<double> tempInput;
	for(int i = 0; i < inverseRate; i++){
	  tempInput.push_back(input[inverseRate*(n-1) + i]);
	}
	
	std::vector<int> tempSequence;
	getEncodedSequence(lastState, bit, tempSequence);
	
	double tempPathMetric = pathMetric[(n-1)*stateNum + lastState] + calcEuclideanDistance(tempInput, tempSequence);
	
	// 時刻nの状態currentStateのpathMetricを値を大きい方を採用
	if(tempPathMetric >= pathMetric[n*stateNum + currentState]){
	  pathMetric[n*stateNum + currentState] = tempPathMetric;
	  stateFlow[n*stateNum + currentState] = lastState;
	}
      }
    }
  }

  // トレースバック
  int currentState = 0;
  // 終端処理してある場合は0のステートからスタート
  if(termination){
    currentState = 0;
  }
  else{
    int metric = input.size();
    for(int i = 0; i < stateNum; i++){
      if(pathMetric[input.size()/inverseRate + i] < metric){
	currentState = i;
      }
    }
  }
      
  for(int n = input.size()/inverseRate; n > 0; n--){
    int lastState = stateFlow[n*stateNum + currentState];
    if(lastState < 0 || lastState >= stateNum){
      std::cout << "Error: state is illegal in trace back.\n";
      exit(1);
    }
    output[n-1] = get1bit(currentState);
	  
    currentState = lastState;
  }

  if(termination){
    output = output.left(output.size()-registerSize);
  }
}



// ++++++++++++++++++++ cPuncturedConvolutionalCode +++++++++++++++++
cPuncturedConvolutionalCode::cPuncturedConvolutionalCode(const int constraint,
                                                                const int invRate,
                                                                const int pp,
                                                                const int denom)
{
  setCode(constraint, invRate, pp, denom);
 
}

void cPuncturedConvolutionalCode::setCode(const int constraint,
                                          const int invRate,
                                          const int pp,
                                          const int denom)
{
  cConvolutionalCode::setCode(constraint, invRate);

  if(invRate*pp == denom){
    puncturingPeriod = 1;
    denominatorCodingRate = invRate;
  }
  else{
    puncturingPeriod = pp;
    denominatorCodingRate = denom;
  }
  setPuncturingMatrix();
}			       

// パンクチャリング行列を設定する
void cPuncturedConvolutionalCode::setPuncturingMatrix()
{

  puncturingMatrix.resize(0);
  if(inverseRate == 2){
    for(int i = 0; i < puncturingPeriod * inverseRate; i++){
      puncturingMatrix.push_back(1);
    }
    if(puncturingPeriod == 1){
      assert(denominatorCodingRate == 2); // 
    }
    else{
      if(registerSize == 2){
        if(puncturingPeriod == 2){
          if(denominatorCodingRate == 3){
            puncturingMatrix[2] = 0;
          }
          else{
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 3){
          switch(denominatorCodingRate){
          case 4:
            puncturingMatrix[3] = 0;
          case 5:
            puncturingMatrix[4] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 4){
          switch(denominatorCodingRate){
          case 5:
            puncturingMatrix[7] = 0;
          case 6:
            puncturingMatrix[2] = 0;
          case 7:
            puncturingMatrix[4] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 5){
          switch(denominatorCodingRate){
          case 6:
            puncturingMatrix[2] = 0;
          case 7:
            puncturingMatrix[0] = 0;
          case 8:
            puncturingMatrix[6] = 0;
          case 9:
            puncturingMatrix[4] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }	  
        }
        else if(puncturingPeriod == 6){
          switch(denominatorCodingRate){
          case 7:
            puncturingMatrix[5] = 0;
          case 8:
            puncturingMatrix[9] = 0;
          case 9:
            puncturingMatrix[3] = 0;
          case 10:
            puncturingMatrix[0] = 0;
          case 11:
            puncturingMatrix[6] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 7){
          switch(denominatorCodingRate){
          case 8:
            puncturingMatrix[3] = 0;
          case 9:
            puncturingMatrix[9] = 0;
          case 10:
            puncturingMatrix[10] = 0;
          case 11:
            puncturingMatrix[4] = 0;
          case 12:
            puncturingMatrix[12] = 0;
          case 13:
            puncturingMatrix[6] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else{
          std::cerr << "Error: must choose avilable puncturing period.\n";
          exit(1);
        }
      }
      else if(registerSize == 3){
        if(puncturingPeriod == 2){
          if(denominatorCodingRate == 3){
            puncturingMatrix[1] = 0;
          }
          else{
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 3){
          switch(denominatorCodingRate){
          case 4:
            puncturingMatrix[5] = 0;
          case 5:
            puncturingMatrix[0] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 4){
          switch(denominatorCodingRate){
          case 5:
            puncturingMatrix[3] = 0;
          case 6:
            puncturingMatrix[6] = 0;
          case 7:
            puncturingMatrix[5] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 5){
          switch(denominatorCodingRate){
          case 6:
            puncturingMatrix[4] = 0;
          case 7:
            puncturingMatrix[8] = 0;
          case 8:
            puncturingMatrix[2] = 0;
          case 9:
            puncturingMatrix[1] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }	  
        }
        else if(puncturingPeriod == 6){
          switch(denominatorCodingRate){
          case 7:
            puncturingMatrix[7] = 0;
          case 8:
            puncturingMatrix[4] = 0;
          case 9:
            puncturingMatrix[3] = 0;
          case 10:
            puncturingMatrix[10] = 0;
          case 11:
            puncturingMatrix[1] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 7){
          switch(denominatorCodingRate){
          case 8:
            puncturingMatrix[4] = 0;
          case 9:
            puncturingMatrix[6] = 0;
          case 10:
            puncturingMatrix[9] = 0;
          case 11:
            puncturingMatrix[1] = 0;
          case 12:
            puncturingMatrix[13] = 0;
          case 13:
            puncturingMatrix[11] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else{
          std::cerr << "Error: must choose avilable puncturing period.\n";
          exit(1);
        }
      }
      else if(registerSize == 4){
        if(puncturingPeriod == 2){
          if(denominatorCodingRate == 3){
            puncturingMatrix[1] = 0;
          }
          else{
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 3){
          switch(denominatorCodingRate){
          case 4:
            puncturingMatrix[2] = 0;
          case 5:
            puncturingMatrix[3] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 4){
          switch(denominatorCodingRate){
          case 5:
            puncturingMatrix[7] = 0;
          case 6:
            puncturingMatrix[1] = 0;
          case 7:
            puncturingMatrix[5] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 5){
          switch(denominatorCodingRate){
          case 6:
            puncturingMatrix[3] = 0;
          case 7:
            puncturingMatrix[8] = 0;
          case 8:
            puncturingMatrix[1] = 0;
          case 9:
            puncturingMatrix[5] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }	  
        }
        else if(puncturingPeriod == 6){
          switch(denominatorCodingRate){
          case 7:
            puncturingMatrix[0] = 0;
          case 8:
            puncturingMatrix[9] = 0;
          case 9:
            puncturingMatrix[2] = 0;
          case 10:
            puncturingMatrix[11] = 0;
          case 11:
            puncturingMatrix[5] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 7){
          switch(denominatorCodingRate){
          case 8:
            puncturingMatrix[10] = 0;
          case 9:
            puncturingMatrix[7] = 0;
          case 10:
            puncturingMatrix[12] = 0;
          case 11:
            puncturingMatrix[1] = 0;
          case 12:
            puncturingMatrix[9] = 0;
          case 13:
            puncturingMatrix[5] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else{
          std::cerr << "Error: must choose avilable puncturing period.\n";
          exit(1);
        }
      }
      else if(registerSize == 5){
        if(puncturingPeriod == 2){
          if(denominatorCodingRate == 3){
            puncturingMatrix[2] = 0;
          }
          else{
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 3){
          switch(denominatorCodingRate){
          case 4:
            puncturingMatrix[4] = 0;
          case 5:
            puncturingMatrix[2] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 4){
          switch(denominatorCodingRate){
          case 5:
            puncturingMatrix[4] = 0;
          case 6:
            puncturingMatrix[6] = 0;
          case 7:
            puncturingMatrix[2] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 5){
          switch(denominatorCodingRate){
          case 6:
            puncturingMatrix[2] = 0;
          case 7:
            puncturingMatrix[6] = 0;
          case 8:
            puncturingMatrix[1] = 0;
          case 9:
            puncturingMatrix[4] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }	  
        }
        else if(puncturingPeriod == 6){
          switch(denominatorCodingRate){
          case 7:
            puncturingMatrix[5] = 0;
          case 8:
            puncturingMatrix[1] = 0;
          case 9:
            puncturingMatrix[7] = 0;
          case 10:
            puncturingMatrix[3] = 0;
          case 11:
            puncturingMatrix[9] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 7){
          switch(denominatorCodingRate){
          case 8:
            puncturingMatrix[13] = 0;
          case 9:
            puncturingMatrix[7] = 0;
          case 10:
            puncturingMatrix[10] = 0;
          case 11:
            puncturingMatrix[1] = 0;
          case 12:
            puncturingMatrix[5] = 0;
          case 13:
            puncturingMatrix[9] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else{
          std::cerr << "Error: must choose avilable puncturing period.\n";
          exit(1);
        }
      }
      else if(registerSize == 6){
        if(puncturingPeriod == 2){
          if(denominatorCodingRate == 3){
            puncturingMatrix[3] = 0;
          }
          else{
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 3){
          switch(denominatorCodingRate){
          case 4:
            puncturingMatrix[4] = 0;
          case 5:
            puncturingMatrix[3] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 4){
          switch(denominatorCodingRate){
          case 5:
            puncturingMatrix[5] = 0;
          case 6:
            puncturingMatrix[7] = 0;
          case 7:
            puncturingMatrix[3] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 5){
          switch(denominatorCodingRate){
          case 6:
            puncturingMatrix[9] = 0;
          case 7:
            puncturingMatrix[5] = 0;
          case 8:
            puncturingMatrix[3] = 0;
          case 9:
            puncturingMatrix[7] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }	  
        }
        else if(puncturingPeriod == 6){
          switch(denominatorCodingRate){
          case 7:
            puncturingMatrix[11] = 0;
          case 8:
            puncturingMatrix[8] = 0;
          case 9:
            puncturingMatrix[7] = 0;
          case 10:
            puncturingMatrix[0] = 0;
          case 11:
            puncturingMatrix[5] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else if(puncturingPeriod == 7){
          switch(denominatorCodingRate){
          case 8:
            puncturingMatrix[9] = 0;
          case 9:
            puncturingMatrix[11] = 0;
          case 10:
            puncturingMatrix[5] = 0;
          case 11:
            puncturingMatrix[3] = 0;
          case 12:
            puncturingMatrix[13] = 0;
          case 13:
            puncturingMatrix[7] = 0;
            break;
          default:
            std::cerr << "Error: must choose available denominator of coding rate.\n";
            exit(1);
          }
        }
        else{
          std::cerr << "Error: must choose avilable puncturing period.\n";
          exit(1);
        }
      }
      else{
        std::cerr << "Error: must choose available register size.\n";
        exit(1);
      }
    } // else of if(puncturingPeriod == 1)
  }
  else{
    std::cerr << "Error: must choose available coding rate.\n";
    exit(1);
  }
}

itpp::bvec cPuncturedConvolutionalCode::puncturing(const itpp::bvec &input)
{
  //  std::cout << "## puncturing called.\n";

  itpp::bvec output(0);

  int p = 0;
  for(int i = 0; i < input.size(); i++){
    if(puncturingMatrix[p] == 1){
      output = itpp::concat(output, input[i]);
    }
    p++;
    if(p == inverseRate * puncturingPeriod){
      p = 0;
    }
  }
      
  return output;
}

itpp::vec cPuncturedConvolutionalCode::fillingUp(const itpp::vec &input)
{
  int unitOfDenom = input.size()/denominatorCodingRate;
  int firstLength = unitOfDenom * denominatorCodingRate;
  int surplusLength = input.size() - firstLength; // denomで割り切れない分

  assert(surplusLength < denominatorCodingRate);
  
  itpp::vec output(unitOfDenom * puncturingPeriod * inverseRate);

  // まず割り切れる部分だけダミーシンボルを挿入していく
  int p = 0;
  int j = 0;
  for(int i = 0; i < output.size(); i++){
    if(puncturingMatrix[p] == 1){
      output[i] = input[j];
      j++;
    }
    else{
      output[i] = 0;
    }
    p++;
    if(p == static_cast<int>(puncturingMatrix.size())){
      p = 0;
    }
  }
  
  // ここから余りの部分にダミーシンボルを入れる
  int count = 0;
  p = 0;
  while(count < surplusLength){
    for(int row = 0; row < inverseRate; row++){
      if(puncturingMatrix[p] == 1){
        count++;
        output = itpp::concat(output, input[j]);
        j++;
      }
      else{
        output = itpp::concat(output, 0.0);
      }
      p++;
      assert(p <= static_cast<int>(puncturingMatrix.size()));
    }
  }

  return output;

}
  
void cPuncturedConvolutionalCode::encode(const itpp::bvec &input, itpp::bvec &output)
{
  // std::cout << "## PCC::encode is called.\n";

  itpp::bvec temp;

  cConvolutionalCode::encode(input, temp);
  
  output = puncturing(temp);
  
}

void cPuncturedConvolutionalCode::encode_term(const itpp::bvec &input, itpp::bvec &output)
{
  //  std::cout << "## PCC::encode_term is called.\n";
  itpp::bvec temp;
  
  cConvolutionalCode::encode_term(input, temp);

  //  std::cout << "## temp.size() = " << temp.size() << std::endl;

  output = puncturing(temp);
}
  
void cPuncturedConvolutionalCode::decode_sd(const itpp::vec &input, itpp::bvec &output)
{
  itpp::vec filled = fillingUp(input);
 

  cConvolutionalCode::decode_sd(filled, output);

}

itpp::bvec cRCPCCode::encode(const std::vector<itpp::bvec> &input)
{
  assert(input.size() == denomRate_.size());
  
  cPuncturedConvolutionalCode PCC;

  itpp::bvec output(0);
  unsigned shiftRegister = 0;
  for(int i = 0; i < static_cast<int>(denomRate_.size()); i++){
    itpp::bvec temp;
    PCC.setCode(constraint_, 2, puncturingPeriod_, denomRate_[i]);
    PCC.en_shiftRegister = shiftRegister;
    if(i == static_cast<int>(denomRate_.size())-1){
      PCC.encode_term(input[i], temp);
    }
    else{
      PCC.encode(input[i], temp);
    }
    codeLength_[i] = temp.size();
    shiftRegister = PCC.en_shiftRegister;

    output = itpp::concat(output, temp);
  }

   return output;
}

std::vector<itpp::bvec> cRCPCCode::decode(const itpp::vec &input)
{
  int totalCodeLength = 0;
  for(int i = 0; i < static_cast<int>(denomRate_.size()); i++){
    totalCodeLength += codeLength_[i];
  }

  assert(totalCodeLength == input.size());

  std::vector<itpp::bvec> output(denomRate_.size());

  cPuncturedConvolutionalCode PCC;
  
  itpp::vec depunctured(0);
  int startIndex = 0;
  std::vector<int> inputLength(denomRate_.size());
  for(int i = 0; i < static_cast<int>(denomRate_.size()); i++){
    itpp::vec temp = input.mid(startIndex, codeLength_[i]);
    startIndex += codeLength_[i];
    
    itpp::vec t_depunctured;

    PCC.setCode(constraint_, 2, puncturingPeriod_, denomRate_[i]);    

    t_depunctured = PCC.fillingUp(temp);

    depunctured = itpp::concat(depunctured, t_depunctured);
    inputLength[i] = t_depunctured.size()/2;
  }

  cConvolutionalCode CC(constraint_, 2);
  CC.termination = true;

  itpp::bvec temp;
  CC.decode_sd(depunctured, temp);
  
  startIndex = 0;
  for(int i = 0; i < static_cast<int>(denomRate_.size()); i++){
    if(i == static_cast<int>(denomRate_.size())-1){
      output[i] = temp.get(startIndex, -1);
    }
    else{
      output[i] = temp.mid(startIndex, inputLength[i]);
    }
    startIndex += inputLength[i];
  }

  return output;
}

