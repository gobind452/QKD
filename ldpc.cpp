#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include <iterator>
#include <set>
#include <cmath>
#include <unordered_map>
#include <cstdlib>
#define MAX 1000000000000000000
using namespace std;

default_random_engine gen(0);
uniform_real_distribution<float> dist(0.0,1.0);

class Matrix{
public:
	int rows;
	int columns;
	vector<vector<int>> values;
	Matrix(int rows=0,int columns=0):rows(rows),columns(columns){
		for(int i=0;i<columns;i++){
			values.push_back(vector<int>());
		}
	}
	void generateRandomMatrix(float density){
		for(int i=0;i<columns;i++){
			for(int j=0;j<rows;j++){
				float rand = dist(gen);
				if(rand<density){
					this->values[i].push_back(j);
				}
			}
		}
	}
	Matrix getTranspose(){
		Matrix temp = Matrix(columns,rows);
		for(int i=0;i<rows;i++){
			temp.values.push_back(vector<int>());
		}
		for(int i=0;i<columns;i++){
			for(int j=0;j<this->values[i].size();j++){
				temp.values[this->values[i][j]].push_back(i);
			}
		}
		return temp;
	}
	void printMatrix(){
		for(int i=0;i<columns;i++){
			for(int j=0;j<this->values[i].size();j++){
				cout << i << " " << this->values[i][j] << endl;
			}
		}
	}

	void print_full_matrix() {
	  int temp1[rows][columns];
	  for(int ii=0;ii<rows;ii++){
	    for(int jj=0;jj<columns;jj++){
	      temp1[ii][jj] = 0;
	    }
	  }
	  for(int i=0;i<columns;i++){
			for(int j=0;j<this->values[i].size();j++){
				temp1[this->values[i][j]][i] = 1;
			}
		}
		for(int ii=0;ii<rows;ii++){
	    for(int jj=0;jj<columns;jj++){
	      cout << temp1[ii][jj] << ' ';
	    }
	    cout << endl;
	  }
	}

};

void multiply(vector<int>& key,vector<int>& toCalculate,Matrix checkMatrix){
	for(int i=0;i<key.size();i++){
		if(key[i] == 1){
			for(int j=0;j<checkMatrix.values[i].size();j++){
				toCalculate[checkMatrix.values[i][j]] = toCalculate[checkMatrix.values[i][j]] ^ 1;
			}
		}
	}
	return;
}

uniform_int_distribution<unsigned long long int> distribution(0,MAX);

Matrix Gallager(int n, int wc, int wr) {
	int m = (n/wr)*wc;
	int num_of_rows_in_submatrix = n/wr;
	Matrix temp = Matrix(m,n);
	int arr[n] = {0};
	for(int i=0;i<n;i++) {
		int t = i/wr;
		arr[i] = t;
	}
	for(int i=0;i<wc;i++){
		for(int j=n-1;j>=1;j--){
			unsigned long long int num = distribution(gen);
			int ind = num%(j+1);                // random_number_in_range_0_to_j
			int temp_1 = arr[j];
			arr[j] = arr[ind];
			arr[ind] = temp_1;
		}
		for(int k=0;k<n;k++){
			temp.values[k].push_back(i*num_of_rows_in_submatrix + arr[k]);
		}
	}
	return temp;
}

Matrix MacKay_Neal(int n, int m, vector<float> v, vector<float> h,int wc=-1,int wr=-1) {
	Matrix mat = Matrix(m,n);
	int cd[n] = {};
	int rd[n] = {};
	int cd_ind = 0;
	int rd_ind = 0;
	for(int i=0;i<v.size();i++) {
		int temp2= v[i]*n;
		for(int j=0;j<temp2;j++) {
			cd[cd_ind] = i;
			cd_ind += 1;
		}
	}
	for(int i=0;i<h.size();i++) {
		int temp2= h[i]*m;
		for(int j=0;j<temp2;j++) {
			rd[rd_ind] = i;
			rd_ind += 1;
		}
	}
	if(wc>-1){
		for(int j=0;j<n;j++){
			cd[j] = wc;
		}
	}
	if(wr>-1){
		for(int j=0;j<m;j++){
			rd[j] = wr;
		}
	}
	int col_ind = 0;
	while (true) {
		int remaining_num_of_cols = n - col_ind;
		int col_weight = cd[col_ind];
		int to_select[col_weight] = {0};
		int to_select_remaining = col_weight;
		int last_ind = col_weight - 1;
		int temp[m] = {0};
		int temp_size = 0;
		for(int i=0;i<m;i++) {
			if (rd[i] > 0) {
				if (rd[i] == remaining_num_of_cols) {
					to_select[last_ind] = i;
					last_ind -= 1;
					to_select_remaining -= 1;
				}
				else {
					temp[temp_size] = i;
					temp_size += 1;
				}
			}
		}
		for(int i=0;i<to_select_remaining;i++) {
			to_select[i] = temp[i];
		}
		for(int i=to_select_remaining;i<temp_size;i++) {
			unsigned long long int num = distribution(gen);
			int ind = num%(i+1);               // random_number_in_range_0_to_i
			if (ind < to_select_remaining) {
				to_select[ind] = temp[i];
			}
		}
		for(int i=0;i<col_weight;i++) {
			rd[to_select[i]] -= 1;
			mat.values[col_ind].push_back(to_select[i]);
		}
		col_ind += 1;
		if (col_ind >= n) {
			break;
		}
	}
	return mat;
}


class TannerGraph{
private:
	int keyLength;
	int syndromeLength;
	vector<float> allBitNodes;
	vector<int> allParityNodes;
	vector<unordered_map<int,float>> parityEstimates;
	vector<float> bitEstimates; // Not needed for hard decoding
	float LLR; // Initial prob for soft decoding
	vector<int> targetKey;
public:
	Matrix checkMatrix;
	Matrix transposeCheckMatrix;
	TannerGraph(int keyLength,int syndromeLength,float density):keyLength(keyLength),syndromeLength(syndromeLength){
		this->checkMatrix = Matrix(syndromeLength,keyLength);
		this->checkMatrix.generateRandomMatrix(density);
		this->allBitNodes = vector<float>(keyLength,0);
		this->allParityNodes = vector<int>(syndromeLength,0);
		this->transposeCheckMatrix = this->checkMatrix.getTranspose();
		this->parityEstimates = vector<unordered_map<int,float>>(syndromeLength,unordered_map<int,float>());
		this->bitEstimates = vector<float>();
		this->targetKey = vector<int>(keyLength,0);
	}
	template <typename T>
	void initKey(vector<T> key){
		for(int i=0;i<keyLength;i++){
			allBitNodes[i] = static_cast<float>(key[i]);
		}
	}
	void initParity(vector<int> syndrome){
		for(int i=0;i<syndromeLength;i++){
			allParityNodes[i] = static_cast<int>(syndrome[i]);
		}
	}
	template <typename T>
	void initTarget(vector<T> key){
		for(int i=0;i<keyLength;i++){
			targetKey[i] = static_cast<int>(key[i]);
		}
	}
	void initSoft(float prob){ // Prob of an error
		float onePriorLLR = log(prob/(1-prob)); // If bit is 1 then likelihood of 0 is prob
		float zeroPriorLLR = log((1-prob)/prob);
		for(int i=0;i<keyLength;i++){
			if(allBitNodes[i] == 1){
				bitEstimates.push_back(onePriorLLR);
			}
			else bitEstimates.push_back(zeroPriorLLR);
		}
		this->LLR = zeroPriorLLR;
	}
	template <typename T>
	void initSystem(vector<T> key,vector<int> syndrome,bool flag,float prob=0,vector<T> target=vector<int>()){
		this->initKey(key);
		if(flag == 1){
			this->initSoft(prob);
		}
		this->initParity(syndrome);
		for(int i=0;i<syndromeLength;i++){
			for(int j=0;j<transposeCheckMatrix.values[i].size();j++){
				int temp = transposeCheckMatrix.values[i][j];
				parityEstimates[i][temp] = 0;
			}
		}
		if(targetKey.size()>0){
			initTarget(target);
		}
	}
	bool checkSolved(bool flag){
		if(flag == 0){
			cout << "Hey" << endl;
			for(int i=0;i<syndromeLength;i++){
				int currSum = 0;
				for(int j=0;j<transposeCheckMatrix.values[i].size();j++){
					int temp = transposeCheckMatrix.values[i][j];
					cout << allBitNodes[temp] << " " << endl;
					currSum = currSum ^ static_cast<int>(allBitNodes[temp]);
				}
				if(currSum != allParityNodes[i]){
					return false;
				}
			}
			return true;
		}
		for(int i=0;i<syndromeLength;i++){
			int currSum = 0;
			for(int j=0;j<transposeCheckMatrix.values[i].size();j++){
				float temp = transposeCheckMatrix.values[i][j];
				temp = bitEstimates[static_cast<int>(temp)];
				int hardDecision = 0;
				if(temp<=0){
					hardDecision = 1;
				}
				currSum = currSum ^ hardDecision;
			}
			if(currSum != allParityNodes[i]){
				return false;
			}
		}
		return true;
	}
	void transferParityBit(bool flag){
		if(flag == 0){
			for(int i=0;i<syndromeLength;i++){
				int currSum = 0;
				for(int j=0;j<transposeCheckMatrix.values[i].size();j++){
					int temp = transposeCheckMatrix.values[i][j];
					currSum = currSum ^ static_cast<int>(allBitNodes[temp]);
				}
				currSum = currSum ^ static_cast<int>(allParityNodes[i]);
				for(int j=0;j<transposeCheckMatrix.values[i].size();j++){
					int temp = transposeCheckMatrix.values[i][j];
					parityEstimates[i][temp] = currSum ^ static_cast<int>(allBitNodes[temp]);
				}
			}
			return;
		}
		for(int i=0;i<syndromeLength;i++){
			vector<float> tanhValues; // Values to be sent to the corresponding bits
			vector<int> signs; // Signs of the tanh values
			int zeroCount = 0; // Number of zeros in the values
			float logProd = 0; // LogProduct of the array
			int signProd = 1;
			for(int j=0;j<transposeCheckMatrix.values[i].size();j++){
				int temp = transposeCheckMatrix.values[i][j];
				float value = tanh((bitEstimates[temp]-parityEstimates[i][temp])/2); // LLR of the bit got by other check nodes
				tanhValues.push_back(abs(value)); // Push back the absolute value
				if(value == 0){ // If it is zero
					zeroCount++;
				}
				else logProd = logProd + log(abs(value)); // LogProduct of the array
				if(value>=0){
					signs.push_back(1); // Store the sign of the tanh
				}
				else signs.push_back(-1);
				signProd = signProd*signs[signs.size()-1];
			}
			for(int j=0;j<transposeCheckMatrix.values[i].size();j++){
				if(transposeCheckMatrix.values[i].size() == 1){ // This is the only connected bit
					tanhValues[j] = 0; // No other bit to get information from
				}
				else if(zeroCount>1){ // There are two zeros
					tanhValues[j] = 0; 
				}
				else if(zeroCount == 1 and tanhValues[j]!=0){ // There is one other zero
					tanhValues[j] = 0;
				}
				else if(tanhValues[j] == 0){ // This is the only zero
					tanhValues[j] = signs[j]*exp(logProd);
				}
				else{
					tanhValues[j] = (signProd*signs[j])*exp(logProd-log(tanhValues[j])); 
				}
				int temp = transposeCheckMatrix.values[i][j];
				parityEstimates[i][temp] = 2*atanh(tanhValues[j]);
				if(allParityNodes[i] == 1){
					parityEstimates[i][temp] = -1*parityEstimates[i][temp];
				}
				if(isinf(parityEstimates[i][temp]) or isnan(parityEstimates[i][temp])){
					cout << parityEstimates[i][temp] <<  " " << tanhValues[j] << " " << logProd << " " << bitEstimates[temp] << " " << zeroCount << endl;;
				}
			}
		}
		return;
	}
	void updateBit(bool flag){
		if(flag == 0){
			for(int i=0;i<keyLength;i++){
				if(checkMatrix.values[i].size() == 0){
					continue;
				}
				int currSum = 0;
				for(int j=0;j<checkMatrix.values[i].size();j++){
					int temp = checkMatrix.values[i][j];
					currSum = currSum + parityEstimates[temp][i];
				}
 				float final = float(currSum)/ checkMatrix.values[i].size();
				if(final>=0.5){
					allBitNodes[i] = 1;
				}
				else allBitNodes[i] = 0;
			}
			return;
		}
		for(int i=0;i<keyLength;i++){
			float initial = bitEstimates[i];
			bitEstimates[i] = 0;
			for(int j=0;j<checkMatrix.values[i].size();j++){
				int temp = checkMatrix.values[i][j];
				bitEstimates[i] = bitEstimates[i] + parityEstimates[temp][i];
			}
			if(allBitNodes[i] == 1){
				bitEstimates[i] = bitEstimates[i] - LLR;
			}
			else bitEstimates[i] = bitEstimates[i] + LLR;
			if(initial*bitEstimates[i]<0){
				cout << "Bit Flipped " << i << endl;
			}
		}
	}
	void decode(bool flag,int iterations,bool targetEnabled=false){
		if(flag == 0){
			decodeHard(iterations,targetEnabled);
			return;
		}
		decodeSoft(iterations,targetEnabled);
	}
	void decodeHard(int iterations,bool targetEnabled){
		for(int iter=0;iter<iterations;iter++){
			if(targetEnabled){
				cout << "Current Errors to Target " << getErrorsToTarget(0) << endl;
			} 
			if(checkSolved(0) == true){
				cout << "Solved" << endl;
				return;
			}
			transferParityBit(0);
			updateBit(0);
		}
	}
	void decodeSoft(int iterations,bool targetEnabled){
		for(int iter=0;iter<iterations;iter++){
			if(targetEnabled){
				cout << "Current Errors to Target " << getErrorsToTarget(1) << endl; 
			}
			if(checkSolved(1) == true){
				cout << "Solved" << endl;
				break;
			}
			transferParityBit(1);
			updateBit(1);
		}
		for(int i=0;i<keyLength;i++){
			if(bitEstimates[i]>0){
				allBitNodes[i] = 0;
			}
			else allBitNodes[i] = 1;
		}
	}
	vector<int> extractKey(bool flag){
		vector<int> temp;
		for(int i=0;i<keyLength;i++){
			temp.push_back(static_cast<int>(allBitNodes[i]));
		}
		return temp;
	}
	int getErrorsToTarget(bool flag){
		if(flag == 0){
			int errors = 0;
			for(int i=0;i<keyLength;i++){
				if(allBitNodes[i]!=this->targetKey[i]){
					errors++;
				}
			}
			return errors;
		}
		int errors = 0;
		for(int i=0;i<keyLength;i++){
			int hardDecision = 0;
			if(bitEstimates[i]<=0){
				hardDecision = 1;
			}
			if(hardDecision!=this->targetKey[i]){
				errors++;
			}
		}
		return errors;
	}
};


int getErrors(vector<int> aliceKey,vector<int> bobKey){
	int sum = 0;
	for(int i=0;i<aliceKey.size();i++){
		sum = sum + int(aliceKey[i]!=bobKey[i]);
	}
	return sum;
}


int main(int argc,char** argv){
	int keyLength = atoi(argv[1]);
	float errorFraction = atof(argv[2]);
	float syndromeFraction = atof(argv[3]);
	float matrixDensity = atof(argv[4]);
	int flag = atoi(argv[5]);
	int iterations = atoi(argv[6]);
	int syndromeLength = int(syndromeFraction*keyLength);
	vector<int> aliceKey;
	vector<int> bobKey;
	for(int i=0;i<keyLength;i++){
		float rand = dist(gen);
		if(rand<0.5)aliceKey.push_back(0);
		else aliceKey.push_back(1);
	}
	for(int i=0;i<keyLength;i++){
		float rand = dist(gen);
		if(rand<errorFraction){
			bobKey.push_back(aliceKey[i]^1);
		}
		else {
			bobKey.push_back(aliceKey[i]);
		}
	}
	int totalErrors = getErrors(aliceKey,bobKey);
	cout << "Keys" << endl;
	for(int i=0;i<keyLength;i++){
		cout << aliceKey[i] << " " << bobKey[i] << endl;
 	}
	vector<int> syndrome = vector<int>(syndromeLength,0);
	TannerGraph tannerGraph = TannerGraph(keyLength,syndromeLength,matrixDensity);
	multiply(aliceKey,syndrome,tannerGraph.checkMatrix);
	cout << "Syndrome " << endl;
	for(int i=0;i<syndromeLength;i++){
		cout << syndrome[i] << endl;
	}
	tannerGraph.initSystem(bobKey,syndrome,flag,errorFraction,aliceKey);
	cout << "Initial Errors " << totalErrors << endl;
	tannerGraph.decode(flag,iterations,false);
	bobKey = tannerGraph.extractKey(flag);
	for(int i=0;i<keyLength;i++){
		cout << aliceKey[i] << " " << bobKey[i] << endl;
 	}
	cout << "Final Errors " << getErrors(aliceKey,bobKey) << endl;
	return 0;
}
