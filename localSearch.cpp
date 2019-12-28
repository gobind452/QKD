#include <iostream>
#include <vector>
#include <ctime>
#include <queue>
#include <algorithm>
#include <random>
#include <cmath>
#include <utility>
#include <cstdlib>
#include <unordered_set>
#include <iterator>
using namespace std;

default_random_engine gen(time(0));
uniform_real_distribution<float> dist(0.0,1.0);

int keyLength;
float errorFraction; // Error Probability for the Bernoulli Process
float syndromeFraction;
float matrixDensity;
int syndromeLength;

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
};

Matrix checkMatrix;

vector<bool> aliceKey;
vector<bool> bobKey;

vector<int> aliceSyndrome;
vector<int> bobSyndrome;
int currDistance = 0;
int currDepth = 0;

unordered_set<int> flipsDone;
int statesVisited = 0;

void multiply(vector<bool>& key,vector<int>& toCalculate){
	for(int i=0;i<key.size();i++){
		if(key[i] == true){
			for(int j=0;j<checkMatrix.values[i].size();j++){
				toCalculate[checkMatrix.values[i][j]]++;
			}
		}
	}
	return;
}

bool compare(const pair<int, int>&i, const pair<int, int>&j)
{
    return i.first < j.first;
}

bool reverseCompare(const pair<int,int>&i,const pair<int,int>&j){
	return i.first > j.first;
}

int getDistance(){
	float norm = 0;
	for(int i=0;i<syndromeLength;i++){
		norm = norm + pow(aliceSyndrome[i]-bobSyndrome[i],2);
	}
	return norm;
}

int getNewDistance(int flipPosition,bool flip){
	int distance = currDistance;
	int flipType = -int(flip)+int(1-flip);
	for(int i=0;i<checkMatrix.values[flipPosition].size();i++){
		int row = checkMatrix.values[flipPosition][i];
		distance = distance - pow(bobSyndrome[row]-aliceSyndrome[row],2);
		distance = distance + pow(bobSyndrome[row]+flipType-aliceSyndrome[row],2);
	}
	return distance;
}

void changeSyndrome(int flipPosition,bool flip){
	int flipType = -int(flip)+int(1-flip);
	for(int i=0;i<checkMatrix.values[flipPosition].size();i++){
		int row = checkMatrix.values[flipPosition][i];
		bobSyndrome[row] = bobSyndrome[row] + flipType;
	}
}

int getErrors(){
	int sum = 0;
	for(int i=0;i<keyLength;i++){
		sum = sum + int(aliceKey[i]!=bobKey[i]);
	}
	return sum;
}

void makeMove(int flip,int distance,int reverse);

bool performDepthFirstSearch(int depthLeft){
	if(currDistance == 0){
		return true;
	}
	statesVisited++;
	vector<pair<int,int>> distances;
	int prevDistance = currDistance;
	for(int i=0;i<keyLength;i++){
		if(flipsDone.find(i) == flipsDone.end()){
			distances.push_back(pair<int,int>(getNewDistance(i,bobKey[i]),i));
		}
	}
	sort(distances.begin(),distances.end(),compare);
	for(int i=0;i<distances.size();i++){
		pair<int,int> elem = distances[i];
		makeMove(elem.second,elem.first,0);
		if(currDistance == 0){
			return true;
		}
		bool done = false;
		if(depthLeft > 0){
			done = performDepthFirstSearch(depthLeft-1);
		}
		if(done == true){
			return done;
		}
		makeMove(elem.second,prevDistance,1);
	}
	return false;
}

void makeMove(int flip,int distance,int reverse=0){
	changeSyndrome(flip,bobKey[flip]);
	bobKey[flip] = !bobKey[flip];
	if(reverse == 1){
		flipsDone.erase(flip);
	}
	else flipsDone.insert(flip);
	currDistance = distance;
}

bool performLocalSearch(float randomProb){
	while(currDistance != 0){
		statesVisited++;
		pair<int,int> greedyChoice = pair<int,int>(currDistance,-1);
		for(int i=0;i<keyLength;i++){
			if(flipsDone.find(i) == flipsDone.end()){
				pair<int,int> curr = pair<int,int>(getNewDistance(i,bobKey[i]),i);
				if(curr.first<greedyChoice.first){
					greedyChoice.first = curr.first;
					greedyChoice.second = curr.second;
				}
			}
		}
		if(greedyChoice.second == -1){
			unordered_set<int> :: iterator itr;
			vector<int> toReverse;
			for(itr = flipsDone.begin();itr!=flipsDone.end();itr++){
				if(dist(gen)<randomProb){
					toReverse.push_back(*itr);
				}
			}
			for(int i=0;i<toReverse.size();i++){
				pair<int,int> curr = pair<int,int>(getNewDistance(toReverse[i],bobKey[toReverse[i]]),toReverse[i]);
				makeMove(curr.second,curr.first,1);
			}
		}
		else{
			makeMove(greedyChoice.second,greedyChoice.first);
		}
	}
}

int main(int argc,char** argv){
	keyLength = atoi(argv[1]);
	errorFraction = atof(argv[2]);
	syndromeFraction = atof(argv[3]);
	matrixDensity = atof(argv[4]);
	syndromeLength = int(syndromeFraction*keyLength);
	checkMatrix = Matrix(syndromeLength,keyLength);
	checkMatrix.generateRandomMatrix(matrixDensity);
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
	int totalErrors = getErrors();
	cout << "Initial errors " << totalErrors << endl;
	for(int i=0;i<syndromeLength;i++){
		aliceSyndrome.push_back(0);
		bobSyndrome.push_back(0);
	}
	multiply(aliceKey,aliceSyndrome);
	multiply(bobKey,bobSyndrome);
	currDistance = getDistance();
	cout << "Initial Distance " << currDistance << endl;
	//performDepthFirstSearch(totalErrors);
	performLocalSearch(0.1);
	cout << "States Visited " << statesVisited << endl;
	cout << "Final Distance " << getDistance() << endl;
	cout << "Final Errors " << getErrors() << endl;
	return 0;
}
