#ifndef BOUND_FORMULA_H_
#define BOUND_FORMULA_H_

#include "sketch.h"

struct BoundFormula {
	int index;
	vector<Sketch*> uncList;
	vector<Sketch*> conList;
	vector<int> activeL; //global attribute ids
	vector<int> hash_sizes; //global attribute id to hash sizes

	BoundFormula(int _index, vector<Sketch*>& _uncList, vector<Sketch*>& _conList, vector<int>& _activeL, vector<int>& _hash_sizes) :
		index(_index), uncList(_uncList), conList(_conList), activeL(_activeL), hash_sizes(_hash_sizes) {}

	long execute(vector<int> indices) {
		long prod = 1;
		for (Sketch* s : uncList)
			prod *= s->access(index, -1, indices);
		for (int i = 0; i < conList.size(); i++)
			prod *= conList[i]->access(index, activeL[i], indices);
		return prod;
	}

	void print() {
		cout << "Bound Index: " << index << endl;
		cout << "Unconditional Sketchs: ";
		for (Sketch* s : uncList) {
			if (ZeroDimensionalSketchUnc* v = dynamic_cast<ZeroDimensionalSketchUnc*>(s)) {
				cout << s->t << " {0D}, ";
			}
			else if (OneDimensionalSketchUnc* v = dynamic_cast<OneDimensionalSketchUnc*>(s)) {
				cout << s->t << " {1D}, ";
			}
			else if (TwoDimensionalSketchUnc* v = dynamic_cast<TwoDimensionalSketchUnc*>(s)) {
				cout << s->t << " {2D}, ";
			}
		}
		cout << endl;
		cout << "Conditional Sketchs: ";
		for (int i = 0; i < conList.size(); i++) {
			if (ZeroDimensionalSketchCon* v = dynamic_cast<ZeroDimensionalSketchCon*>(conList[i])) {
				cout << conList[i]->t << " {0D}[" << activeL[i] << "], ";
			}
			else if (OneDimensionalSketchCon* v = dynamic_cast<OneDimensionalSketchCon*>(conList[i])) {
				cout << conList[i]->t << " {1D}[" << activeL[i] << "], ";
			}
			else if (TwoDimensionalSketchCon* v = dynamic_cast<TwoDimensionalSketchCon*>(conList[i])) {
				cout << conList[i]->t << " {2D}[" << activeL[i] << "], ";
			}
		}
		cout << endl;
		cout << "Hash Sizes" << endl;
		for (int h : hash_sizes)
			cout << h << ", "; 
		cout << endl;
	}
};

struct CrossProductIterator {
	vector<int> hashSizes;
	vector<int> indices;
	int totalBuckets;
	int count;
	vector<int> notOne;

	CrossProductIterator(vector<int>& _hashSizes) : 
		hashSizes(_hashSizes), indices(_hashSizes.size()) 
	{
		totalBuckets = 1;
		for (int h : hashSizes)
			totalBuckets *= h;
		count = 0;

		int numNotOnes = 0;
		for (int h : hashSizes)
			if (h != 1)
				numNotOnes++;
		notOne.resize(numNotOnes);
		int currPos = 0;
		for (int i = hashSizes.size() - 1; i > -1; i--) {
			if (hashSizes[i] > 1) {
				notOne[currPos] = i;
				currPos++;
			}
		}
		if (notOne.size() > 0)
			indices[notOne[0]] = -1;
	}

	bool hasNext() {
		return count < totalBuckets;
	}

	vector<int> next() {
		for (int i : notOne) {
			if (indices[i] < hashSizes[i]-1) {
				indices[i]++;
				count++;
				break;
			}
			else
				indices[i] = 0;
		}
		return indices;
	}
};

#endif
