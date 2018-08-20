/*
 * comms.h
 *
 *  Created on: 2015-09-23
 *      Author: kassy
 */

#ifndef COMMS_H_
#define COMMS_H_

using namespace arma;
using namespace std;


long genRandomSeed();

void readParityFile(umat &H, umat &G, int &n, int &k);

void checkForValidFilename(ifstream &inputFile);

void getUserDecodingScheme(bool &SPA, bool &GAA);

void getSNRValues(vector<double> &dBSNR);

void getPuncturingFile(vector<int> &bitsToPuncture, int &totalNumBitsToPuncture);

const char* createOutputFileName(int n, int k);

enum State {

	GenerateInfo = 0,
	Encode = 1,
	Puncture = 2,
	Modulate = 3,
	Channel = 4,
	Unpuncture = 5,
	Decode = 6,
	Resend = 7
};

#endif /* COMMS_H_ */
