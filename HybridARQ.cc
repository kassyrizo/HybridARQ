#include <iostream>
#include <armadillo>
#include <time.h>
#include <numeric>
#include <iomanip>

#include "comms.h"
#include "decoders.h"

using namespace std;
using namespace arma;

const long maxWordsPerSNR = 1000000;
const long minErrsPerSNR = 100;
const int  maxNumIterations = 10; // The maximum number of iterations that the sum-product algorithm is permitted to run

int main(int argc, char** argv)
{
	srand(genRandomSeed()); // Set random seed so all random variables are different every time

	// read in the LDPC matrix from binary file
	int  n, k;
	umat H, G;
	readParityFile(H, G, n, k);

	// Display the properties of the parity matrix for the user
	cout << "\nThe code has the following properties >> " << endl;
	cout << "\nn = " << n << " k = " << k <<". The number of check nodes is " << H.n_rows << "." << endl;

	// Get decoding algorithm
	bool SPA, GAA;//BF;
	getUserDecodingScheme(SPA, GAA);

	// Puncturing variables
	int totalBitsToPuncture = n - k; // total number of bits that can be punctured from the code
	int rate0_7Bits = n - ( (double) k / 0.7 ); // number of bits punctured corresponding to rate 0.7

	// Set the SNR values to run this simulation for
	vector<double> dBSNR;
	getSNRValues(dBSNR);

	// Create CSV file for output
	ofstream dataFile;
	const char* fileName = createOutputFileName(n,k);
	dataFile.open(fileName);

	// Display the start of the simulation`
	time_t ctt = time(0);
	cout << "\n-- Hybrid ARQ Simulation for (" << n << ", " << k << ") LDPC Code --" << endl;
	cout << "Simulation start on " << asctime(localtime(&ctt)) << endl;

	// Simulation variables
	unsigned long int numWords = 0;
	unsigned long int erasures = 0;
	int numBitsToPuncture = totalBitsToPuncture;
	int numBitsToSend = k;
	int numIterations = 0;
	int bitStepSize = 1;
	int totalBits = 0;
	int tempNumBits = 0;
	int bits = 0;
	bool wordError = false;
	bool insertBits = false;
	double codeRate = 1.0;
	double numMults = 0.0;
	double numLogics = 0.0;
	double numAdds = 0.0;
	double numCacheReads = 0.0;
	double numCacheWrites = 0.0;
	double numRamReads = 0.0;
	double numRamWrites = 0.0;
	State state = GenerateInfo;
	vector<double> ratesForSNR;
	vector<double> multsPerSNR;
	vector<double> logicsPerSNR;
	vector<double> addsPerSNR;
	vector<double> cacheReadsPerSNR;
	vector<double> cacheWritesPerSNR;
	vector<double> ramReadsPerSNR;
	vector<double> ramWritesPerSNR;
	vector<double> iterationsPerSNR;

	// Communication variables
	mat  	informationBits; // The information bits to be encoded (1 x k bit matrix)
	int* 	codeWord = new int[n](); // The encoded bits (1 x n bit matrix)
	int*	bitsToSend = new int[n]();
	int*	nextBitsToSend = new int[bitStepSize]();
	int* 	modulatedBitsToSend = new int[n]();
	mat  	noiseVector; // The error vector (i.e. noise) (1 x n bit matrix)
	double*	softReceivedVector = new double[n](); // The received vector distorted by channel noise (1 x n bit matrix)
	double* softReceivedBits = new double[bitStepSize]();
	int*   	estimatedCodeWord = new int[n]();
	int* 	estimatedInformationBits = new int[k]();

	//Helper Objects
	BlockEncoder blockEncoder 	(&G);
	CodePuncturer codePuncturer (n, k);
	SumProductDecoder sumProductDecoder (&H, maxNumIterations);
	//BitFlipDecoder bitFlipDecoder       (&H, maxNumIterations);
	GallagerADecoder gallagerADecoder(&H, maxNumIterations);

	// Calculate the noise power for each SNR
	mat linearSNR = exp10( conv_to<mat>::from(dBSNR)/10 );
	mat noiseVariance = 1 / (linearSNR) / 2;
	mat noiseStd = sqrt(noiseVariance);

	// Begin Simulation, run for each SNR value
	for(int snr = 0; snr < (int)dBSNR.size(); snr++)
	{
		printf("\n-- SNR: %g dB -- \n", dBSNR.at(snr));
	    printf(" Number of erasures:        ");

	    // Reset the simulation variables to start counting again for new SNR
	    numWords = 0;
	    erasures = 0;

	    numBitsToPuncture = totalBitsToPuncture; // Initially puncture all bits
	    numBitsToSend = k;	//Initially send only k bits

	    // Reset the state to start sending a new word, as well as all state variables
	    state = GenerateInfo;
	    insertBits = false;
	    numIterations = 0;
	    bitStepSize = 1;
	    totalBits = 0;
	    bits = 0;

		ratesForSNR.clear();
		multsPerSNR.clear();
		logicsPerSNR.clear();
		addsPerSNR.clear();
		cacheReadsPerSNR.clear();
		cacheWritesPerSNR.clear();
		ramReadsPerSNR.clear();
		ramWritesPerSNR.clear();
		iterationsPerSNR.clear();


	    while(numWords < maxWordsPerSNR && erasures < minErrsPerSNR)
	    {
	    	switch(state)
	    	{
	    		case GenerateInfo	: // Data Source, generate digital info
	    		{
	                printf("\b\b\b\b\b\b\b%7ld", erasures);
	                fflush(0);

	                wordError = false;
	    			insertBits = false;
	    			numIterations = 0;
					bitStepSize = 1;
					totalBits = 0;
					bits = 0;

					sumProductDecoder.ResetComplexity();
					gallagerADecoder.ResetComplexity();

	    			informationBits.randu(1, k);
	    			informationBits = round(informationBits);

	    			/*cout << "\n" << std::setw(15) << "GenerateInfo: ";
	    			for(int i = 0; i < (int)informationBits.size(); i++)
	    				cout << std::setw(3) << informationBits(i) << " ";
	    			cout << endl;*/

	    			state = Encode; // the next step is to encode the information bits
	    			break;
	    		}

	    		case Encode	:	// Encoder, encode any generated information bits
	    		{
	    			blockEncoder.Encode( conv_to<umat>::from(informationBits), codeWord);

	    			/*cout << std::setw(15) << "Encode: ";
	    			for(int i = 0; i < n; i++)
	    				cout << std::setw(3) <<  codeWord[i] << " ";
	    			cout << endl;*/

	    			state = Puncture;
	    			break;
	    		}

	    		case Puncture	:
	    		{
	    			// These vectors need resizing, delete the old versions and re-populate with a new random puncturing scheme
	    			delete [] bitsToSend;
					numBitsToSend = n - numBitsToPuncture;
					bitsToSend = new int[n]();

					for(int i = 0; i < n; i++)
						bitsToSend[i] = codeWord[i]; // need to preserve the codeWord, so create a copy that can be punctured

					codePuncturer.ChangePuncturedBits(numBitsToPuncture, n, k);
					codePuncturer.Puncture(bitsToSend, numBitsToPuncture);

					/*cout << std::setw(15) << "Puncture: ";
					for(int i = 0; i < numBitsToSend; i++)
						cout << std::setw(3) << bitsToSend[i] << " ";
					cout << endl;*/

					state = Modulate;
	    			break;
	    		}

	    		case Modulate	:
	    		{
	    			delete[] modulatedBitsToSend;

	    			modulatedBitsToSend = new int[numBitsToSend]();

	    			for(int i = 0; i < numBitsToSend; i++)
	    				modulatedBitsToSend[i] = -2*bitsToSend[i] + 1;

	    			/*cout << std::setw(15) << "Modulate: ";
					for(int i = 0; i < numBitsToSend; i++)
						cout << std::setw(3) << modulatedBitsToSend[i] << " ";
					cout << endl;*/

	    			state = Channel;
	    			break;
	    		}

	    		case Channel	:
	    		{
	    			noiseVector.randn(1, numBitsToSend);
	    			noiseVector = noiseVector*noiseStd(snr);

	    			if(!insertBits)
	    			{
	    				delete[] softReceivedVector;
	    				softReceivedVector = new double[numBitsToSend]();

	    				for(int i = 0; i < numBitsToSend; i++)
	    					softReceivedVector[i] = (double)modulatedBitsToSend[i] + noiseVector[i];

	    				/*cout << std::setw(15) << "Channel: ";
						for(int i = 0; i < numBitsToSend ; i++)
							cout << std::setw(3) << softReceivedVector[i] << " ";
						cout << endl;*/
	    			}
	    			else if(insertBits)
	    			{
	    				delete[] softReceivedBits;
						softReceivedBits = new double[numBitsToSend]();

						for(int i = 0; i < numBitsToSend; i++)
							softReceivedBits[i] = (double)modulatedBitsToSend[i] + noiseVector[i];

						/*cout << std::setw(15) << "Channel: ";
						for(int i = 0; i < numBitsToSend ; i++)
							cout << std::setw(3) << softReceivedBits[i] << " ";
						cout << endl;*/
	    			}

	    			totalBits += numBitsToSend;
					codeRate = (double)k / (double)totalBits;

					//cout << "totalBits: " << totalBits << " codeRate: " << codeRate << endl;

					state = Unpuncture;
	    			break;
	    		}

	    		case Unpuncture	:
	    		{
	    			if(numBitsToPuncture > 0 && !insertBits)
	    			{
	    				codePuncturer.UnPuncture(softReceivedVector, numBitsToPuncture);

	    				/*cout << std::setw(15) << "Unpuncture: ";
						for(int i = 0; i < n; i++)
							cout << std::setw(3) << softReceivedVector[i] << " ";
						cout << endl;*/
	    			}
	    			else if(insertBits)
	    			{
	    				codePuncturer.InsertBits(softReceivedVector, softReceivedBits, nextBitsToSend, numBitsToSend);

	    				/*cout << std::setw(15) << "InsertBits: ";
						for(int i = 0; i < n; i++)
							cout << std::setw(3) << softReceivedVector[i] << " ";
						cout << endl;*/
	    			}

	    			state = Decode;
	    			break;
	    		}

	    		case Decode	:
	    		{
	    			if(codeRate >= 1.0)
					{
						for(int i = 0; i < n; i++)
							estimatedCodeWord[i] = softReceivedVector[i] > 0 ? 1 : -1;
					}
	    			else if(SPA) // If using the sum-product algorithm
					{
						sumProductDecoder.Decode(softReceivedVector, noiseVariance(snr), estimatedCodeWord);
						numIterations += sumProductDecoder.GetNumIterations();
					}
					else if(GAA /*&& (numBitsToPuncture != totalBitsToPuncture || insertBits)*/ ) // If using the GAA algorithm, and haven't punctured all bits
					{
						gallagerADecoder.DecodeWithErasures(softReceivedVector, estimatedCodeWord, k);
						numIterations += gallagerADecoder.GetNumIterations();
					}

	    			// Demodulate
					for(int i = n-k; i < n; i++)
						estimatedInformationBits[i-n+k] = -1*estimatedCodeWord[i] > 0;

					wordError = false;
					for(int i = 0; i < k; i++)
					{
						if(informationBits[i] != estimatedInformationBits[i])
						{
							wordError = true;
							break;
						}
					}

					/*cout << std::setw(15) << "Decode: ";
					for(int i = 0; i < k ; i++)
						cout << std::setw(3) << estimatedInformationBits[i] << " ";
					cout << endl;*/


					if(wordError)
						state = Resend;
					else
					{
						// If the codeword worked the first time, add one to the number of bits to puncture
						// If we are already at the highest possible code rate, we can't puncture any more
						tempNumBits = 0;
						tempNumBits = n - (double)k/codeRate;

						/*if( (tempNumBits == numBitsToPuncture) && (numBitsToPuncture != totalBitsToPuncture) )
						{
							numBitsToPuncture++;
							if( numBitsToPuncture > totalBitsToPuncture || (SPA && numBitsToPuncture < totalBitsToPuncture && numBitsToPuncture > rate0_7Bits) )
								numBitsToPuncture = totalBitsToPuncture;
						}
						else
							numBitsToPuncture = tempNumBits;*/
						numBitsToPuncture = totalBitsToPuncture;

						if(SPA)
							sumProductDecoder.GetComplexity(numMults, numLogics, numAdds, numCacheReads, numCacheWrites, numRamReads, numRamWrites);
						else if(GAA)
							gallagerADecoder.GetComplexity(numMults, numLogics, numAdds, numCacheReads, numCacheWrites, numRamReads, numRamWrites);

						multsPerSNR.push_back(numMults);
						logicsPerSNR.push_back(numLogics);
						addsPerSNR.push_back(numAdds);
						cacheReadsPerSNR.push_back(numCacheReads);
						cacheWritesPerSNR.push_back(numCacheWrites);
						ramReadsPerSNR.push_back(numRamReads);
						ramWritesPerSNR.push_back(numRamWrites);

						iterationsPerSNR.push_back(numIterations);

						ratesForSNR.push_back(codeRate);

						numWords++;
						state = GenerateInfo;
					}
	    			break;
	    		}

	    		case Resend	:
	    		{
	    			if( bits >= numBitsToPuncture )
	    			{
	    				erasures++;
	    				numWords++;
	    				state = GenerateInfo;

	    				//cout << "erasures: " << erasures << endl;
	    				break;
	    			}
	    			else if( numBitsToPuncture > rate0_7Bits && bits == 0 && SPA)
	    				numBitsToSend = totalBitsToPuncture - rate0_7Bits; //
	    			else
	    				numBitsToSend = 1;

	    			codePuncturer.BitsToPuncture(nextBitsToSend, bits, numBitsToSend);

	    			delete [] bitsToSend;
					bitsToSend = new int[numBitsToSend]();

					for(int i = 0; i < numBitsToSend; i++)
						bitsToSend[i] = codeWord[ nextBitsToSend[i] ];

					bits += numBitsToSend;

					/*cout << std::setw(15) << "Resend: ";
					for(int i = 0; i < numBitsToSend ; i++)
						cout << std::setw(3) << bitsToSend[i] << " ";
					cout << endl;

					cout << std::setw(15) << "numBitsToSend: " <<  numBitsToSend << endl;
					cout << std::setw(15) << "numBitsToPuncture: " << numBitsToPuncture << endl;*/
					insertBits = true;

					state = Modulate;
	    			break;
	    		}

	    		default:
	    			break;
	    	}
	    }


	    long double sumOfMults = accumulate( multsPerSNR.begin(), multsPerSNR.end(), 0.0 );
		double averageMults = sumOfMults / (double)multsPerSNR.size();

		long double sumOfLogic = accumulate( logicsPerSNR.begin(), logicsPerSNR.end(), 0.0 );
		double averageLogic = sumOfLogic / (double)logicsPerSNR.size();

		long double sumOfAdds  = accumulate( addsPerSNR.begin(), addsPerSNR.end(), 0.0 );
		double averageAdds  = sumOfAdds /  (double)addsPerSNR.size();

		long double sumOfCacheReads  = accumulate( cacheReadsPerSNR.begin(), cacheReadsPerSNR.end(), 0.0 );
		double averageCacheReads = sumOfCacheReads /  (double)cacheReadsPerSNR.size();

		long double sumOfCacheWrites  = accumulate( cacheWritesPerSNR.begin(), cacheWritesPerSNR.end(), 0.0 );
		double averageCacheWrites = sumOfCacheWrites /  (double)cacheWritesPerSNR.size();

		long double sumOfRamReads  = accumulate( ramReadsPerSNR.begin(), ramReadsPerSNR.end(), 0.0 );
		double averageRamReads = sumOfRamReads /  (double)ramReadsPerSNR.size();

		long double sumOfRamWrites  = accumulate( ramWritesPerSNR.begin(), ramWritesPerSNR.end(), 0.0 );
		double averageRamWrites = sumOfRamWrites /  (double)ramWritesPerSNR.size();

		double sumOfRates = accumulate(ratesForSNR.begin(), ratesForSNR.end(), 0.0);
		double averageRate = sumOfRates/(double)ratesForSNR.size();

		long double sumOfIterations = accumulate( iterationsPerSNR.begin(), iterationsPerSNR.end(), 0.0 );
		double averageIterations = sumOfIterations/ (double) iterationsPerSNR.size();

		double erasureRate = (double)erasures / (double) numWords;

		printf("\nAverage rate: %g, Erasure Rate: %g\n", averageRate, erasureRate);
		dataFile << ( 10*log10( (double) linearSNR(snr) / averageRate) ) << "," << averageRate << "," << erasureRate << "," << averageIterations << "," << averageMults << "," << averageLogic << "," << averageAdds << "," << averageCacheReads << averageCacheWrites << averageRamReads << averageRamWrites << "\n";
	}

	dataFile.close();
	ctt = time(0);

	cout << "\nSimulation completed on " << asctime(localtime(&ctt)) << endl;

	// Clean-up and delete any dynamic memory
	delete [] codeWord;
	delete [] modulatedBitsToSend;
	delete [] bitsToSend;
	delete [] nextBitsToSend;
	delete [] softReceivedVector;
	delete [] softReceivedBits;
	delete [] estimatedCodeWord;
	delete [] estimatedInformationBits;

	return 0;
}



