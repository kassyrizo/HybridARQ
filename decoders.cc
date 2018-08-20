#include <armadillo>
#include <time.h>
#include <sys/time.h>
#include <stdexcept>
#include <iomanip>
#include "decoders.h"


////////////////////////////////
// Puncturing Class Functions //
////////////////////////////////

CodePuncturer::CodePuncturer(int n, int k)
{
	sizeOfCodeWord = n;
	totalNumberOfPuncturedBits = n-k;
	bitsThatCanBePunctured = new int[totalNumberOfPuncturedBits]();

	return;
}

CodePuncturer::~CodePuncturer()
{
	delete [] bitsThatCanBePunctured;
}

void CodePuncturer::Puncture(int* &codeWord, int numberOfBitsToPuncture)
{
	if(numberOfBitsToPuncture > totalNumberOfPuncturedBits /*|| numberOfBitsToPuncture < 0*/)
		throw std::invalid_argument("cannot puncture %d that many bits! ");

	int* newPuncturedCodeWord = new int[sizeOfCodeWord - numberOfBitsToPuncture]();
	int* bitsToPuncture = new int[numberOfBitsToPuncture]();

	for(int i = 0; i < numberOfBitsToPuncture; i++)
		bitsToPuncture[i] = bitsThatCanBePunctured[i];

	sortBits(bitsToPuncture, numberOfBitsToPuncture);

	int* puncturePtr = bitsToPuncture;
	int  bitsPunctured = 0;
	for(int i = sizeOfCodeWord - 1; i >= 0; i--)
	{
		if((*puncturePtr) != i)
			newPuncturedCodeWord[i - numberOfBitsToPuncture + bitsPunctured] = codeWord[i];
		else
		{
			bitsPunctured++;
			if(bitsPunctured != numberOfBitsToPuncture)
				puncturePtr++;
		}
	}

	delete [] bitsToPuncture;
	delete [] codeWord;
	codeWord = newPuncturedCodeWord;

	return;
}

void CodePuncturer::UnPuncture(double* &puncturedWord, int numOfBitsPunctured)
{
	if(numOfBitsPunctured > totalNumberOfPuncturedBits)
			throw std::invalid_argument("cannot unpuncture that many bits!");

	double* newUnpuncturedCodeWord = new double[sizeOfCodeWord]();
	int* bitsToUnpuncture = new int[numOfBitsPunctured]();

	for(int i = 0; i < numOfBitsPunctured; i++)
		bitsToUnpuncture[i] = bitsThatCanBePunctured[i];

	sortBits(bitsToUnpuncture, numOfBitsPunctured);

	int* puncturePtr = (&bitsToUnpuncture[numOfBitsPunctured - 1]);
	int  bitsPutBack = 0;
	for(int i = 0; i < sizeOfCodeWord; i++)
	{
		if( (*puncturePtr) == i)
		{
			newUnpuncturedCodeWord[i] = 0;
			bitsPutBack++;
			puncturePtr--;
		}
		else
			newUnpuncturedCodeWord[i] = puncturedWord[i - bitsPutBack];
	}

	delete [] bitsToUnpuncture;
	delete [] puncturedWord;
	puncturedWord = newUnpuncturedCodeWord;

	return;
}

void CodePuncturer::BitsToPuncture(int* &bitsToSend, int index, int size)
{
	// Fill bitsToSend with the indices of the bits requested, sorted in order
	int* newBitsToSend = new int[size]();

	for(int i = 0; i < size; i++)
		newBitsToSend[i] = bitsThatCanBePunctured[i + index];

	sortBits(newBitsToSend, size);

	delete [] bitsToSend;
	bitsToSend = newBitsToSend;

	/*cout << std::setw(15) << "newBitsToSend: ";
	for(int i = 0; i < size; i++)
		cout << std::setw(3) << newBitsToSend[i] << " ";
	cout << endl;*/

	return;
}

void CodePuncturer::InsertBits(double* &receivedVector, double* softReceivedBits, int* bitsToInsert, int numBits)
{
	double* newReceivedVector = new double[sizeOfCodeWord + numBits]();

	int* puncturePtr = (&bitsToInsert[numBits - 1]);
	int  bitsPutBack = 0;

	/*cout << std::setw(15) << "bitsToInsert: ";
	for(int i = 0; i < numBits; i++)
		cout << std::setw(3) << bitsToInsert[i] << " ";
	cout << endl;*/

	for(int i = 0; i < sizeOfCodeWord; i++)
	{
		if( (*puncturePtr) == i )
		{
			newReceivedVector[i] = softReceivedBits[numBits - 1 - bitsPutBack];
			bitsPutBack++;
			puncturePtr--;
		}
		else
			newReceivedVector[i] = receivedVector[i];
	}

	delete [] receivedVector;
	receivedVector = newReceivedVector;

	return;
}

void CodePuncturer::ChangePuncturedBits(int numBitsToPuncture, int n, int k)
{
	totalNumberOfPuncturedBits = numBitsToPuncture;

	delete [] bitsThatCanBePunctured;

	bitsThatCanBePunctured = new int[numBitsToPuncture]();

	int number = 0;
	bool newBit = true;

	for(int i = 0; i < numBitsToPuncture; i++)
	{
		newBit = true;
		number = std::rand() % (n - k);

		for(int p = 0; p < i; p++)
		{
			if(bitsThatCanBePunctured[p] == number)
			{
				newBit = false;
				break;
			}
		}

		if(newBit)
			bitsThatCanBePunctured[i] = number;
		else
			--i;
	}

	/*cout << std::setw(15) << "PuncturedBits: ";
	for(int i = 0; i < numBitsToPuncture; i++)
		cout <<  std::setw(3) << bitsThatCanBePunctured[i] << " ";
	cout << endl;*/

	return;
}

void CodePuncturer::sortBits(int* bitsToSort, int numBitsToSort)
{
	// sort the bits in ascending order
	int temp;
	for(int i = 0; i < numBitsToSort; i++)
	{
		for(int j = 0; j < numBitsToSort - 1; j++)
		{
			if(bitsToSort[j] < bitsToSort[j+1])
			{
				temp = bitsToSort[j+1];
				bitsToSort[j+1] = bitsToSort[j];
				bitsToSort[j] = temp;
			}
		}
	}

	return;
}

////////////////////////////////////
/// Block Encoder Class Functions //
////////////////////////////////////
BlockEncoder::BlockEncoder(umat* Gptr)
{
	rows = (*Gptr).n_rows;
	cols = (*Gptr).n_cols;
	colSupport = new int*[cols]();
	colOnes = new int[cols]();

	// Allocate the memory for the columns in G
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find( (*Gptr).col(i) );
		colSupport[i] = new int[ones.size()]();
		colOnes[i] = ones.size();
		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}
}

BlockEncoder::~BlockEncoder()
{
	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];

	delete [] colSupport;
	delete [] colOnes;
}

void BlockEncoder::Encode(umat u, int* &codeWord)
{
	// The encoded word is basically the xor (mod 2 addition) of
	// the bits in u with row support in each column in G
	for(int i = 0; i < cols - rows; i++)
	{
		codeWord[i] = 0; // reset the bit value

		for(int j = 0; j < colOnes[i]; j++)
			codeWord[i] ^= u[ colSupport[i][j] ];
	}

	// Concat u on the end, since this is a systematic code
	for(int i = cols - rows; i < cols; i++)
		codeWord[i] = u[ i - (cols - rows) ];

	return;
}

///////////////////////////////////////////
// Sum Product Algorithm Class Functions //
///////////////////////////////////////////

SumProductDecoder::SumProductDecoder(umat* Hptr, int maxNumIterations)
{
	// Max number of iterations for the sum product algorithm
	totalIterations = maxNumIterations;
	numIterations = 0;

	// Parity matrix dimensions
	rows = (*Hptr).n_rows;
	cols = (*Hptr).n_cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	// Allocate the memory for the rows in H
	for(int i = 0; i < rows; i++)
	{
		uvec ones = find((*Hptr).row(i));
		rowSupport[i] = new int[ones.size()]();
		rowOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			rowSupport[i][j] = ones[j];
	}

	// Allocate the memory for the columns in  H
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find((*Hptr).col(i));
		colSupport[i] = new  int[ones.size()]();
		colOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}

	hardReceivedVector = new int[cols]();
	syndrome = new int[rows]();

	r = new double[cols]();
	eps_tot = new double[cols]();
	Y = new double*[rows]();
	Z = new double*[rows]();
	E = new double*[rows]();
	epsilon = new double*[rows]();


	//Memory allocation and initialization for SPA variables
	for(int i = 0; i < cols; i++)
	{
		r[i]   = 0;
		eps_tot[i] = 0;
	}

	for(int i = 0; i < rows; i++)
	{
		Y[i] = new double[rowOnes[i]]();
		Z[i] = new double[rowOnes[i]]();
		E[i] = new double[rowOnes[i]]();
		epsilon[i] = new double[rowOnes[i]];

		for(int j = 0; j < rowOnes[i]; j++)
		{
			Y[i][j] = 0;
			Z[i][j] = 0;
			E[i][j] = 0;
			epsilon[i][j] = 0;
		}
	}

	ResetComplexity();
}

SumProductDecoder::~SumProductDecoder()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
		delete [] rowSupport[i];
	delete [] rowSupport;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];
	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	// Delete memory allocations
	for(int i = 0; i < rows; i++)
	{
		delete [] Y[i];
		delete [] Z[i];
		delete [] E[i];
		delete [] epsilon[i];
	}

	delete [] hardReceivedVector;
	delete [] syndrome;

	delete [] Y;
	delete [] Z;
	delete [] E;
	delete [] epsilon;
	delete [] r;
	delete [] eps_tot;
}

void SumProductDecoder::Decode(double* receivedVector, double noiseVariance, int* &estimatedCodeWord)
{
	pSoftReceivedVector = receivedVector; // set the soft Received vector to point at the received vector

	/////////////////////////////////////
	// Begin the Sum Product Algorithm //
	/////////////////////////////////////

	//////////////////////////////////////////
	// Initialize the Sum Product Algorithm //
	//////////////////////////////////////////
	numIterations = 0;

	for(int i = 0; i < cols; i++)
	{
		r[i] = 2 * pSoftReceivedVector[i] / noiseVariance;
		eps_tot[i] = 0;
		estimatedCodeWord[i]  = r[i] > 0 ? 1 : -1;
		numLogic++;
	}

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < rowOnes[i]; j++)
		{
			Y[i][j] = r[rowSupport[i][j]];
			Z[i][j] = Y[i][j];
			E[i][j] = 0;
			epsilon[i][j] = 0;
		}
	}

	// Read H matrix, rows*6 numbers per row, * 16 bits for short int
	// Cache Paging here for reading H matrix
	totalHNumbers = (rows*6)*16;
	while(totalHNumbers > 0)
	{
		numCacheReads++;
		totalHNumbers -= 64; // We can read 64 bits at a time from Cache to Registers.
	}
	// RAM paging here for reading H matrix
	totalHNumbers = (rows*6)*16 - (16*1024); // subtract 16 kB, because that is how much can be stored in Cache
	while(totalHNumbers > 0)
	{
		numRamReads++;
		totalHNumbers -= 128;
	}

	// Read r, cols * 32 bits for SP float
	// Read Y, rows*6 * 32 bits for SP float
	// Cache paging for reading/setting
	totalVarNumbers = (rows*6)*32 + cols*32;
	while(totalVarNumbers > 0)
	{
		numCacheReads++;
		totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
	}
	// RAM paging
	totalVarNumbers = (rows*6)*32 + cols*32 -(16*1024) ; // Subtract 16 kB because that is how much can be stored in Cache
	while(totalVarNumbers > 0)
	{
		numRamReads++;
		totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
	}

	// Write r, cols * 32 bits for SP float
	// Write Y, rows*6 * 32 bits for SP float
	// Write Z, rows*6 * 32 bits for SP float
	totalVarNumbers = (rows*6)*32*2 + cols*32;
	while(totalVarNumbers > 0)
	{
		numCacheWrites++;
		totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
	}
	// RAM paging
	totalVarNumbers = (rows*6)*32*2 + cols*32 -(16*1024) ; // Subtract 16 kB because that is how much can be stored in Cache
	while(totalVarNumbers > 0)
	{
		numRamWrites++;
		totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
	}

	while(numIterations < totalIterations)
	{
		// Check the syndrome to see if it is all 0's, if it is, don't do the SPA
		for(int i = 0; i < cols; i++)
			hardReceivedVector[i] = -1*estimatedCodeWord[i] > 0;

		if(calculateSyndrome())
			break;

		// For each row, calculate the extrinsic information provided by the OTHER bits with support in the row, not including the bit itself
		double eps = 1.0;
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
			{
				eps = 1.0;
				for(int l = 0; l < rowOnes[i]; l++)
				{
					if(rowSupport[i][j] != rowSupport[i][l])
					{
						eps *= tanh( Z[i][l]/2 );
						numMults++;
						numLogic++;
					}
					numLogic++;
				}

				eps = log( (1+eps)/(1-eps) );
				epsilon[i][j] = eps;
				numAdds += 2;
				numLogic++;
			}
		}

		// Cache Paging here for reading H matrix
		totalHNumbers = (rows*6)*16; // Storing ints as short ints (16 bits), assume no bits stored in Regsiters, must read everything
		while(totalHNumbers > 0)
		{
			numCacheReads++;
			totalHNumbers -= 64; // We can read 64 bits at a time from Cache to Registers.
		}
		// RAM paging here for reading H matrix
		totalHNumbers = (rows*6)*16 - (16*1024); // Storing ints as short ints (16 bits), assume 16 kB in cache, the rest is stored in RAM
		while(totalHNumbers > 0)
		{
			numRamReads++;
			totalHNumbers -= 128; // We can read 16 kB at a time from Ram to Cache.(This probably needs to change, not sure exactly to what)
		}

		// Read Z, rows*6 *32 bits
		// Write epsilon, rows*6 *32 bits
		// Cache paging
		totalVarNumbers = (rows*6)*32;
		while(totalVarNumbers > 0)
		{
			numCacheWrites++;
			numCacheReads++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6)*32 - 16*1024; // Set epsilon
		while(totalVarNumbers > 0)
		{
			numRamWrites++;
			numRamReads++;
			totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
		}

		// Calculate the Extrinsic Matrix
		double sum = 0;
		for(int i = 0; i < cols; i++)
		{
			for(int j = 0; j < colOnes[i]; j++)
			{
				sum = 0;
				for(int l = 0; l < colOnes[i]; l++)
				{
					numLogic++;
					if(colSupport[i][j] != colSupport[i][l])
					{
						for(int t = 0; t < rowOnes[ colSupport[i][l] ]; t++)
						{
							numLogic++;
							if (rowSupport[ colSupport[i][l] ][t] == i)
							{
								sum += epsilon[ colSupport[i][l] ][t];
								numAdds++;
								break;
							}
						}
					}
				}

				for(int t = 0; t < rowOnes[ colSupport[i][j] ]; t++)
				{
					if(rowSupport[ colSupport[i][j]][t] == i)
					{
						E[ colSupport[i][j] ][t] = sum;
						break;
					}
				}
			}
		}

		// Cache Paging here for reading H matrix
		totalHNumbers = (rows*6)*16; // Storing ints as short ints (16 bits), assume no bits stored in Regsiters, must read everything
		while(totalHNumbers > 0)
		{
			numCacheReads++;
			totalHNumbers -= 64; // We can read 64 bits at a time from Cache to Registers.
		}
		// RAM paging here for reading H matrix
		totalHNumbers = (rows*6)*16 - (16*1024); // Storing ints as short ints (16 bits), assume 16 kB in cache, the rest is stored in RAM
		while(totalHNumbers > 0)
		{
			numRamReads++;
			totalHNumbers -= 128; // We can read 16 kB at a time from Ram to Cache.(This probably needs to change, not sure exactly to what)
		}

		// Read epsilon, rows*6 * 32 bits
		// Write E, rows*6 *32 bits
		// Cache paging
		totalVarNumbers = (rows*6)*32;
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6)*32 - 16*1024; // Read epsilon, set E
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			numRamWrites++;
			totalVarNumbers -= 128;
		}

		// Calculate epsilon_total
		for(int i = 0; i < cols; i++)
		{
			for(int t = 0; t < rowOnes[ colSupport[i][0] ]; t++)
			{
				numLogic++;
				if(rowSupport[ colSupport[i][0] ][t] == i)
				{
					eps_tot[i] = E[ colSupport[i][0] ][t] + epsilon[ colSupport[i][0] ][t];
					numAdds++;
					break;
				}
			}
		}

		// Cache Paging here for reading H matrix
		totalHNumbers = (rows*6)*16; // Storing ints as short ints (16 bits), assume no bits stored in Regsiters, must read everything
		while(totalHNumbers > 0)
		{
			numCacheReads++;
			totalHNumbers -= 64; // We can read 64 bits at a time from Cache to Registers.
		}
		// RAM paging here for reading H matrix
		totalHNumbers = (rows*6)*16 - (16*1024); // Storing ints as short ints (16 bits), assume 16 kB in cache, the rest is stored in RAM
		while(totalHNumbers > 0)
		{
			numRamReads++;
			totalHNumbers -= 128;
		}

		// Read E, rows*6 * 32 bits
		// Read epsilon, rows*6 * 32 bits
		// Cache paging
		totalVarNumbers = (rows*6)*32*2;
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6)*32*2 - 16*1024;
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			totalVarNumbers -= 128;
		}

		// Write eps_tot, cols * 32 bits
		// Cache paging
		totalVarNumbers = cols*32;
		while(totalVarNumbers > 0)
		{
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = cols*32 - 16*1024;
		while(totalVarNumbers > 0)
		{
			numRamWrites++;
			totalVarNumbers -= 128;
		}

		////////////////////////////////////////////////
		// Step 2: Form the next Z, r, and z matrices //
		////////////////////////////////////////////////

		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
			{
				Z[i][j] = Y[i][j]+E[i][j];
				numAdds++;
			}
		}

		for(int i = 0; i < cols; i++)
		{
			r[i] += eps_tot[i];
			estimatedCodeWord[i] = r[i] > 0 ? 1 : -1;

			numAdds++;
			numLogic++;
		}

		// Cache Paging here for reading H matrix
		totalHNumbers = (rows*6)*16; // Storing ints as short ints (16 bits), assume no bits stored in Regsiters, must read everything
		while(totalHNumbers > 0)
		{
			numCacheReads++;
			totalHNumbers -= 64; // We can read 64 bits at a time from Cache to Registers.
		}
		// RAM paging here for reading H matrix
		totalHNumbers = (rows*6)*16 - (16*1024); // Storing ints as short ints (16 bits), assume 16 kB in cache, the rest is stored in RAM
		while(totalHNumbers > 0)
		{
			numRamReads++;
			totalHNumbers -= 128;
		}

		// Read Y, rows*6 * 32 bits
		// Read E, rows*6 * 32 bits
		// Read eps_tot, cols * 32 bits
		// Cache paging
		totalVarNumbers = (rows*6)*32*2 + cols*32*2;
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6)*32*2 + cols*32*2 - 16*1024;
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			totalVarNumbers -= 128;
		}

		// Write Z, rows*6 * 32 bits
		// Write r, cols * 32 bits
		// Cache paging
		totalVarNumbers = (rows*6)*32 + cols*32;
		while(totalVarNumbers > 0)
		{
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6)*32 + cols*32 - 16*1024;
		while(totalVarNumbers > 0)
		{
			numRamWrites++;
			totalVarNumbers -= 128;
		}


		numIterations++;
	}

	return;
}

void SumProductDecoder::GetComplexity(double &mults, double &logic, double &adds, double &cacheReads, double &cacheWrites, double &ramReads, double &ramWrites)

{
	mults = numMults;
	logic = numLogic;
	adds  = numAdds;
	cacheReads  = numCacheReads;
	cacheWrites = numCacheWrites;
	ramReads = numRamReads;
	ramWrites = numRamWrites;
}

void SumProductDecoder::ResetComplexity()
{
	numMults = 0.0;
	numLogic = 0.0;
	numAdds  = 0.0;
	numCacheReads  = 0.0;
	numCacheWrites = 0.0;
	numRamReads = 0.0;
	numRamWrites = 0.0;
}

int SumProductDecoder::GetNumIterations()
{
	return numIterations;
}

bool SumProductDecoder::calculateSyndrome()
{
	bool syndromeCheck = true;

	// The syndrome is basically the xor (mod 2 addition) of
	// the bits codeword with support in each row in H
	for(int i = 0; i < rows; i++)
	{
		syndrome[i] = 0; // reset the bit value

		for(int j = 0; j < rowOnes[i]; j++)
		{
			syndrome[i] ^= hardReceivedVector[ rowSupport[i][j] ];
		}

		if(syndrome[i] == 1)
			syndromeCheck = false;
	}

	return syndromeCheck;
}

////////////////////////////////////////////
// Gallager's Algorithm A Class Functions //
////////////////////////////////////////////

GallagerADecoder::GallagerADecoder(umat* Hptr, int maxNumIterations)
{
	// Max number of iterations for the sum product algorithm
	totalIterations = maxNumIterations;
	numIterations = 0;

	// Parity matrix dimensions
	rows = (*Hptr).n_rows;
	cols = (*Hptr).n_cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	// Allocate the memory for the rows in H
	for(int i = 0; i < rows; i++)
	{
		uvec ones = find((*Hptr).row(i));
		rowSupport[i] = new int[ones.size()]();
		rowOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			rowSupport[i][j] = ones[j];
	}

	// Allocate the memory for the columns in  H
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find((*Hptr).col(i));
		colSupport[i] = new  int[ones.size()]();
		colOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}

	initHardDecision = new int[cols]();
	hardDecisionVector = new int[cols]();
	erasedBits = new int[cols]();
	syndrome = new int[rows]();

	varNodeVals  = new int*[rows]();
	parityChecks = new int*[rows]();

	for(int i = 0; i < rows; i++)
	{
		varNodeVals[i]  = new int[rowOnes[i]]();
		parityChecks[i] = new int[rowOnes[i]]();
	}

	peelingVector = new double[cols]();

	numMults = 0;
	numLogic = 0;
	numAdds = 0;
	numCacheReads = 0;
	numCacheWrites = 0;
	numRamReads = 0;
	numRamWrites = 0;

	/*for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < rowOnes[i]; j++)
			cout << rowSupport[i][j] << " ";
		cout << endl;
	}
	cout << endl;

	for(int i = 0; i < cols; i++)
	{
		for(int j = 0; j < colOnes[i]; j++)
			cout << colSupport[i][j] << " ";
		cout << endl;
	}
	cout << endl;*/

}

GallagerADecoder::~GallagerADecoder()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
		delete [] rowSupport[i];
	delete [] rowSupport;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];

	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	for(int i = 0; i < rows; i++)
	{
		delete [] parityChecks[i];
		delete [] varNodeVals[i];
	}

	delete [] initHardDecision;
	delete [] hardDecisionVector;
	delete [] erasedBits;
	delete [] syndrome;
	delete [] parityChecks;
	delete [] varNodeVals;
	delete [] peelingVector;
}

void GallagerADecoder::Decode(double* receivedVectorPtr, int* &estimatedCodeWord)
{
	for(int i = 0; i < cols; i++)
	{
		initHardDecision[i] = receivedVectorPtr[i] < 0;
		hardDecisionVector[i] = receivedVectorPtr[i] < 0;
		numLogic++;
	}

	// Write initHardDecision, cols bits (no complexity for hardDecision vector, since that will be identical for SPA)
	// Cache paging
	totalVarNumbers = cols;
	while(totalVarNumbers > 0)
	{
		numCacheWrites++;
		totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
	}
	// RAM paging
	totalVarNumbers = cols - 16*1024; // Read epsilon, read E, set eps_tot
	while(totalVarNumbers > 0)
	{
		numRamWrites++;
		totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
	}

	/*cout << "hrdDecision: ";
	for(int i = 0; i < cols; i++)
		cout << hardDecisionVector[i] << " ";
	cout << endl;*/

	// Initialize the algorithm using the hard decision guess
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < rowOnes[i]; j++)
			varNodeVals[i][j] = initHardDecision[rowSupport[i][j]];

	// Cache Paging here for reading H matrix
	totalHNumbers = (rows*6)*16; // rows*6 numbers in matrix * 16 bits (short int)
	while(totalHNumbers > 0)
	{
		numCacheReads++;
		totalHNumbers -= 64; // We can read/write 64 bits at a time from Cache to Registers.
	}
	// RAM paging here for reading H matrix
	totalHNumbers = (rows*6)*16 - (16*1024); // 16 kB can be stored in Cache, so remove from RAM calculation
	while(totalHNumbers > 0)
	{
		numRamReads++;
		totalHNumbers -= 128; // Can read/write 256 bits at a time from RAM
	}

	// Read initHardDecision, cols bits
	// Cache paging
	totalVarNumbers = cols;
	while(totalVarNumbers > 0)
	{
		numCacheReads++;
		totalVarNumbers -= 64; // Can read/write 64 bits at a time from Cache to Registers
	}
	// RAM paging
	totalVarNumbers = cols - 16*1024; // 16 kB can be stored in Cache, remove from RAM calculation
	while(totalVarNumbers > 0)
	{
		numRamReads++;
		totalVarNumbers -= 128; // Can read/write 256 bits at a time from RAM
	}

	// Write varNodeVals, rows*6 bits
	// Cache paging
	totalVarNumbers = rows*6;
	while(totalVarNumbers > 0)
	{
		numCacheWrites++;
		totalVarNumbers -= 64; // Can read/write 64 bits at a time from Cache to Registers
	}

	// RAM paging
	totalVarNumbers = (rows*6) - 16*1024; // 16 kB can be stored in Cache, remove from RAM calculation
	while(totalVarNumbers > 0)
	{
		numRamWrites++;
		totalVarNumbers -= 128; // Can read/write 256 bits at a time from RAM
	}

	/*cout << "varNodeVals: ";
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < rowOnes[i]; j++)
			cout << varNodeVals[i][j] << " ";
		cout << "\n             ";
	}
	cout << endl;*/

	for(numIterations = 0; numIterations < totalIterations; numIterations++)
	{
		if(calculateSyndrome()) // no complexity for syndrome, since that will be identical to SPA
			break;

		// calculate the parity check equations (i.e. extrinsic information)
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
			{
				parityChecks[i][j] = 0;
				numLogic++;
				for(int l = 0; l < rowOnes[i]; l++)
				{
					if(rowSupport[i][j] != rowSupport[i][l])
					{
						parityChecks[i][j] ^= varNodeVals[i][l];
						numLogic++;
					}
					numLogic++;
				}
			}
		}

		// Cache paging for H matrix
		totalHNumbers = (rows*6)*16; // Storing ints as short ints (16 bits), assume no bits stored in Regsiters, must read everything
		while(totalHNumbers > 0)
		{
			numCacheReads++;
			totalHNumbers -= 64; // We can read/write 64 bits at a time from Cache to Registers.
		}

		// RAM paging here for reading H matrix
		totalHNumbers = (rows*6)*16 - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalHNumbers > 0)
		{
			numRamReads++;
			totalHNumbers -= 128; // Read 256 bits at a time from RAM
		}

		// Read varNodeVals, row*6 bits
		// Write parityChecks, row*6 bits
		// Only need one instance of paging here
		// Cache paging
		totalVarNumbers = (rows*6); // Read parityChecks and varNodeVals, store parityChecks
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6) - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			numRamWrites++;
			totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
		}

		/*cout << "parityChecks: ";
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
				cout << parityChecks[i][j] << " ";
			cout << "\n              ";
		}
		cout << endl;*/

		// Check the check node information for each variable node
		int value = 0;

		for(int i = 0; i < cols; i++) // Iterate through the columns (i.e. variable nodes)
		{
			for(int j = 0; j < colOnes[i]; j++) // Iterate through the one values in the columns
			{
				value = 100; // set to an initial value, as an indication that this hasn't been set yet.
				numLogic++;

				for(int l = 0; l < colOnes[i]; l++)
				{
					numLogic++;
					if(j != l) // Make sure we're only sending information from all OTHER check nodes
					{
						int columnIndex = 0;
						// populate a new vector with all the values, then check if all values are equal
						for(int t = 0; t < rowOnes[ colSupport[i][l] ]; t++) // initially set value to paritycheck[first 1 in column][column we're looking at]
						{
							numLogic++;
							if (rowSupport[ colSupport[i][l] ][t] == i)
							{
								columnIndex = t;
								break;
							}
						}

						numLogic+=2;
						if(value == 100)
							value = parityChecks[ colSupport[i][l] ][columnIndex]; // set value if it hasn't been set before
						else if(value != parityChecks[ colSupport[i][l] ][columnIndex])
						{
							value = initHardDecision[i];
							break;
						}
					}
				}

				// Update the variable node value
				for(int t = 0; t < rowOnes[ colSupport[i][j] ]; t++)
				{
					numLogic++;
					if(rowSupport[ colSupport[i][j]][t] == i)
					{
						varNodeVals[ colSupport[i][j] ][t] = value;
						numLogic++;
						break;
					}
				}
			}
		}

		// Cache Paging here for reading H matrix
		totalHNumbers = (rows*6)*16; // Storing ints as short ints (16 bits), assume no bits stored in Regsiters, must read everything
		while(totalHNumbers > 0)
		{
			numCacheReads++;
			totalHNumbers -= 64; // We can read 64 bits at a time from Cache to Registers.
		}
		// RAM paging here for reading H matrix
		totalHNumbers = (rows*6)*16 - (16*1024); // Storing ints as short ints (16 bits), assume 16 kB in cache, the rest is stored in RAM
		while(totalHNumbers > 0)
		{
			numRamReads++;
			totalHNumbers -= 128; // We can read 16 kB at a time from Ram to Cache.(This probably needs to change, not sure exactly to what)
		}

		// Read parityChecks, rows*6 bits
		// Read initHardDecision, cols bits
		// Cache paging
		totalVarNumbers = (rows*6) + cols;
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6) + cols - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
		}

		// Write varNodeVals
		// Cache paging
		totalVarNumbers = (rows*6);
		while(totalVarNumbers > 0)
		{
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6) - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamWrites++;
			totalVarNumbers -= 128;
		}

		/*cout << "varNodeVals: ";
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
				cout << varNodeVals[i][j] << " ";
			cout << "\n             ";
		}
		cout << endl;*/

		// hard decision here

		int numOnes = 0;
		int numZeros = 0;

		for(int i = 0; i < cols; i++)
		{
			numOnes = 0;
			numZeros = 0;

			for(int j = 0; j < colOnes[i]; j++)
			{
				for(int t = 0; t < rowOnes[ colSupport[i][j] ]; t++)
				{
					numLogic++;
					if( (rowSupport[ colSupport[i][j] ][t]) == i )
					{
						numLogic++;
						if(varNodeVals[ colSupport[i][j] ][t] == 0)
							numZeros++;
						else
							numOnes++;
						numAdds++;
						break;
					}
				}
			}

			numLogic++;
			if( (colOnes[i] % 2) == 0 ) // if the degree of the variable node is even, also take into account the initial channel information
			{
				numLogic++;
				if(initHardDecision[i] == 0)
					numZeros++;
				else
					numOnes++;
				numAdds++;
			}

			hardDecisionVector[i] = numZeros > numOnes? 0 : 1;
			numLogic++;
		}

		/*cout << std::setw(15) << "hardDecision: ";
		for(int i = 0; i < cols; i++)
			cout << std::setw(3) << hardDecisionVector[i] << " ";
		cout << endl;*/

		// Cache Paging here for reading H matrix
		totalHNumbers = (rows*6)*16; // Storing ints as short ints (16 bits), assume no bits stored in Regsiters, must read everything
		while(totalHNumbers > 0)
		{
			numCacheReads++;
			totalHNumbers -= 64; // We can read 64 bits at a time from Cache to Registers.
		}
		// RAM paging here for reading H matrix
		totalHNumbers = (rows*6)*16 - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalHNumbers > 0)
		{
			numRamReads++;
			totalHNumbers -= 128;
		}

		// Read varNodeVals, rows*6 bits
		// Read initHardDecision, cols bits
		// Cache paging
		totalVarNumbers = (rows*6) + cols;
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = (rows*6) + cols - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			totalVarNumbers -= 128;
		}

		// Write hardDecisionVector, cols bits
		// Cache paging
		totalVarNumbers = cols;
		while(totalVarNumbers > 0)
		{
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = cols - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamWrites++;
			totalVarNumbers -= 128;
		}
	}

	for(int i = 0; i < cols; i++)
		estimatedCodeWord[i] = -2*hardDecisionVector[i]+1;

	return;
}

void GallagerADecoder::DecodeWithErasures(double* receivedVectorPtr, int* &estimatedCodeWord, int k)
{
	// In the erasure case, the message alphabet is {1, 0, -1}
	// First, use the Peeling Decoder to fill in the missing (erased) bits
	numErasedBits = 0;
	bool erased = true;

	// Fill in the hard decision vector with erasures
	for(int i = 0; i < cols; i++)
	{
		if(receivedVectorPtr[i] > 0)
			peelingVector[i] = 1.0;
		else if(receivedVectorPtr[i] < 0)
			peelingVector[i] = -1.0;
		else if(receivedVectorPtr[i] == 0) // check if the bit is erased
			peelingVector[i] = 0.0;
	}

	// Write peelingVector, cols*2 bits (signed 1, -1, or 0)
	// Cache paging
	totalVarNumbers = cols*2;
	while(totalVarNumbers > 0)
	{
		numCacheWrites++;
		totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
	}
	// RAM paging
	totalVarNumbers = cols*2 - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
	while(totalVarNumbers > 0)
	{
		numRamWrites++;
		totalVarNumbers -= 128;
	}

	while(erased)
	{
		// Fill in the erased bits vector first (will only be in parity bits, this code is in systematic form)
		numErasedBits = 0;

		for(int i = 0; i < cols - k; i++)
		{
			numLogic++;
			if(peelingVector[i] == 0.0) // check if the bit is erased
			{
				numErasedBits++;
				numAdds++;
				int* temp = new int[numErasedBits];

				for(int j = 0; j < numErasedBits-1; j++)
					temp[j] = erasedBits[j];
				temp[numErasedBits-1] = i;

				delete [] erasedBits;
				erasedBits = temp;
			}
		}

		// Read peelingVector, cols*2 bits
		// Cache paging
		totalVarNumbers = cols*2;
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = cols*2 - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			totalVarNumbers -= 128;
		}

		// Write erasedBits, numErasedBits*16 bits (short ints)
		// Write numErasedBits, 16 bits (1 short int)
		// Cache paging
		totalVarNumbers = numErasedBits*16 + 16;
		while(totalVarNumbers > 0)
		{
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = numErasedBits*16 + 16 - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamWrites++;
			totalVarNumbers -= 128;
		}

		/*cout << std::setw(15) << "hardDecision: ";
		for(int i = 0; i < cols; i++)
			cout << std::setw(3) << peelingVector[i] << " ";
		cout << endl;

		cout << "Erased Bits: ";
		for(int i = 0; i < numErasedBits; i++)
			cout << erasedBits[i] << " ";
		cout << endl;*/

		// Create a matrix to store the bit guesses
		peelingGuesses = new int*[numErasedBits]();

		for(int i = 0; i < numErasedBits; i++)
			peelingGuesses[i] = new int[ 3 ](); // change this from hard coding, basically increase the positive, negative or erasure counters

		for(int i = 0; i < numErasedBits; i++) // for each erased bit
		{
			for(int j = 0; j < colOnes[ erasedBits[i] ]; j++) // find the index of the parity check equations for each erased bit
			{
				int product = 1;

				for(int k = 0; k < rowOnes[ colSupport[ erasedBits[i] ][j] ]; k++) // find the ones in the parity check equation
				{
					numLogic++;
					if(rowSupport[ colSupport[ erasedBits[i] ][j] ][k] != erasedBits[i])
					{
						product *= peelingVector[ rowSupport[ colSupport[ erasedBits[i] ][j] ][k] ];
						numMults++;
					}

					numLogic++;
					if(product == 0)
						break;
				}

				numLogic++;
				if(product > 0)
					peelingGuesses[i][POSITIVE]++;
				else if(product < 0)
				{
					numLogic++;
					peelingGuesses[i][NEGATIVE]++;
				}
				else
				{
					numLogic++;
					peelingGuesses[i][ERASURE]++;
				}
				numAdds++;
			}
		}

		// Cache Paging here for reading H matrix
		totalHNumbers = (rows*6)*16; // Storing ints as short ints (16 bits), assume no bits stored in Regsiters, must read everything
		while(totalHNumbers > 0)
		{
			numCacheReads++;
			totalHNumbers -= 64; // We can read 64 bits at a time from Cache to Registers.
		}
		// RAM paging here for reading H matrix
		totalHNumbers = (rows*6)*16 - (16*1024); // 16 kB stored in Cache, so remove from RAM calculation
		while(totalHNumbers > 0)
		{
			numRamReads++;
			totalHNumbers -= 128; // We can read 16 kB at a time from Ram to Cache.(This probably needs to change, not sure exactly to what)
		}

		// Read erasedBits, numErasedBits*16 bits (short ints)
		// Read peelingVector, cols*2 bits
		// Read peelingGuesses, numErasedBits*3*16 (3 short ints for each erased bit)
		// Cache paging
		totalVarNumbers = numErasedBits*16 + cols*2 + numErasedBits*3*16;
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = numErasedBits*16 + cols*2 + numErasedBits*3*16 - 16*1024; // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
		}

		// Write peelingGuesses, numErasedBits*3*16 bits (short ints)
		// Cache paging
		totalVarNumbers = numErasedBits*3*16;
		while(totalVarNumbers > 0)
		{
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = numErasedBits*3*16 - 16*1024; // 16 kB stored in Cache, so remove from RAM calculation
		while(totalVarNumbers > 0)
		{
			numRamWrites++;
			totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
		}

		/*cout << "peeling: \n";
		for(int i = 0; i < numErasedBits; i++)
		{
			for(int j = 0; j < colOnes[ erasedBits[i] ]; j++)
				cout << peelingGuesses[i][j] << " ";
			cout << endl;
		}
		cout << endl;*/

		// Make the guesses here based on the peeling matrix
		bitsPutBack = 0;

		for(int i = 0; i < numErasedBits; i++)
		{
			positives = peelingGuesses[i][POSITIVE];
			negatives = peelingGuesses[i][NEGATIVE];
			erasures  = peelingGuesses[i][ERASURE];

			// This decision criteria can probably be tweaked and played with
			// If the number of positives is greater than the number

			numLogic++;
			if(positives > negatives/*&& !negatives*/)
			{
				peelingVector[ erasedBits[i] ] = 1.0;
				bitsPutBack++;
				numLogic+=2;
			}
			else if(negatives > positives/*&& !positives*/)
			{
				peelingVector[ erasedBits[i] ] = -1.0;
				bitsPutBack++;
				numLogic++;;
			}
		}

		// Read peelingGuesses, numErasedBits*3*16 bits
		// Read erasedBits, numErasedBits*16 bits
		// Cache paging
		totalVarNumbers = numErasedBits*3*16 + numErasedBits*16;
		while(totalVarNumbers > 0)
		{
			numCacheReads++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = numErasedBits*3*16 + numErasedBits*16 - 16*1024;
		while(totalVarNumbers > 0)
		{
			numRamReads++;
			totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
		}

		// Write peelingVector, cols*2 bits
		// Cache paging
		totalVarNumbers = cols*2;
		while(totalVarNumbers > 0)
		{
			numCacheWrites++;
			totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
		}
		// RAM paging
		totalVarNumbers = cols*2 - 16*1024;
		while(totalVarNumbers > 0)
		{
			numRamWrites++;
			totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
		}

		// If all of the bits have been erased, stop the loop, or if nothing was updated, end the loop
		if( (numErasedBits - bitsPutBack) == 0 || (bitsPutBack == 0) )
			erased = false;
		numLogic++;

		for(int i = 0; i < numErasedBits; i++)
			delete [] peelingGuesses[i];
		delete [] peelingGuesses;
	}

	/*cout << std::setw(15) << "hardDecision: ";
	for(int i = 0; i < cols; i++)
		cout << std::setw(3) << peelingVector[i] << " ";
	cout << endl;*/

	// Set any bits that have not been peeled to zero
	if( bitsPutBack == 0)
	{
		// Since no bits were put back, the erasedBits haven't changed. Fill them with 1's
		for(int i = 0; i < numErasedBits; i++)
			peelingVector[ erasedBits[i] ] = 1.0;

		/*cout << std::setw(15) << "hardDecision: ";
		for(int i = 0; i < cols; i++)
			cout << std::setw(3) << peelingVector[i] << " ";
		cout << endl;*/
	}
	numLogic++;

	// Read erasedBits, numErasedBits*16 bits
	// Cache paging
	totalVarNumbers = numErasedBits*16;
	while(totalVarNumbers > 0)
	{
		numCacheReads++;
		totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
	}
	// RAM paging
	totalVarNumbers = numErasedBits*16 - 16*1024;
	while(totalVarNumbers > 0)
	{
		numRamReads++;
		totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
	}

	// Write peelingVector, cols*2 bits
	// Cache paging
	totalVarNumbers = cols*2;
	while(totalVarNumbers > 0)
	{
		numCacheWrites++;
		totalVarNumbers -= 64; // Can read 64 bits at a time from Cache to Registers
	}
	// RAM paging
	totalVarNumbers = cols*2 - 16*1024;
	while(totalVarNumbers > 0)
	{
		numRamWrites++;
		totalVarNumbers -= 128; // Not sure how many bits Cache can read from RAM yet, probably need to change this.
	}

	Decode(peelingVector, estimatedCodeWord);

	return;
}

void GallagerADecoder::GetComplexity(double &mults, double &logic, double &adds, double &cacheReads, double &cacheWrites, double &ramReads, double &ramWrites)
{
	mults = numMults;
	logic = numLogic;
	adds  = numAdds;
	cacheReads = numCacheReads;
	cacheWrites = numCacheWrites;
	ramReads = numRamReads;
	ramWrites = numRamWrites;
}

void GallagerADecoder::ResetComplexity()
{
	numMults = 0.0;
	numLogic = 0.0;
	numAdds  = 0.0;
	numCacheReads  = 0.0;
	numCacheWrites = 0.0;
	numRamReads = 0.0;
	numRamWrites = 0.0;
}

int GallagerADecoder::GetNumIterations()
{
	return numIterations;
}

bool GallagerADecoder::calculateSyndrome()
{
	bool syndromeCheck = true;

	// The syndrome is basically the xor (mod 2 addition) of
	// the bits codeword with support in each row in H
	for(int i = 0; i < rows; i++)
	{
		syndrome[i] = 0; // reset the bit value

		for(int j = 0; j < rowOnes[i]; j++)
		{
			syndrome[i] ^= hardDecisionVector[ rowSupport[i][j] ];
		}

		if(syndrome[i] == 1)
			syndromeCheck = false;
	}

	return syndromeCheck;
}

////////////////////////////////////////
// Bit Flip Algorithm Class Functions //
////////////////////////////////////////

BitFlipDecoder::BitFlipDecoder(umat* Hptr, int maxNumIterations)
{
	// Max number of iterations for the sum product algorithm
	totalIterations = maxNumIterations;

	// Parity matrix dimensions
	rows = (*Hptr).n_rows;
	cols = (*Hptr).n_cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	// Allocate the memory for the rows in H
	for(int i = 0; i < rows; i++)
	{
		uvec ones = find((*Hptr).row(i));
		rowSupport[i] = new int[ones.size()]();
		rowOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			rowSupport[i][j] = ones[j];
	}

	// Allocate the memory for the columns in  H
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find((*Hptr).col(i));
		colSupport[i] = new  int[ones.size()]();
		colOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}

	failedChecks = new int[cols]();

	hardReceivedVector = new int[cols]();
	syndrome = new int[rows]();

	ResetComplexity();
}

BitFlipDecoder::~BitFlipDecoder()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
		delete [] rowSupport[i];
	delete [] rowSupport;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];
	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	delete [] failedChecks;

	delete [] hardReceivedVector;
	delete [] syndrome;
}

void BitFlipDecoder::Decode(double* receivedVectorPtr, int* &estimatedCodeWord)
{
	pSoftReceivedVector = receivedVectorPtr;

	numMems++;

	for(int i = 0; i < cols; i++)
	{
		hardReceivedVector[i] = pSoftReceivedVector[i] < 0;
		numMems += 3;
		numLogic++;
	}

	for(iterations = 0; iterations < totalIterations; iterations++)
	{
		numMems++;
		if( calculateSyndrome() )
			break;

		// If the syndrome is non-zero, find the  number of failed check sums for each bit
		for(int i = 0; i < cols; i++)
		{
			failedChecks[i] = 0;
			numMems += 2;
		}

		parityCheck = 0;
		numMems++;

		for(int row = 0; row < rows; row++)
		{
			numMems++;
			for(int col = 0; col < rowOnes[row]; col++)
			{
				parityCheck = 0;
				numMems += 2;

				for(int parCheck = 0; parCheck < rowOnes[row]; parCheck++)
				{
					numMems++;
					if(parCheck != col)
					{
						parityCheck ^= hardReceivedVector[ rowSupport[row][ parCheck] ];
						numMems += 3;
						numLogic++;
					}
				}

				if(parityCheck != (int)hardReceivedVector[ rowSupport[row][col] ])
				{
					failedChecks[ rowSupport[row][col] ]++;
					numMems += 2;
					numAdds++;
				}
				numMems += 3;
			}
		}

		// Find the maximum number of failed checks, and flip all bits in that set
		for(int i = 0; i < cols; i++)
		{
			if(i == 0 || maxFailedChecks < failedChecks[i])
			{
				maxFailedChecks = failedChecks[i];
				numMems += 2;
			}
			numMems += 3;
		}


		for(int i = 0; i < cols; i++)
		{
			if(failedChecks[i] == maxFailedChecks)
			{
				hardReceivedVector[i] = 1 - hardReceivedVector[i];
				numMems += 2;
				numAdds++;
			}
			numMems += 3;
		}
	}

	for(int i = 0; i < cols; i++)
	{
		estimatedCodeWord[i] = -2*hardReceivedVector[i]+1;
		numMems += 3;
		numAdds++;
		numMults++;
	}
	return;
}

void BitFlipDecoder::GetComplexity(double &mults, double &logic, double &adds, double &mems)
{
	mults = numMults;
	logic = numLogic;
	adds  = numAdds;
	mems  = numMems;
}

void BitFlipDecoder::ResetComplexity()
{
	numMults = 0.0;
	numLogic = 0.0;
	numAdds  = 0.0;
	numMems  = 0.0;
}

int BitFlipDecoder::GetNumIterations()
{
	return iterations;
}

bool BitFlipDecoder::calculateSyndrome()
{
	bool syndromeCheck = true;

	// The syndrome is basically the xor (mod 2 addition) of
	// the bits codeword with support in each row in H
	for(int i = 0; i < rows; i++)
	{
		syndrome[i] = 0; // reset the bit value
		numMems += 2;

		for(int j = 0; j < rowOnes[i]; j++)
		{
			syndrome[i] ^= hardReceivedVector[ rowSupport[i][j] ];
			numMems += 4;
			numLogic++;
		}

		if(syndrome[i] == 1)
		{
			syndromeCheck = false;
			numMems++;
		}
		numMems++;
	}

	return syndromeCheck;
}

// Timing functions

double get_wall_time(){
	struct timeval time;

	if(gettimeofday(&time,NULL)){
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}
