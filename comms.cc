/*
 * comms.cc
 *
 *  Created on: 2015-09-23
 *      Author: kassy
 */

#include <armadillo>
#include <fcntl.h>
#include <iostream>
#include <iomanip>


#include "comms.h"

long genRandomSeed()
{
  long randnum = 0;
  int fd = open ("/dev/urandom", O_RDONLY);
  if (fd != -1) {
    read( fd, &randnum, 4 );
    close(fd);
  }
  else{
  	printf("ERROR: Can't read /dev/urandom.\n");
  	exit(1);
  }

  if(randnum>0) randnum*=-1;

  if(!randnum){
  	printf("ERROR: Zero seed.\n");
  	exit(1);
  }

  return randnum;
}

void readParityFile(umat &H, umat &G, int &n, int &k)
{
    cout << "Please enter the parity check matrix binary filename >> ";

    // check for valid filename
    ifstream inputFile;
    checkForValidFilename(inputFile);

    // Read in n and k
    inputFile.read( (char*)&n, sizeof(int) );
    inputFile.read( (char*)&k ,sizeof(int) );

    // Read in H
    H.zeros(k,n);
    for( int i = 0; i < k; i++ )
        for( int j = 0; j < n; j++ )
            inputFile.read( (char*)&H(i,j), sizeof(char) );

    inputFile.read( (char*)&k,sizeof(int) );

    // Read in G
    G.zeros(k,n);
    for( int i = 0; i < k; i++ )
        for( int j = 0; j < n; j++ )
            inputFile.read( (char*)&G(i,j),sizeof(char) );

    inputFile.read( (char*)&k, sizeof(int) );
    inputFile.close();

    return;
}

void checkForValidFilename(ifstream &inputFile)
{
    string fileName;

    while(1)
    {
        cin >> fileName;
        inputFile.open( fileName.c_str(), ios::binary | ios::in );
        if( inputFile )
            break;
        cout << "You have entered an invalid filename. Please try again >> " << flush;
    }

    return;
}

void getUserDecodingScheme(bool &SPA, bool &GAA)
{
    cout << "\nWhich decoding scheme would you like to use; SPA, or GAA? >> ";

    string decodingScheme;

    SPA = false;
    GAA  = false;

    while(1)
    {
        cin >> decodingScheme;

        if( (decodingScheme.compare("SPA") == 0) )
            SPA = true;
        else if( decodingScheme.compare("GAA") == 0 )
            GAA = true;
        else
        {
            cout << "You have entered an invalid decoding scheme, please try again >> " << flush;
            continue;
        }

        return;
    }
}

void getSNRValues(vector<double> &dBSNR)
{
    cout << "\nPlease enter the SNRs that you wish to simulate. When you are done, enter D/d >> ";

    while(1)
    {
        string input;
        cin >> input;

        if( (input.compare("D") == 0) || (input.compare("d") == 0 ) )
            break;

        double snr = atof(input.c_str());

        // if atof fails, (i.e. no floating point number, it will return 0.0
        // check if the input was not 0.0 before putting out an error
        if( snr == 0.0 && !( (input.compare("0.0") == 0) || (input.compare("0") == 0) ) )
        {
            cout << "You have entered an invalid SNR value. Please try again >> ";
            continue;
        }

        dBSNR.push_back(snr);
    }

    // sort the SNR's, and remove any duplicates
    sort( dBSNR.begin(), dBSNR.end() );
    dBSNR.erase( unique( dBSNR.begin(), dBSNR.end() ), dBSNR.end() );

    return;
}

void getPuncturingFile(vector<int> &bitsToPuncture, int &totalBitsToPuncture)
{
    cout << "\nPlease enter the name of the puncturing file >> ";

    ifstream puncturingFile;
    checkForValidFilename(puncturingFile);

    // Read in the total number of punctured bits
    puncturingFile.read( (char*)&totalBitsToPuncture, sizeof(int) );
    bitsToPuncture.resize(totalBitsToPuncture);

    // Read in the punctured bits
    for( int i = 0; i < totalBitsToPuncture; i++ )
        puncturingFile.read( (char*)&bitsToPuncture[i], sizeof(int) );

    // Close the puncturing file
    puncturingFile.close();

    cout << "\nThere are " << totalBitsToPuncture << " bits that can be punctured. \n";

    return;
}

const char* createOutputFileName(int n, int k)
{
    string file;
    cout << "\nPlease enter the name of the file that you wish to save the output to >> ";

    cin >> file;
    const char* fileName = file.c_str();

    return fileName;
}

