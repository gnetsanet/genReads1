#include "StringRef.h"
#include "miscFunctions.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include <unistd.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>		// this isn't required on my comp, but is needed for compiling on RCF

#include <string>
#include <vector>
#include <tr1/unordered_map>

using namespace std;

const char MAGIC[4] = {'B','A','M','\1'};
const unsigned int MAGIC_LEN = 4;

// lookup table describing number of Ns present in each of the possible 8bit representations of 2 adjacent SEQ characters
const int numberOfNs[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2};

const string seqChar("=ACMGRSVTWYHKDBN");

const int MAX_LEN = 4096;		// maximum read length we expect to encounter
const int MAX_TAG_STR = 4096;	// maximum length we expect string elements of alignment tags to be
const int MAX_MAPQ = 4096;		// maximum value we expect mapping-quality to be
const int MAX_IMD = 10000;		// maximum value we expect inner mate distance to be (discard if higher)

// Only roughly estimate coverage (ignores structural variation specified in cigar strings)
const bool ROUGH_COV = true;
// What is a reasonable estimate for the maximum coverage we're likely to see in an alignment?
const unsigned int INIT_MAXCOV = 1000;

// within this many nucleotides of origin is still "correct"
const int WITHIN_OK = 10;


int getI()
{
	char str[4];
	cin.read(str,4);
	int i = *(int*)str;		// evil bit-hackery, dangerous! change me later.
	return i;				// also I think this assumes machine is little-endian...
}

unsigned int getUI()
{
	char str[4];
	cin.read(str,4);
	unsigned int i = *(unsigned int*)str;
	return i;
}

string getString(const unsigned int length)
{
	char * buff = new char [length];
	cin.read(buff,length);
	string str(buff);
	free(buff);
	return str;
}

string getStringUntilDelim(const char delim)
{
	char buff = cin.get();
	string str;
	while (buff != delim)
	{
		str += buff;
		buff = cin.get();
	}
	return str;
}

bool str_is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

void printMap(countMap& m)
{
	vector<string> keys;
	for ( countMap::iterator it = m.begin(); it != m.end(); ++it )
	{
		keys.push_back(it->first);
	}
	sort(keys.begin(),keys.end());
	for (unsigned int i=0; i<keys.size(); i++)
		cout << keys[i] << " : " << m[keys[i]] << endl;
}

int main(int argc, char *argv[])
{
	int c;
	char *outFile = NULL;

	while ((c = getopt(argc, argv, "o:")) != -1)
	{
		switch (c)
		{
			case 'o':
				outFile = optarg;
				break;
			case '?':
			default:
				fprintf(stderr, "Usage: %s [-o outFile]\n",argv[0]);
				exit(1);
		}
	}

	bool needReadName = true;//false;
	bool needCigar    = true;
	bool needSeq      = false;
	bool needQual     = false;
	bool needTags     = true;

	bool printNHWarning = true;	// only print "NH tag not found" warning once

	cin.sync_with_stdio(false);

	bool getHeaderInfo = true;
	unsigned long int nRecords = 0;
	unsigned long int uniqueReads = 0;
	unsigned long int nUnmapped = 0;
	unsigned long int nIncorrectlyMapped = 0;
	unsigned long int nIncorrectAndMulti = 0;

	countMap refLengthByReference;
	double multiMappedReadCount = 0.0;
	int nReferences = 0;
	int prevPos  = 0;
	int prevLseq = 0;
	vector<string> refList;
	vector<unsigned long int> readCountByReference;
	vector<unsigned long int> baseCountByReference;
	vector<int> N_Frequency(MAX_LEN,0);
	vector<unsigned long int> readLengthDistribution(MAX_LEN,0);
	vector<unsigned long int> IMD_Distribution(2*MAX_IMD,0);
	vector<unsigned short int> perBaseCoverage;
	tr1::unordered_map< string, vector<unsigned long int> > allCoverages;	// by Ref
	tr1::unordered_map< string, vector<unsigned long int> > allMappingQual;	// by Ref
	string sortingType;
	string myAligner;
	string myVersion;
	string prevRef;

	while (cin)
	{
		// first time through? let's grab header info
		if (getHeaderInfo)
		{
			getHeaderInfo = false;

			string magic = getString(MAGIC_LEN);
			for (unsigned int i=0; i<MAGIC_LEN; i++)
			{
				if (magic[i] != MAGIC[i])
				{
					cerr << "Error[0]: Magic string not found" << endl;
					exit(1);
				}
			}

			int l_header  = getI();
			string header = getString(l_header);

			//cout << header << endl;
			vector<StringRef> const header_splt = split3( header, '\n' );
			for (unsigned int i=0; i<header_splt.size(); i++)
			{
				string header_line(header_splt[i].begin(),header_splt[i].end());
				if (header_line.size() < 5)
				{
					cerr << "Warning[1]: Incomplete line in header." << endl;
					continue;
				}
				vector<StringRef> const line_splt = split3( header_line, '\t' );
				//cout << i << " " << header_line << endl;
				string tag(line_splt[0].begin(),line_splt[0].end());


				if (tag == "@HD")
				{
					bool VN_ok = false;
					bool SO_ok = true;
					for (unsigned int j=1; j<line_splt.size(); j++)
					{
						string tag2(line_splt[j].begin(), line_splt[j].begin()+3);
						string val2(line_splt[j].begin()+3, line_splt[j].end());	// I have no idea why I need this +1, hmmm...

						if (tag2 == "VN:")
							VN_ok = true;
						else if (tag2 == "SO:")
						{
							if (val2 == "unknown" || val2 == "unsorted" || val2 == "queryname" || val2 == "coordinate")
								sortingType = val2;
							else
								SO_ok = false;
						}
						else if (tag2 == "GO:");	// other misc tags that are ok.
						else
						{
							cerr << "Error[2]: Unknown tag '" << tag2 << "' found in header line: " << header_line << endl;
							exit(1);
						}
					}
					if (!VN_ok)
					{
						cerr << "Error[3]: 'VN' tag not found in header" << endl;
						exit(1);
					}
					if (!SO_ok)
					{
						cerr << "Error[4]: Unknown value for 'SO' tag in header line: " << header_line << endl;
						exit(1);
					}
				}


				else if (tag == "@SQ")
				{
					string SN_ok;
					unsigned long int LN_ok = 0;
					for (unsigned int j=1; j<line_splt.size(); j++)
					{
						string tag2(line_splt[j].begin(), line_splt[j].begin()+3);
						string val2(line_splt[j].begin()+3, line_splt[j].end());

						if (tag2 == "SN:")
							SN_ok.assign(val2);
						if (tag2 == "LN:" && str_is_number(val2))
							LN_ok = strtoul(val2.c_str(), NULL, 0);
					}
					if (SN_ok != "" and LN_ok > 0)
					{
						if (refLengthByReference.find(SN_ok) != refLengthByReference.end())
						{
							cerr << "Error[5]: Duplicate 'SN' tags: " << SN_ok << endl;
							exit(1);
						}
						refLengthByReference[SN_ok] = LN_ok;
						nReferences += 1;
						refList.push_back(SN_ok);
						readCountByReference.push_back(0);
						baseCountByReference.push_back(0);
						vector<unsigned long int> temp(MAX_MAPQ,0);
						allMappingQual[SN_ok] = temp;
					}
					else
					{
						cerr << "Error[6]: 'SN' and/or 'LN' tag not correct in header line: " << header_line << endl;
						exit(1);
					}
				}


				else if (tag == "@RG")
				{
					//
				}


				else if (tag == "@PG")
				{
					for (unsigned int j=1; j<line_splt.size(); j++)
					{
						string tag2(line_splt[j].begin(), line_splt[j].begin()+3);
						string val2(line_splt[j].begin()+3, line_splt[j].end());

						if (tag2 == "ID:")
						{
							transform(val2.begin(), val2.end(), val2.begin(), ::tolower);
							myAligner.assign(val2);
						}
						else if (tag2 == "VN:")
							myVersion.assign(val2);
					}
				}


				else if (tag == "@CO")
				{
					//
				}


				else
				{
					cerr << "Error[7]: Unknown tag in header: " << header_line << endl;
					exit(1);
				}
			}

			int n_ref = getI();
			for (int i=0; i<n_ref; i++)
			{
				int l_refName  = getI();
				string refName = getString(l_refName);
				int refLen     = getI();

				// let's make it unsigned to get rid of annoying compiler warnings
				unsigned int refLen_ui = 0;
				if (refLen > 0)
					refLen_ui = (unsigned int)refLen;
				else
				{
					cerr << "Error[8]: Reference length in dictionary is negative." << endl;
					exit(1);
				}

			}
		}

		if (sortingType != "coordinate")
		{
			cerr << "Error! Please sort your bam by coordinate, thanks!" << endl;
			exit(1);
		}


		// read in record
		int l_record  = getI();

		int refID = getI();
		int pos = getI();

		unsigned int bin_mq_nl = getUI();
		unsigned int bin = bin_mq_nl>>16;
		unsigned int mapQ = (bin_mq_nl-(bin<<16))>>8;
		unsigned int l_readName = bin_mq_nl-(bin<<16)-(mapQ<<8);

		unsigned int flag_nc = getUI();
		unsigned int samFlag = flag_nc>>16;
		unsigned int nCigarOp = flag_nc-(samFlag<<16);

		int l_seq = getI();
		prevLseq = l_seq;
		int next_refID = getI();
		int next_pos = getI();
		int tlen = getI();

		if (!cin)	// if we started reading over into the end of the stream, get out!
			break;

		if (l_seq >= MAX_LEN)
		{
			// we encountered reads longer than we anticipated, extend some vectors..
			while((unsigned int)l_seq > N_Frequency.size()+1)
			{
				N_Frequency.push_back(0);
				readLengthDistribution.push_back(0);
			}
		}
		readLengthDistribution[(unsigned int)l_seq] += 1;

		string readName;
		if (needReadName)
		{
			readName = getString(l_readName);
		}
		else
			cin.ignore(l_readName);

		//cout << "rawr: " << readName << ' ' << samFlag << ' ' << refID << ' ' << pos << endl;


		if (needCigar)
		{
			unsigned int rawCigar[nCigarOp];
			for (unsigned int i=0; i<nCigarOp; i++)
				rawCigar[i] = getUI();
		}
		else
			cin.ignore(4*nCigarOp);



		if (needSeq)
		{
			string rawSeq = getString((l_seq+1)/2);

			int N_count = 0;
			for (int i=0; i<(l_seq+1)/2; i++)
				N_count += numberOfNs[(unsigned char)rawSeq[i]];

			if (N_count > 0)
				N_Frequency[N_count] += 1;
		}
		else
			cin.ignore((l_seq+1)/2);



		
		if (needQual)
			string qual = getString(l_seq);
		else
			cin.ignore(l_seq);

		

		// parse alignment tags
		unsigned int NH_value = 0;

		int tagBytes = l_record-32-l_readName-4*nCigarOp-((l_seq+1)/2)-l_seq;
		if (needTags)
		{
			string allTagData = getString(tagBytes);
			for (int i=0; i<tagBytes-4; i++)
			{
				if (allTagData[i] == 'N' && allTagData[i+1] == 'H')
				{
					if (allTagData[i+2] == 'C')
						NH_value = (unsigned int)(allTagData[i+3]);
					else if (allTagData[i+2] == 'S')
						NH_value = (unsigned int)( (allTagData[i+3]<<8) + allTagData[i+4]);
					break;
				}
			}
		}
		else
			cin.ignore(tagBytes);

		//
		// parse SAM flag
		//
		bool isPaired     = (samFlag&1);
		bool isUnmapped   = (samFlag&4);
		bool mateUnmapped = (samFlag&8);
		bool isFirstRead  = (samFlag&64);
		bool isSecondRead = (samFlag&128);
		bool isNotPrimary = (samFlag&256);
		bool isSupplement = (samFlag&2048);


		//
		// make sure our reference ID is valid before we start using it for stuff
		//
		string refName;
		if (refID >= nReferences)
		{
			cerr << "Error[13]: Invalid reference ID in record " << readName << " (#" << nRecords << ")" << endl;
			exit(1);
		}
		else if (isUnmapped)
		{
			refName = "*";
		}
		else
			refName = refList[refID];


		//
		// did we move on to a new reference? If so, reinitialize coverage/pos variables
		//
		bool changedRef = false;
		if (!isUnmapped && refName != prevRef && sortingType == "coordinate")
		{
			if (prevRef != "")
			{
				vector<unsigned long int> coverageHistogram(INIT_MAXCOV,0);
				for (unsigned int i=0; i<perBaseCoverage.size(); i++)
				{
					unsigned short int cov = perBaseCoverage[i];
					if (cov <= coverageHistogram.size())
						coverageHistogram[cov] += 1;
					else
					{
						while (coverageHistogram.size() <= cov)
							coverageHistogram.push_back(0);
						coverageHistogram[cov] += 1;
					}
				}
				coverageHistogram[0] -= l_seq;
				allCoverages[prevRef] = coverageHistogram;
			}

			if (!isUnmapped)
			{
				perBaseCoverage.assign(refLengthByReference[refName]+l_seq,0);
			}

			//cout << prevRef << " --> " << refName << endl;
			prevRef = refName;
			prevPos = 0;
			changedRef = true;
		}


		//
		// try to determine is read is multi-mapped (incorporate aligner info here?)
		//
		bool isMultimapped = false;
		if (isUnmapped)
			nUnmapped += 1;
		else
		{
			string alignerPrefix = myAligner.substr(0,3);
			if (isNotPrimary || isSupplement)
				isMultimapped = true;
			else if (NH_value > 1)
				isMultimapped = true;
			// a special case for BWA:
			else if (alignerPrefix == "bwa" && mapQ <= 3)
				isMultimapped = true;
		}
		
		//
		// if multi-mapped, tally up "partial" reads. i.e. this read was mapped 3 times --> add 1/3 to read count
		//
		if (isMultimapped)
		{
			if (NH_value > 1 && printNHWarning)
			{
				// with too many reads the precision of double might start causing issues...
				multiMappedReadCount += 1.0/NH_value;
			}
			else
			{
				// I don't know what the multiplicity is, ohwell...
				if (printNHWarning)
				{
					cerr << "Warning[14]: Multi-mapped read multiplicity unknown. We're going to assume no NH tag means read was mapped uniquely, adjust your faith in the counted totals accordingly!" << endl;
					printNHWarning = false;
				}
				multiMappedReadCount += 1.0;
			}
		}
		else if (!isUnmapped)
			uniqueReads += 1;

		//
		// count # of reads and # bases that cover each reference
		//
		if (!isUnmapped)
		{
			readCountByReference[refID] += 1;
			baseCountByReference[refID] += l_seq;
		}

		//
		// count mapping quality value for each reference
		//
		if (!isUnmapped)
		{
			// expand vector if somehow we come across weird bigger mapping quality score
			while (mapQ >= allMappingQual[refName].size())
				allMappingQual[refName].push_back(0);
			allMappingQual[refName][mapQ] += 1;
		}

		//
		// count inner-mate distance if applicable
		//
		if (isPaired && !isUnmapped && !mateUnmapped && refID == next_refID)
		{
			//because this is just an estimation, let's assume mate read length is same as ours
			int imd = abs(next_pos - pos) - l_seq;
			if (abs(imd) < MAX_IMD)	// discard outliers
				IMD_Distribution[imd + MAX_IMD] += 1;
		}

		//
		// compute per-base coverage
		//
		if (!isUnmapped)
		{
			// rough coverage estimation (ignore cigar strings)
			if (ROUGH_COV)
			{
				if (pos < 0 || (unsigned int)pos+l_seq >= perBaseCoverage.size())
				{
					cerr << "Warning[15]: Read mapped out of bounds: " << readName << " (#" << nRecords << ")" << endl;
				}
				else
				{
					for (int i=0; i<l_seq; i++)
					{
						perBaseCoverage[i+pos] += 1;
					}
				}
			}
		}

		//
		// ensure sorting is correct
		//
		if (sortingType == "coordinate" && !isUnmapped && !changedRef && pos < prevPos)
		{
			cerr << "Error[16]: Record " << nRecords << " is not sorted properly" << endl;
			exit(1);
		}
		prevPos = pos;

		//
		// so, did the read map to the correct position?
		//
		bool isIncorrectlyMapped = false;
		if(!isUnmapped)
		{
			if (isPaired)
			{
				vector<StringRef> const line_splt = split3( readName, '-' );
				int correctPos = 0;
				if (isFirstRead)
				{
					string correctPosStr(line_splt[1].begin(),line_splt[1].end());
					correctPos = atoi(correctPosStr.c_str());
				}
				else if (isSecondRead)
				{
					string correctPosStr(line_splt[2].begin(),line_splt[2].end());
					correctPos = atoi(correctPosStr.c_str());
				}
				else
				{
					cerr << "Error[17]: Paired read is not flagged as first or second: " << readName << endl;
					exit(1);
				}

				if (abs(correctPos - pos - 1) <= WITHIN_OK)
				{
					//
				}
				else
				{
					isIncorrectlyMapped = true;
				}
			}
			else
			{
				vector<StringRef> const line_splt = split3( readName, '-' );
				int correctPos = 0;
				string correctPosStr(line_splt[1].begin(),line_splt[1].end());
				correctPos = atoi(correctPosStr.c_str());
				if (abs(correctPos - pos - 1) <= WITHIN_OK)
				{
					//
				}
				else
				{
					isIncorrectlyMapped = true;
				}
			}
		}
		if (isIncorrectlyMapped)
		{
			nIncorrectlyMapped += 1;
			if (isMultimapped)
			{
				nIncorrectAndMulti += 1;
			}
		}

		nRecords += 1;
		if (nRecords % 1000000 == 0)
			cout << nRecords << endl;

	}


	// if we ended in something other than unmapped, add coverage hitogram to allCoverages
	if (prevRef != "*")
	{
		vector<unsigned long int> coverageHistogram(INIT_MAXCOV,0);
		for (unsigned int i=0; i<perBaseCoverage.size(); i++)
		{
			unsigned short int cov = perBaseCoverage[i];
			if (cov <= coverageHistogram.size())
				coverageHistogram[cov] += 1;
			else
			{
				while (coverageHistogram.size() <= cov)
					coverageHistogram.push_back(0);
				coverageHistogram[cov] += 1;
			}
		}
		coverageHistogram[0] -= prevLseq;
		allCoverages[prevRef] = coverageHistogram;
	}

	long unsigned int multimapped_int = (long unsigned int)ceil(multiMappedReadCount);
	long unsigned int uniqueReadIDs   = (uniqueReads + nUnmapped + multimapped_int);
	float unmapPercent = 100.0 * nUnmapped / (float)nRecords;
	float incorrectPercent = 100.0 * nIncorrectlyMapped / (float)nRecords;
	float incorMultPercent = 100.0 * nIncorrectAndMulti / (float)nRecords;

	cout << "\n";
	cout << ">>RECORDS_READ\t\t" << nRecords << "\n";
	cout << ">>UNIQUELY_MAPPED\t" << uniqueReads << "\n";
	cout << ">>MULTI_MAPPED_IDS\t" << multimapped_int << "\n";
	cout << ">>UNMAPPED\t\t" << nUnmapped << "\t(" << unmapPercent << " %)\n";
	cout << ">>INCORRECTLY_MAPPED\t" << nIncorrectlyMapped << "\t(" << incorrectPercent << " %)\n";
	cout << ">>INCORRECT_AND_MULTI_MAPPED\t" << nIncorrectAndMulti << "\t(" << incorMultPercent << " %)\n";
	cout << ">>TOTAL_UNIQUE_IDS\t" << uniqueReadIDs << "\n";
	cout << "\n";

	//
	// WRITE OUTFILE
	//
	if (outFile != NULL)
	{
		ofstream outF;
		outF.open(outFile);
		
		//	readCount stats
		//
		outF << ">>ALIGNER\t" << myAligner << "\t" << myVersion << "\n";
		outF << ">>RECORDS_READ\t" << nRecords << "\n";
		outF << ">>UNIQUELY_MAPPED\t" << uniqueReads << "\n";
		outF << ">>MULTI_MAPPED_IDS\t" << multimapped_int << "\n";
		outF << ">>UNMAPPED\t" << nUnmapped << "\n";
		outF << ">>INCORRECTLY_MAPPED\t" << nIncorrectlyMapped << "\n";
		outF << ">>INCORRECT_AND_MULTI_MAPPED\t" << nIncorrectAndMulti << "\n";
		outF << ">>TOTAL_UNIQUE_IDS\t" << uniqueReadIDs << "\n";

		//	ref-Lengths
		//
		outF << ">>REF_LENGTHS\n";
		for (unsigned int i=0; i<refList.size(); i++)
			outF << refList[i] << "\t" << refLengthByReference[refList[i]] << "\n";

		//	readCount per reference
		//
		outF << ">>READ_COUNT_PER_REF\tnRecords\tnBases\n";
		for (unsigned int i=0; i<refList.size(); i++)
			outF << refList[i] << "\t" << readCountByReference[i] << "\t" << baseCountByReference[i] << "\n";

		//	read-length distribution
		outF << ">>READ_LENGTH\n";
		for (unsigned int i=0; i<readLengthDistribution.size(); i++)
		{
			if (readLengthDistribution[i] > 0)
				outF << i << "\t" << readLengthDistribution[i] << "\n";
		}

		//	N-frequency
		//
		outF << ">>N_FREQUENCY\n";
		for (unsigned int i=0; i<N_Frequency.size(); i++)
		{
			if (N_Frequency[i] > 0)
				outF << i << "\t" << N_Frequency[i] << "\n";
		}

		//	inner mate distance histogram
		//
		outF << ">>IMD_HISTOGRAM\n";
		for (int i=0; (unsigned int)i<IMD_Distribution.size(); i++)
		{
			if (IMD_Distribution[i] > 0)
				outF << i - MAX_IMD << "\t" << IMD_Distribution[i]/2 << "\n";
		}

		//	MAPQ histogram
		//
		outF << ">>MAPQ_HISTOGRAM\n";
		for (unsigned int i=0; i<refList.size(); i++)
		{
			string key = refList[i];
			outF << ">" << key << "\n";
			for (unsigned int j=0; j<allMappingQual[key].size(); j++)
			{
				if (allMappingQual[key][j] > 0)
					outF << j << "\t" << allMappingQual[key][j] << "\n";
			}
		}

		//	coverage histogram
		//
		outF << ">>COVERAGE_HISTOGRAM\n";
		for (unsigned int i=0; i<refList.size(); i++)
		{
			string key = refList[i];
			outF << ">" << key << "\n";
			for (unsigned int j=0; j<allCoverages[key].size(); j++)
			{
				if (allCoverages[key][j] > 0)
					outF << j << "\t" << allCoverages[key][j] << "\n";
			}
		}

		outF.close();
	}

	return 0;
}

