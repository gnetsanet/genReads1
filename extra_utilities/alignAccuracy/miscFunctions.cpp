#include "miscFunctions.h"
#include "StringRef.h"

#include <vector>
#include <string>

using namespace std;

vector<StringRef> split3( string const& str, char delimiter = ' ' )
{
	vector<StringRef>   result;

	enum State { inSpace, inToken };

	State state = inSpace;
	char const*     pTokenBegin = 0;    // Init to satisfy compiler.
	for( string::const_iterator it = str.begin(); it != str.end(); ++it )
	{
		State const newState = (*it == delimiter? inSpace : inToken);
		if( newState != state )
		{
			switch( newState )
			{
			case inSpace:
				result.push_back( StringRef( pTokenBegin, &*it - pTokenBegin ) );
				break;
			case inToken:
				pTokenBegin = &*it;
			}
		}
		state = newState;
	}
	if( state == inToken )
	{
		result.push_back( StringRef( pTokenBegin, &*str.rbegin() - pTokenBegin + 1 ) );
	}
	return result;
}

void defaultHumanChr( countMap &refLens )
{
	refLens.insert(make_pair("chr1",249250621));
	refLens.insert(make_pair("chr2",243199373));
	refLens.insert(make_pair("chr3",198022430));
	refLens.insert(make_pair("chr4",191154276));
	refLens.insert(make_pair("chr5",180915260));
	refLens.insert(make_pair("chr6",171115067));
	refLens.insert(make_pair("chr7",159138663));
	refLens.insert(make_pair("chr8",146364022));
	refLens.insert(make_pair("chr9",141213431));
	refLens.insert(make_pair("chr10",135534747));
	refLens.insert(make_pair("chr11",135006516));
	refLens.insert(make_pair("chr12",133851895));
	refLens.insert(make_pair("chr13",115169878));
	refLens.insert(make_pair("chr14",107349540));
	refLens.insert(make_pair("chr15",102531392));
	refLens.insert(make_pair("chr16",90354753));
	refLens.insert(make_pair("chr17",81195210));
	refLens.insert(make_pair("chr18",78077248));
	refLens.insert(make_pair("chr19",59128983));
	refLens.insert(make_pair("chr20",63025520));
	refLens.insert(make_pair("chr21",48129895));
	refLens.insert(make_pair("chr22",51304566));
	refLens.insert(make_pair("chrX",155270560));
	refLens.insert(make_pair("chrY",59373566));
	refLens.insert(make_pair("chrM",16569));
	refLens.insert(make_pair("chrGL000207.1",4262));
	refLens.insert(make_pair("chrGL000226.1",15008));
	refLens.insert(make_pair("chrGL000229.1",19913));
	refLens.insert(make_pair("chrGL000231.1",27386));
	refLens.insert(make_pair("chrGL000210.1",27682));
	refLens.insert(make_pair("chrGL000239.1",33824));
	refLens.insert(make_pair("chrGL000235.1",34474));
	refLens.insert(make_pair("chrGL000201.1",36148));
	refLens.insert(make_pair("chrGL000247.1",36422));
	refLens.insert(make_pair("chrGL000245.1",36651));
	refLens.insert(make_pair("chrGL000197.1",37175));
	refLens.insert(make_pair("chrGL000203.1",37498));
	refLens.insert(make_pair("chrGL000246.1",38154));
	refLens.insert(make_pair("chrGL000249.1",38502));
	refLens.insert(make_pair("chrGL000196.1",38914));
	refLens.insert(make_pair("chrGL000248.1",39786));
	refLens.insert(make_pair("chrGL000244.1",39929));
	refLens.insert(make_pair("chrGL000238.1",39939));
	refLens.insert(make_pair("chrGL000202.1",40103));
	refLens.insert(make_pair("chrGL000234.1",40531));
	refLens.insert(make_pair("chrGL000232.1",40652));
	refLens.insert(make_pair("chrGL000206.1",41001));
	refLens.insert(make_pair("chrGL000240.1",41933));
	refLens.insert(make_pair("chrGL000236.1",41934));
	refLens.insert(make_pair("chrGL000241.1",42152));
	refLens.insert(make_pair("chrGL000243.1",43341));
	refLens.insert(make_pair("chrGL000242.1",43523));
	refLens.insert(make_pair("chrGL000230.1",43691));
	refLens.insert(make_pair("chrGL000237.1",45867));
	refLens.insert(make_pair("chrGL000233.1",45941));
	refLens.insert(make_pair("chrGL000204.1",81310));
	refLens.insert(make_pair("chrGL000198.1",90085));
	refLens.insert(make_pair("chrGL000208.1",92689));
	refLens.insert(make_pair("chrGL000191.1",106433));
	refLens.insert(make_pair("chrGL000227.1",128374));
	refLens.insert(make_pair("chrGL000228.1",129120));
	refLens.insert(make_pair("chrGL000214.1",137718));
	refLens.insert(make_pair("chrGL000221.1",155397));
	refLens.insert(make_pair("chrGL000209.1",159169));
	refLens.insert(make_pair("chrGL000218.1",161147));
	refLens.insert(make_pair("chrGL000220.1",161802));
	refLens.insert(make_pair("chrGL000213.1",164239));
	refLens.insert(make_pair("chrGL000211.1",166566));
	refLens.insert(make_pair("chrGL000199.1",169874));
	refLens.insert(make_pair("chrGL000217.1",172149));
	refLens.insert(make_pair("chrGL000216.1",172294));
	refLens.insert(make_pair("chrGL000215.1",172545));
	refLens.insert(make_pair("chrGL000205.1",174588));
	refLens.insert(make_pair("chrGL000219.1",179198));
	refLens.insert(make_pair("chrGL000224.1",179693));
	refLens.insert(make_pair("chrGL000223.1",180455));
	refLens.insert(make_pair("chrGL000195.1",182896));
	refLens.insert(make_pair("chrGL000212.1",186858));
	refLens.insert(make_pair("chrGL000222.1",186861));
	refLens.insert(make_pair("chrGL000200.1",187035));
	refLens.insert(make_pair("chrGL000193.1",189789));
	refLens.insert(make_pair("chrGL000194.1",191469));
	refLens.insert(make_pair("chrGL000225.1",211173));
	refLens.insert(make_pair("chrGL000192.1",547496));
}