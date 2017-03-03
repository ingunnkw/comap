#include <stoken.h>

using namespace std;

namespace skn
{

vector<string> tokenize(const string & line, const string & sep)
{
	using namespace std;
	vector<string> result;
	int i = line.find_first_not_of(sep);
	int j = line.find_first_of(sep, i);
	while(i != string::npos)
	{
		result.push_back(line.substr(i,j-i));
		i = line.find_first_not_of(sep, j);
		j = line.find_first_of(sep,i);
	}
	return result;
}

vector<string> tokenize_full(const string & line, const string & sep)
{
	using namespace std;
	vector<string> result;
	if(line.empty()) return result;
	size_t i = 0, j = line.find_first_of(sep, i);
	while(j != string::npos)
	{
		result.push_back(line.substr(i,j-i));
		i = j+1;
		j = line.find_first_of(sep, i);
	}
	result.push_back(line.substr(i,j-1));
	return result;
}


}
