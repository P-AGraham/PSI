
#ifndef __HU_TIMER
#define __HU_TIMER

//#include <sys/time.h>
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <cmath>



namespace hu{

using namespace std;

template<typename T> 
std::ostream & 
operator<<(std::ostream & s, std::vector<T> const& obj)
{
	for(auto &i : obj)
		s << i << " ";
	s << "\n";
	return s;
}


struct timeInteval
{
	double sec = 0.0;
	double cpusec = 0.0;
	double per = 0.0;
	double tper = 0.0;
	//timeval t1;
	std::chrono::time_point<std::chrono::high_resolution_clock> t1;
	std::clock_t clock1;
	size_t num = 0;
	bool hast1 = false;
};

string join(const vector<string> tag);

class timer
{
	using name = vector<string>;
	using ump = unordered_map<string, timeInteval>;
private:
	static name tag_;
	static vector<name> tagv_;
	static ump timeManager_;
	void showTime(const name tag, unsigned deepth=0) const;
	void showTime(const name tag, unsigned deepth, timeInteval i) const;
	void insertOthers(unsigned& it, unsigned deepth=0);
	unsigned l_, maxdp_ = 0;
	unsigned templ_ = 0;
	unsigned ltimes_ = 0;
public:
	timer() { };
	~timer() { };
	void start(const string & s);
	void stop( const string & s);
	void stop();
	void showAllSorted();
};


class cmptagv
{
private:
	const unordered_map<string, timeInteval> & timeManager_;
public:
	cmptagv() = delete;
	cmptagv(const unordered_map<string, timeInteval> & timeManager): timeManager_(timeManager) {};
	~cmptagv() { };
	bool operator()(const vector<string> & a, const vector<string> & b);
};
} //end hu


#endif