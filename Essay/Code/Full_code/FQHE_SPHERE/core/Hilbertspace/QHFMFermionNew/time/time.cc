
#include <iostream>
#include "time.h"
#include <unistd.h>
#include <iomanip>

/*int main()
{
	using hu::timer;
	timer t1;
	t1.start("test");
	sleep(2);
	t1.stop("test");
	t1.showTime("test");
	return 0;
}*/

namespace hu{

timer::name timer::tag_;
vector<timer::name> timer::tagv_;
timer::ump timer::timeManager_;

void timer::
start(const string & s)
{
	tag_.push_back(s);
	maxdp_ = (maxdp_<tag_.size())?tag_.size():maxdp_;
	ump::iterator it = timeManager_.find(join(tag_));
	if( it!=timeManager_.end() )
	{
		//gettimeofday(&(it->second.t1), NULL);
		it->second.t1 = std::chrono::high_resolution_clock::now();
		it->second.clock1 = std::clock();
		it->second.hast1 = true;
	}
	else
	{
		timeInteval ti;
		it = timeManager_.insert({join(tag_), ti}).first;
		tagv_.push_back(tag_);
		//gettimeofday(&(it->second.t1), NULL);
		it->second.t1 = std::chrono::high_resolution_clock::now();
		it->second.clock1 = std::clock();
		it->second.hast1 = true;
	}
}

void timer::
stop(const string & s)
{
	if(tag_.back()!=s)
	{
		cout << "Warning\n";
		cout << "tag_ = " << tag_ << ", but s = " << s << endl;
	}
	else
		timer::stop();
}

void timer::
stop()
{
	//timeval t2;
	//gettimeofday(&t2, NULL);
	auto t2 = std::chrono::high_resolution_clock::now();
	ump::iterator it = timeManager_.find(join(tag_));
	if( it!=timeManager_.end() )
	{
		if(it->second.hast1)
		{
			it->second.hast1 = false;
			//it->second.sec += (t2.tv_sec - it->second.t1.tv_sec) + (double)(t2.tv_usec - it->second.t1.tv_usec)/1000000.0;
			std::chrono::duration<double> diff = t2 - it->second.t1;
			it->second.sec += diff.count();
			it->second.cpusec += (double)(std::clock()-it->second.clock1)/CLOCKS_PER_SEC;
			it->second.num++;
		}
		tag_.pop_back();
	}
}


void timer::
showTime(const name tag, unsigned deepth) const
{
	showTime(tag, deepth, timeManager_.find(join(tag))->second);
}


void timer::
showTime(const name tag, unsigned deepth, timeInteval i) const
{
	std::string ll = "";
	for(int i=deepth;i>2;--i)
		ll += "| ";
	if(deepth!=1) ll+= "|-";
	{
		stringstream temp;
		temp.setf(ios::fixed);
		temp << left << std::setprecision(2) << setw(5) << i.per*100.0 << "%";
		cout << setiosflags(ios::left)
		 	 << setw(l_+10) << ll+tag.back()+"["+temp.str()+"]:";
	}
	stringstream temp;
	temp.setf(ios::fixed);
	temp << setiosflags(ios::left) << right;
	temp << std::setprecision(4) << setw(10) << setfill(' ') << i.sec << "s ";
	temp << setiosflags(ios::right) << setfill(' ');
	temp << "(" 
		 << std::setprecision(2) << setw(5) << double((i.cpusec/i.sec>100)?(99):(i.cpusec/i.sec)) << " "<< setfill(' ')
		 << std::setprecision(2) << setw(5) << i.tper*100.0 << "%) ";
	cout << temp.str();
	cout << "    run " << right << setw(ltimes_) << i.num << " times,   " << left;
	cout << "a: " << left << fixed
		 << std::setprecision(4) << setw(10) << i.sec/i.num << "s" << endl;;
	cout.unsetf(ios::floatfield);
	std::setprecision(std::numeric_limits<double>::digits10 + 1);
}

void timer:: 
insertOthers(unsigned& it, unsigned deepth)
{
	++deepth;
	if(it>=templ_ || tagv_.at(it).size()!=deepth)
		return;
	timeInteval iup;
	name tag;
	double walltime = 0.0;
	if(deepth>1)
	{
		tag = tagv_.at(it-1);
		iup = timeManager_.find(join(tagv_.at(it-1)))->second;
		walltime = iup.sec;
	}
	while( it<templ_ && tagv_.at(it).size()==deepth )
	{
		if(deepth==1)
			walltime += timeManager_.find(join(tagv_.at(it)))->second.sec;
		else if(deepth>1 && tagv_.at(it).back()!="others")
		{
			auto ss = timeManager_.find(join(tagv_.at(it)));
			iup.sec -= ss->second.sec;
			iup.cpusec -= ss->second.cpusec;
			ss->second.per = ss->second.sec/walltime;
		}
		++it;
		insertOthers(it, deepth);
	}
	if(deepth>1 && tag.back()!="others")
	{
		iup.hast1 = false;
		iup.per = iup.sec/walltime;
		tag.push_back("others");
		maxdp_ = (maxdp_<tag.size())?tag.size():maxdp_;
		ump::iterator i = timeManager_.find(join(tag));
		if( i!=timeManager_.end() )
			i->second = iup;
		else
		{
			timeManager_.insert({join(tag), iup});
			tagv_.push_back(tag);
		}
	}
	if(deepth==1)
		for(auto i : tagv_)
		{
			auto ss = timeManager_.find(join(i));
			if(i.size()==1)
				ss->second.per = ss->second.sec/walltime;
			ss->second.tper = ss->second.sec/walltime;
		}
}

void timer::
showAllSorted()
{
	templ_ = tagv_.size();
	std::sort(tagv_.begin(), tagv_.end(), cmptagv(timeManager_));
	unsigned it = 0;
	insertOthers(it,0);
	std::sort(tagv_.begin(), tagv_.end(), cmptagv(timeManager_));
	l_=0;
	for(auto & j : tagv_)
	{
		l_ = (l_<j.back().length())?j.back().length():l_;
		maxdp_ = (maxdp_<j.size())?j.size():maxdp_;
	}
	l_ += maxdp_*2;
	for(auto & j : timeManager_)
	{
		auto tt = unsigned(std::log(j.second.num)/std::log(10.0))+1;
		ltimes_ = (tt>ltimes_)?tt:ltimes_;
	}
	//for(auto & j : tagv_)
	//	showTime(j, j.size());
	for(auto & j : tagv_)
		showTime(j, j.size());
}

bool cmptagv::
operator()(const vector<string> & a, const vector<string> & b)
{
	int l1, l2, l;
	l1 = a.size();
	l2 = b.size();
	l = (l1<l2)?l1:l2;
	string head1{""}, head2{""};
	for(int i=0;i<l;++i)
	{
		head1 += a[i];
		head2 += b[i];
		if(a[i]!=b[i])
			return timeManager_.at(head1).sec>timeManager_.at(head2).sec;
	}
	return l1<l2;
}

string join(const vector<string> tag)
{
	string out = "";
	for(auto const& i : tag)
		out += i;
	return out;
}

} // end hu
