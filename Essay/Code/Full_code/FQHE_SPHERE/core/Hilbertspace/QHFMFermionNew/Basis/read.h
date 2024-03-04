
#ifndef __READ_WRITE__H
#define __READ_WRITE__H


#include <iostream> 
#include <fstream> 
#include <vector> 
#include <string> 

template<class T>
std::vector<T> read_vec(std::string &fname)
{
    std::ifstream s(fname.c_str(),std::ios::binary); 
    if(!s.good()) 
        throw std::runtime_error("Couldn't open file \"" + fname + "\" for reading");
    T data;
    size_t length=0;
    while(s.read( (char*)&data, sizeof(data) ))
    {
        ++length;
        //std::cout << data << "\t";
    }
    s.clear();
    s.seekg(0, std::ios::beg);

    std::vector<T> result(length);
    s.read( (char*)result.data(), sizeof(data)*length );
    s.close(); 

    //for(auto i:result)
    //    std::cout << "\n" << i << "\t";
    //std::cout << std::endl;
    return result;
}
 
void
writeStringToFile(const std::string& fname, const std::string& t) 
{ 
    std::ofstream s(fname.c_str()); 
    if(!s.good()) 
        throw std::runtime_error("Couldn't open file \"" + fname + "\" for writing");
    s << t << std::endl;
    s.close();
}

#endif