
#pragma once

#include <iostream>
#include <vector>

template<class T, class T_INS, class T_END>
void print_vec(std::ostream& os, const std::vector<T>& obj, T_INS insert, T_END end)
{
    for(const auto& i: obj)
        std::cout << i << insert;
    std::cout << end;
}

template<class T, class T_INS, class T_END>
void print_vec(const std::vector<T>& obj, T_INS insert, T_END end)
{
    print_vec(std::cout, obj, insert, end);
}

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& obj)
{
    print_vec(os, obj, "\n", "\n");
    return os;
}