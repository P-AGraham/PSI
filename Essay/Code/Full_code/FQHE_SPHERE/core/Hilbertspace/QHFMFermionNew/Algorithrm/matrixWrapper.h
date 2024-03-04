

#pragma once

#include<iostream>

class matrixWrapper
{
public:
    double* matrix_ = nullptr;
    long dimRow_ = 0;
    long dimCol_ = 0;
    matrixWrapper() = default;
    matrixWrapper(double* matrix, long m, long n)
    {
        dimRow_ = m;
        dimCol_ = n;
        matrix_ = matrix;
    }
    ~matrixWrapper() = default;
};

std::ostream& operator<<(std::ostream& os, const matrixWrapper& A)
{
    long m = A.dimRow_;
    long n = A.dimCol_;
    const double* a = A.matrix_;
    for(auto i=0;i<m;++i)
    {
        for(auto j=0;j<n;++j)
            os << a[i*n+j] << "\t";
        os << "\n";
    }
    return os;
}