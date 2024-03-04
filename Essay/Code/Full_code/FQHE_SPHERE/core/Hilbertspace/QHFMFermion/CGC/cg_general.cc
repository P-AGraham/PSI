
#include "cg_general.h"

int normCG2Body_j(std::vector<double>& cg, int j1, int j2, int j)
{
    double c, norm = 1.0;
    double s1 = 0.5*(j1-1);
    double s2 = 0.5*(j2-1);
    double s =  0.5*(j -1);

    int m1 = (j+j1-j2-1)/2;  // s+s1-s2
    int m2 = j2-1;           // 2*s2
    cg[m1*j2+m2] = 1.0;

    while( m1<j1-1&&m2>0 )
    {   
        c = -std::sqrt( (s1-(m1-s1) )*(s1+m1-s1+1)/( (s2-(m2-s2)+1)*(s2+m2-s2) ) );  //std::sqrt((j1-m1)*(j1+m1+1)/((j2-m2+1)*(j2+m2)));
        cg[(m1+1)*j2+m2-1] = c*cg[m1*j2+m2];
        norm += (cg[(m1+1)*j2+m2-1]*cg[(m1+1)*j2+m2-1]);
        ++m1;
        --m2;
    }
    norm = sqrt(norm);

    m1 = (j+j1-j2-1)/2;  // s+s1-s2
    m2 = j2-1;           // 2*s2
    cg[m1*j2+m2] /= norm;
    while( m1<j1-1&&m2>0 )
    {
        ++m1;
        --m2;
        cg[m1*j2+m2] /= norm;
    }
    return 0;
}

void calCG2Body_j(std::vector<double>& cg, int j1, int j2, int j)
{
    for(int x=0;x<j1;++x)// m1 = x - s1;
    {
        for(int y=0;y<j2;++y)// m2 = y - s2;
        {
            calcg(cg, j1, j2, j, x, y);
        }
    }
    return;
}

double calcg(std::vector<double>& cg, int j1, int j2, int j, int x, int y)
{
    double s1 = (j1-1)*0.5;
    double s2 = (j2-1)*0.5;
    double s  = (j -1)*0.5;
    double m1 = x - s1;
    double m2 = y - s2;
    if( ((fabs(m1)-s1)>0.1)||((fabs(m2)-s2)>0.1) )
        return 0.0;
    if(cg[x*j2+y]>-1.1)
        return cg[x*j2+y];
    if(fabs(m1+m2)>0.1+s)
    {
        cg[x*j2+y] = 0.0;
        return 0.0;
    }
    double c  = std::sqrt( (s+m1+m2+1)*(s-m1-m2) );
    double c1 = std::sqrt( (s1+m1+1)*(s1-m1) )/c;
    double c2 = std::sqrt( (s2+m2+1)*(s2-m2) )/c;
    cg[x*j2+y] = calcg(cg, j1, j2, j, x+1, y)*c1 + calcg(cg, j1, j2, j, x, y+1)*c2;
    return cg[x*j2+y];
}

/*class cgj1j2j
{
private:
    int j1_; // stored 2*j1+1
    int j2_; // stored 2*j2+1
    int j_; //  stored 2*j +1
    bool physicalQ_ = false;
    std::vector<double> cg_;
public:
    cgj1j2j()==delete;
    cgj1j2j(int j1, int j2, int j);
    ~cgj1j2j();

    const double& cg(int m1, int m2);
}*/
cgj1j2j::cgj1j2j(int j1, int j2, int j)
{
    j1_ = j1;
    j2_ = j2;
    if( j-1<std::abs(j1-j2) || j>j1+j2 ) // |2*j1-2*j2| = 2*j
        printf("Not physical!!(j1=%d, j2=%d, j=%d)\n", j1, j2, j);
    j_ = j;
    physicalQ_ = true;
    cg_ = std::vector<double>( j1*j2, -5.0 );
    normCG2Body_j(cg_, j1_, j2_, j_);
    calCG2Body_j(cg_, j1_, j2_, j_);
    s1_ = 0.5*(j1_-1);
    s2_ = 0.5*(j2_-1);
    s_  = 0.5*(j_ -1);
}

#ifdef __PRINTCG

std::ostream & operator<<(std::ostream & os, const cgj1j2j & obj)
{
    os << std::setiosflags(std::ios::showpos);
    os << "\n************************************************************************************\n";
    os << "ClebschGordan Coefficients:\t" << std::setprecision(2) ;
    os << "j1 = " << obj.s1_ << ", j2 = " << obj.s2_ << ", j = " << obj.s_;
    os << "\n------------------------------------------------------------------------------------\n       | ";
    os << setiosflags(std::ios::fixed) << std::setprecision(10) ;
    for(int y=0;y<obj.j2_;++y)
        os << y-(0.5*obj.j2_-0.5) << "\t";
    os << "|\n------------------------------------------------------------------------------------\n";
    for(int x=0;x<obj.j1_;++x)
    {
        os << std::defaultfloat << setiosflags(std::ios::fixed) << std::setprecision(4) ;
        os << x-(0.5*obj.j1_-0.5) << "| ";
        os << std::defaultfloat << setiosflags(std::ios::scientific) << std::setprecision(6);
        for(int y=0;y<obj.j2_;++y)
            os << obj.cg_[x*obj.j2_+y] << "\t";
        os << "|\n";
    }
    os << "------------------------------------------------------------------------------------\n\n";
    os << std::defaultfloat;
    return os;
}

std::ostream & operator<<(std::ostream & os, const cgj1j2 & obj)
{
    os << std::defaultfloat << setiosflags(std::ios::fixed) << std::setprecision(4);
    os << "(s1="<<obj.s1_<<", s2="<<obj.s2_<<", smin="<<obj.smin_<<", smax"<<obj.smax_<<")\n" << std::defaultfloat;
    for(auto& i : obj.cglist_)
        os << i;
    return os;
}

#endif

cgj1j2::cgj1j2(int j1, int j2)
{
    j1_ = j1;
    j2_ = j2;
    s1_ = 0.5*(j1_-1);
    s2_ = 0.5*(j2_-1);
    smin_ = std::abs(s1_-s2_);
    smax_ = std::abs(s1_+s2_);
    cglist_ = std::vector<cgj1j2j>( round(smax_-smin_)+1 );
    for(double s=smin_;s<smax_+0.1;++s)
    {
        int ind = round(s-smin_);
        int j = round(2*s+1);
        cglist_.at(ind) = std::move( cgj1j2j( j1_, j2_, j) );
    }
}

const cgj1j2j&
cgj1j2::cgs(double s)
{
    if( s-smin_<-0.1 || s-smax_>0.1 )
    {
        printf("Invalid value (s1=+%lf, s2=+%lf, smin=+%lf, smax=+%lf, s=+%lf), exit!\n", s1_, s2_, smin_, smax_, s);
        exit(10);
    }
    return cglist_.at(round(s-smin_));
}

/*int main()
{
    auto cg = cgj1j2j(9, 6, 4);
    std::cout << cg;
    cg = cgj1j2j(6, 9, 13);
    std::cout << cg;
    auto cglist = cgj1j2(6, 9);
    std::cout << cglist;
    std::cout << cglist.cgs(5.5);
    std::cout << cglist.cgj(10);
    std::cout << cglist.cgj(10).cgm(0.5,-1.0);
    std::cout << round(-1.0+4) << std::endl;
    return 0;
}*/

extern "C" double cg_py(double la, double lb, double l, double m1, double m2, double m)
{
    //printf("cg_py: <%lf,%lf;%lf,%lf|%lf,%lf>\n", la, m1, lb, m2, l, m);
    if( la<-1e-5 || lb<-1e-5 || l<-1e-5)
        return 0;
    if( la+lb-l<-1e-5 || l-std::abs(la-lb)<-1e-5 )
        return 0;
    if( std::abs(m1)>la+1e-5 || std::abs(m2)>lb+1e-5 || std::abs(m)>l+1e-5 )
        return 0;
    if( std::abs(m1+m2-m)>1e-5 )
        return 0;
    auto cgc = cgj1j2j( round(2*la+1), round(2*lb+1), round(2*l+1) );
    //printf("\t = %lf\n", cgc.cgm(m1, m2));
    return cgc.cgm(m1, m2);

}

extern "C" double wigner_3j_py(double la, double lb, double l, double m1, double m2, double m)
{
    double comm = ((int((round(-la+lb-m)))&1)?-1.0:1.0);
    comm /= std::sqrt(2*l+1);
    //printf("wigner_3j_py: <%lf,%lf;%lf,%lf|%lf,%lf> = %lf*%lf = %lf", la, m1, lb, m2, l, m, comm, cg_py(la, lb, l, m1, m2, -m), comm * cg_py(la, lb, l, m1, m2, -m));
    return comm * cg_py(la, lb, l, m1, m2, -m);
}