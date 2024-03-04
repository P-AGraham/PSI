
#include "QHFM.h"
#include "Basis/Basis.h"
#include "Basis/Bipartition.h"
#include "Basis/Basis.cc"
#include "Basis/Bipartition.cc"

#include "Algorithrm/gemm.h"
//--------------------------------------------------------------------------------------------------------------------


int main(int argc, char* argv[])
{	
	hu::timer hutm;
	std::string labels = "N_o\tN_e\ttheta_p\tRES:S(A)\t\terror\t\t\tRES:S(A+B)\t\terror\t\t\tmutual\t\t\tuncertainty\n";
	std::string res = labels;
	if(argc < 2) 
	{ 
	   //reminds us to give an input file if we forget
	   printf("Usage: %s Ne \n",argv[0]); 
	   return 0; 
	}
	size_t N_o = std::stoi(argv[1]);
	size_t N_phi = N_o-1;
	size_t N_e = N_o;
	int L_z = (N_o-1)*N_e/2;
	std::string readpath = "/lab_mem/huliangdong/share/QHFM_ED/";
	std::string outputdir = "/lab_mem/huliangdong/share/QHFM_ED/mutual/";

    std::string command;
    command = "mkdir -p " + outputdir;
    int status = system(command.c_str());
    assert(status==0);

	// read wavefunction 
	std::string path = readpath + "/N" + std::to_string(N_o) + "/vec_0.bin";
	std::cout << path << std::endl;
	std::vector<Data_dtype> vec = read_vec<Data_dtype>(path);
	printf("%lu\n", vec.size());

	// generate Hilbert space
	auto basis = lzBasis(N_o, N_e, L_z);
	basis.generateIndex();

	for(double theta_p=0.02;theta_p<0.49999;theta_p+=0.02)
	{
		hutm.start("mutual_ED");
		auto [SA, errorA, SAB, errorAB] = mutual_ED(theta_p, vec, basis, "dsyevd");
		hutm.stop("mutual_ED");
		double mutual = 2.0*SA-SAB;
		printf("Mutual information: 2*SA-SAB = %.16lf\n", 2*SA-SAB);
		double uA = (errorA<1e-20)?0.0:(-errorA*std::log(errorA));
		double uAB = (errorAB<1e-20)?0.0:(-errorAB*std::log(errorAB));
		// (truncerrAB<1e-20)?0.0:(-truncerrAB*std::log(truncerrAB/numeigsAB));
		std::string temp = boost::str( boost::format("%02d\t%02d\t%.4lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n") % N_o % N_e % theta_p % SA % errorA % SAB % errorAB % (2.0*SA-SAB) % (2*uA+uAB));
		res += temp;
		writeStringToFile(outputdir+"/"+std::to_string(N_o)+".dat", res);
		std::cout << labels+temp;
		hutm.showAllSorted();
	}
	std::cout << res;
	hutm.showAllSorted();
	std::cout << "All is done!\n";
	return 0;
}