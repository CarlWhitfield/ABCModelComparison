#include"ABC_model_selection.h"
#include"ABCsettings.h"
#include<boost/random/uniform_01.hpp>
#include<fstream>
#include<file_manip.h>

double ks_statistic_normalised(const std::vector<double> & v0, const std::vector<double> & v1)
{
	size_t n = v0.size();
	size_t m = v1.size();
	std::vector<double> v = v0;
	std::vector<size_t> i(n,0), i1(n,1);
	std::vector<size_t> ind(n+m);
	std::iota(ind.begin(),ind.end(),0);
	v.insert(v.end(),v1.begin(),v1.end());
	i.insert(i.end(),i1.begin(),i1.end());
	std::sort(ind.begin(),ind.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
	std::vector<double> F0(n+m+1,0.0), F1(n+m+1,0.0);
	double max_diff = 0;
	for(size_t k = 0; k < n+m; k++)
	{
		size_t ik = ind[k];
		if(i[ik]==0)
		{
			F0[k+1] = F0[k] + 1.0/n;
			F1[k+1] = F1[k];
		}
		else
		{
			F0[k+1] = F0[k];
			F1[k+1] = F1[k] + 1.0/m;
		}
		double diff = abs(F1[k+1] - F0[k+1]);
		if(diff > max_diff) max_diff = diff;
	}
	return max_diff*sqrt(n*m/(n+m));
}

double DistanceFunctionBase::distance(const std::vector<double> & measured, 
		                              const std::vector<double> & simulated)
{
	//default distance function: root mean squared residual
	double dist = 0;
	for(int i = 0; i < int(measured.size()); i++)
	{
		dist += (measured[i] - simulated[i])*(measured[i] - simulated[i]);
	}
	return sqrt(dist);
}

bool ModelOutputBase::extra_outputs() const
{
	//should return true in dervied object if there are extra outputs to be returned
	return false;
}