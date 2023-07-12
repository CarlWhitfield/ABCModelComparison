#include"compartmental.h"
#include<algorithm>

DSVolume::DSVolume(const size_t & reserve_units)   //default constructor
{
	this->VolumeElements.reserve(reserve_units);
	this->Volume = 0;
	this->VolumeElements.push_back(std::make_shared<FlexibleVolumeElement>(0, 0, 0));
	this->update_conc_old();
}

DSVolume::DSVolume(std::shared_ptr<FlexibleVolumeElement> el, const size_t & reserve_units )
{
	this->VolumeElements.resize(0);
	this->VolumeElements.reserve(reserve_units);
	this->VolumeElements.push_back(el);
	this->initialise_from_elements();
}

DSVolume::DSVolume(const std::vector<std::shared_ptr<FlexibleVolumeElement>> & VE, 
	               const size_t & reserve_units)    //multiple element constructor
{
	this->VolumeElements = VE;
	if(reserve_units > VE.size())
	{
		this->VolumeElements.reserve(reserve_units);
	}
	this->initialise_from_elements();
}

void DSVolume::initialise_from_elements()
{
	this->Volume = 0;
	for(size_t i=0; i<this->VolumeElements.size(); i++)
	{
		this->Volume += this->VolumeElements[i]->get_volume();
	}
	this->update_conc_old();
}

void DSVolume::add_element(std::shared_ptr<FlexibleVolumeElement> VE, const int &pos)
{
	std::vector<std::shared_ptr<FlexibleVolumeElement>>::iterator it;

	it = this->VolumeElements.begin();
  
	this->VolumeElements.insert(it + pos, VE);
}

void DSVolume::shunt_down(const std::vector<std::shared_ptr<FlexibleVolumeElement>> & NewElements, 
	                      std::vector<std::shared_ptr<FlexibleVolumeElement>> & removed_elements)
{
	this->update_conc_old();

	int pos = 0;
	double VolumeLeft = 0;
	for(size_t i = 0; i < NewElements.size(); i++) //add new elements and count volume
	{
		this->add_element(NewElements[i], pos);
		VolumeLeft += NewElements[i]->get_volume(); 
	}
    while(VolumeLeft > 0)
    {
		if(VolumeElements.back()->get_volume() >= VolumeLeft)  //volume to be removed is smaller than volume of last element
		{
            //new conc at bottom
			double cbnew = this->VolumeElements.back()->get_ctop() 
				         + (this->VolumeElements.back()->get_cbottom() - this->VolumeElements.back()->get_ctop())
						   * (1 - VolumeLeft/this->VolumeElements.back()->get_volume());

			//create new element that is pushed out
			removed_elements.push_back(std::make_shared<FlexibleVolumeElement>(VolumeLeft, 
				                          cbnew, this->VolumeElements.back()->get_cbottom()));
            
			//update old element
            VolumeElements.back()->add_to_volume(-VolumeLeft);
			VolumeElements.back()->update_cbottom(cbnew);
            VolumeLeft = 0;
		}
        else  //volume to be removed is larger than last element
		{
            //push out last element
			removed_elements.push_back(VolumeElements.back());
			//update volume left to be removed
			VolumeLeft -= VolumeElements.back()->get_volume();
			//remove last element
            VolumeElements.pop_back(); 
		}
	}
}

void DSVolume::shunt_up(const std::vector<std::shared_ptr<FlexibleVolumeElement>> & NewElements, 
	                    std::vector<std::shared_ptr<FlexibleVolumeElement>> & removed_elements)
{
	this->update_conc_old();

	double VolumeLeft = 0;
	for(size_t i = 0; i < NewElements.size(); i++)
	{
		VolumeElements.push_back(NewElements[i]);
		VolumeLeft += NewElements[i]->get_volume(); 
	} 
    while(VolumeLeft > 0)
    {
		if(VolumeElements[0]->get_volume() > VolumeLeft)   //volume to be removed is smaller than first element
		{
            //new conc at top
            double ctnew = this->VolumeElements[0]->get_cbottom() 
				           + (this->VolumeElements[0]->get_ctop() - this->VolumeElements[0]->get_cbottom())
						     * (1 - VolumeLeft/this->VolumeElements[0]->get_volume());
           
			//create new element that is pushed out
			removed_elements.push_back(std::make_shared<FlexibleVolumeElement>(VolumeLeft, 
				                                this->VolumeElements[0]->get_ctop(), ctnew));

			//update old top element
            VolumeElements[0]->add_to_volume(-VolumeLeft);
			VolumeElements[0]->update_ctop(ctnew);
            VolumeLeft = 0;
		}
        else    //volume to be removed is larger than first element
		{
			//push out first element
			removed_elements.push_back(VolumeElements[0]);
			//update volume left to be removed
			VolumeLeft -= VolumeElements[0]->get_volume();
			//remove first element
            VolumeElements.erase(VolumeElements.begin());  
		}
	}
}

LungUnit::LungUnit(const double & frc, const double & igvol, const double & x)
{
	this->FRC = frc;
	this->IGvolume = igvol;
	this->volume = frc;
	this->vent_ratio_orig = x;
	this->vent_ratio = x;
	this->conc = igvol / frc;
	this->conc_old = this->conc;
}

double LungUnit::return_volume_change(const double & mean_dvol, const double & dt)
{ 
	if(-this->vent_ratio*mean_dvol > (1 - 1E-06)*this->volume) //about to go into negative volume
	{
		return ((1 - 1E-06)*this->volume);   //prevent this
	}
	else
	{
		return (this->vent_ratio*mean_dvol); 
	}

	
}

void LungUnit::inhale(const std::vector<std::shared_ptr<FlexibleVolumeElement>> & fv, const double & dt)
{
	this->conc_old = this->conc;
	for(size_t n = 0; n < fv.size(); n++)   //add elements in order inhaled
	{
		this->volume += fv[n]->get_volume();
		this->IGvolume += fv[n]->get_igvol();
	}
	this->conc = this->IGvolume / this->volume;   //on inhalation
}

void LungUnit::exhale(const double & dv, std::vector<std::shared_ptr<FlexibleVolumeElement>> & exhaled, 
	                  const double & dt)
{
	this->conc_old = this->conc;
	this->conc = this->IGvolume / this->volume;  //instantaneous mixing
	this->volume -= dv;
	this->IGvolume -= 0.5*dv*(this->conc_old + this->conc);

	exhaled.resize(1);
	exhaled[0] = std::make_shared<FlexibleVolumeElement>(dv, this->conc_old, this->conc);
}

//fv_in vector is vector of vector of elements entering mixing point,
//each vector should be ordered first to last entering
//outputs will be split into private ds volumes, in same order
void MixingPoint::mix_elements(const std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> & fv_in_above, 
				  const std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> & fv_in_below, 
				  const std::vector<double> & vol_split_above,
				  const std::vector<double> & vol_split_below,
				  std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> & fv_out_above,
				  std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> & fv_out_below) const
{
	fv_out_above.resize(vol_split_above.size());
	fv_out_below.resize(vol_split_below.size());
	double totvol = 0;
	for(size_t s = 0; s < vol_split_above.size(); s++)
	{
		totvol += vol_split_above[s];
	}
	for(size_t s = 0; s < vol_split_below.size(); s++)
	{
		totvol += vol_split_below[s];
	}
	size_t Nelements_out = size_t(totvol / this->vol_scale) + 1;   //set scale for mixing
	
	std::vector<double> ct(Nelements_out + 1);
	std::fill(ct.begin(), ct.end(), 0);
	std::vector<double> vol_split(Nelements_out);
	std::fill(vol_split.begin(), vol_split.end(),  this->vol_scale);
	vol_split[Nelements_out-1] = totvol - (Nelements_out-1)* this->vol_scale;
	for(size_t s = 0; s < fv_in_above.size(); s++)
	{
		if(fv_in_above[s].size() > 0)
		{
			double vol = 0;
			for(size_t el = 0; el < fv_in_above[s].size(); el++)
			{
				vol += fv_in_above[s][el]->get_volume();
			}
			size_t el = 0;
			double cumul_vol = fv_in_above[s][el]->get_volume(); //counts volume including el
			ct[0] += vol*fv_in_above[s][el]->get_cbottom(); //first in from above is bottom of element 0
			for(size_t i = 1; i < Nelements_out; i++)
			{
				double vnext = i* this->vol_scale*(vol/totvol);
				while(el < fv_in_above[s].size() && cumul_vol < vnext)  //count until threshold is crossed
				{
					el++;
					cumul_vol += fv_in_above[s][el]->get_volume();  //so solution is in el
				}
				ct[i] += vol*((cumul_vol - vnext)*fv_in_above[s][el]->get_cbottom()
							 + (vnext + fv_in_above[s][el]->get_volume() - cumul_vol)
							   *fv_in_above[s][el]->get_ctop()) / fv_in_above[s][el]->get_volume();  //interpolate
			}
			ct[Nelements_out] += vol*fv_in_above[s].back()->get_ctop(); //first in from above is bottom of element 0
		}
	}
	for(size_t s = 0; s < fv_in_below.size(); s++)
	{
		if(fv_in_below[s].size() > 0)
		{
			double vol = 0;
			for(size_t el = 0; el < fv_in_below[s].size(); el++)
			{
				vol += fv_in_below[s][el]->get_volume();
			}
			size_t el = 0;
			double cumul_vol = fv_in_below[s][el]->get_volume(); //counts volume including el
			ct[0] += vol*fv_in_below[s][el]->get_ctop(); //first in from above is bottom of element 0
			for(size_t i = 1; i < Nelements_out; i++)
			{
				double vnext = i* this->vol_scale*(vol/totvol);
				while(el < fv_in_below[s].size() && cumul_vol < vnext)  //count until threshold is crossed
				{
					el++;
					cumul_vol += fv_in_below[s][el]->get_volume();  //so solution is in el
				}
				ct[i] += vol*((cumul_vol - vnext)*fv_in_below[s][el]->get_ctop()
							 + (vnext + fv_in_below[s][el]->get_volume() - cumul_vol)
							   *fv_in_below[s][el]->get_cbottom()) / fv_in_below[s][el]->get_volume();  //interpolate
			}
			ct[Nelements_out] += vol*fv_in_below[s].back()->get_cbottom(); //first in from above is bottom of element 0
		}
	}

	if(totvol > 0)
	{
		for(size_t i = 0; i <= Nelements_out; i++)
		{
			ct[i] /= totvol;
		}
	}
	for(size_t n = 0; n < vol_split_above.size(); n++)
	{
		fv_out_above[n].resize(Nelements_out);
		for(size_t i = 0; i < Nelements_out; i++)
		{
			if(totvol > 0)
			{
				fv_out_above[n][i] = std::make_shared<FlexibleVolumeElement>(
				           vol_split_above[n]*vol_split[i]/totvol, ct[i], ct[i+1]);
			}
			else
			{
				fv_out_above[n][i] = std::make_shared<FlexibleVolumeElement>(
				                                             0, ct[i], ct[i+1]);
			}
		}
	}
	for(size_t n = 0; n < vol_split_below.size(); n++)
	{
		fv_out_below[n].resize(Nelements_out);
		for(size_t i = 0; i < Nelements_out; i++)
		{
			if(totvol > 0)
			{
				fv_out_below[n][i] = std::make_shared<FlexibleVolumeElement>(
				           vol_split_below[n]*vol_split[i]/totvol, ct[i+1], ct[i]);
			}
			else
			{
				fv_out_below[n][i] = std::make_shared<FlexibleVolumeElement>(
				                                          0, ct[i+1], ct[i]);
			}
		}
	}
}
