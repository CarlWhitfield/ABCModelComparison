#ifndef COMPARTMENTAL_H
#define COMPARTMENTAL_H

#include <vector>
#include <memory>

class FlexibleVolumeElement 
{
protected:
	double Volume, conc_top, conc_bottom;
public:
	FlexibleVolumeElement()
	{
		this->Volume = 0;
		this->conc_top = 0;
		this->conc_bottom = 0;
	}
	FlexibleVolumeElement (const double &V, const double &ct, const double &cb)
	{
		this->Volume = V;
		this->conc_top = ct;
		this->conc_bottom = cb;
	}
	inline double get_volume(){ return (this->Volume); }
	inline double get_ctop(){ return (this->conc_top); }
	inline double get_cbottom(){ return (this->conc_bottom); }
	inline double get_igvol(){ return 0.5*this->Volume*(this->conc_bottom + this->conc_top); }

	inline void add_to_volume(const double & dv){ this->Volume += dv; }
	inline void set_volume(const double & v){ this->Volume = v; }
	inline void update_ctop(const double & c){ this->conc_top = c; }
	inline void update_cbottom(const double & c){ this->conc_bottom = c; }
};

class DSVolume 
{
protected:
	double Volume, conc_top_old, conc_bottom_old;
	void add_element(std::shared_ptr<FlexibleVolumeElement>, const int &);  //function to insert element
	void initialise_from_elements();
public:
	std::vector<std::shared_ptr<FlexibleVolumeElement>> VolumeElements;
	DSVolume(const size_t & reserve_units=20);   //default constructor
	DSVolume(std::shared_ptr<FlexibleVolumeElement>, const size_t & reserve_units=20);
	DSVolume(const std::vector<std::shared_ptr<FlexibleVolumeElement>> &,
		     const size_t & reserve_units=20);   //multiple element constructor

	void shunt_down(const std::vector<std::shared_ptr<FlexibleVolumeElement>> &, 
		                std::vector<std::shared_ptr<FlexibleVolumeElement>> &);  //function to create new element at top
	void shunt_up(const std::vector<std::shared_ptr<FlexibleVolumeElement>> &,
		                std::vector<std::shared_ptr<FlexibleVolumeElement>> &);   //function to create new element at bottom, returns element ejected from top


	inline void update_conc_old()   //stoes old concs from current
	{
		this->conc_top_old = this->get_conc_top();
		this->conc_bottom_old = this->get_conc_bottom();
	}
	inline double get_conc_top() const { return this->VolumeElements[0]->get_ctop(); }
	inline double get_conc_bottom() const { return this->VolumeElements.back()->get_cbottom(); }
	inline double get_volume() const { return this->Volume; }
	inline double sum_ig_volume() const
	{
		double tot_ig_vol = 0;
		for(size_t j = 0; j < this->VolumeElements.size(); j++)
		{
			tot_ig_vol += this->VolumeElements[j]->get_igvol();
		}
		return tot_ig_vol;
	}
};

class LungUnit
{
protected:
	double FRC, volume, IGvolume, vent_ratio, vent_ratio_orig, conc, conc_old;

public:
	LungUnit(){};
	LungUnit(const double & frc, const double & igvol, const double & x);   //lung unit constructor

	virtual void inhale(const std::vector<std::shared_ptr<FlexibleVolumeElement>> & fv, 
		                const double & dt = 0);   //absorb and destroy volume element
	virtual void exhale(const double & dv, 
		     std::vector<std::shared_ptr<FlexibleVolumeElement>> & exhaled, 
			 const double & dt = 0);
	virtual double return_volume_change(const double & mean_dvol, const double & dt);

	inline double get_volume() const { return (this->volume); }
	inline double get_conc() const { return (this->conc); }
	inline double get_FRC_volume() const { return (this->FRC); }
	inline double get_vent_ratio() const { return (this->vent_ratio); }
	inline double get_vent_ratio_original() const { return (this->vent_ratio_orig); }
	inline double get_ig_volume() const { return (this->IGvolume); }

	inline void set_volume(const double & v){ this->volume = v; }
	inline void set_vent_ratio(const double & x){ this->vent_ratio = x; }
	inline void set_orig_vent_ratio(const double & x){ this->vent_ratio_orig = x; }
	inline void set_FRC_volume(const double & f){ this->FRC = f; }
	inline void initialise_conc(const double & c)   //for restart
	{
		this->conc = c;
		this->conc_old = c;
		this->IGvolume = 0.5*this->volume*(this->conc + this->conc_old);
	}

	inline virtual void reset(const double & frc, const double & igvol, const double & x)
	{
		*(this) = LungUnit(frc,igvol,x);  //call copy constructor
	}
};

class MixingPoint
{
private:
	double vol_scale;

public:
	MixingPoint(const double & vs)
	{
		this->vol_scale = vs;
	}
	void mix_elements(const std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> & fv_in_above, 
				      const std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> & fv_in_below, 
				      const std::vector<double> & vol_split_above,
				      const std::vector<double> & vol_split_below,
				      std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> & fv_out_above,
				      std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> & fv_out_below) const;

};

#endif