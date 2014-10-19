/*
 * Hamiltonian.cpp
 * PNT
 *
 * Created on: Sep 30, 2008
 *
 * Copyright (C) 2008 Yi Pan
 */

#include "Hamiltonian.h"

int hamilton_count = 0;
int dh_count = 0;
time_t hamilton_time = 0;

Hamiltonian::Hamiltonian(int n_variables, int n_saves) 
		   : Flags_ptr(), Wave_params_ptr(), Tuning_params_ptr(), 
		     variables(n_variables), saves(n_saves)
{
}

Hamiltonian::~Hamiltonian()
{
}
/**
 * Store the flags in a Hamiltonian object
 */
void Hamiltonian::set_flags_ptr(copy_ptr<Flags> new_flags_ptr /**< Pointer to a Flags object */)
{
	Flags_ptr = new_flags_ptr;
	return;
}

/**
 * Store the waveform parameters in a Hamiltonian object
 */
void Hamiltonian::set_wave_params_ptr(copy_ptr<WaveformParameters> new_params_ptr /**< Pointer to a WaveformParameters object */)
{
	Wave_params_ptr = new_params_ptr;
	return;
}

/**
 * Store the adjustable parameters in a Hamiltonian object
 */
void Hamiltonian::set_tuning_params_ptr(copy_ptr<TuningParameters> new_params_ptr /**< Pointer to a TuningParameters object */)
{
	Tuning_params_ptr = new_params_ptr;
	return;
}

/**
 * Store the configuration variables in a Hamiltonian object
 */
void Hamiltonian::set_variables(std::vector<double> new_variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */)
{
	if (new_variables.size() != variables.size())
	{
		std::cout << "Error! Hamiltonian::set_variables: ";
		std::cout << "wrong number of variables." << std::endl;
		abort();
	}
	variables = new_variables;
	return;
}

/**
 * Store constant variables in a Hamiltonian object to speed up computation
 */
void Hamiltonian::set_saves(std::vector<double> new_saves /**< Array of the constant quantities to be set */)
{
	if (new_saves.size() != saves.size())
	{
		std::cout << "Error! Hamiltonian::set_saves: ";
		std::cout << "wrong number of saves." << std::endl;
		abort();
	}
	saves = new_saves;
	return;
}

/**
 * Compute the A potential for the stored dynamical configuration
 */
double Hamiltonian::evaluate_current_eoba() const
{
	return evaluate_eoba(get_variables());
}

/**
 * Compute the D potential for the stored dynamical configuration
 */
double Hamiltonian::evaluate_current_eobd() const
{
	return evaluate_eobd(get_variables());
}

/**
 * Compute Heff for the stored dynamical configuration
 */
double Hamiltonian::evaluate_current_Heff() const
{
	return evaluate_Heff(get_variables());
}

/**
 * Compute dA/dr for the stored dynamical configuration
 */
double Hamiltonian::evaluate_current_eoba_deriv() const
{
	return evaluate_eoba_deriv(get_variables());
}

/**
 * Compute dD/dr for the stored dynamical configuration
 */
double Hamiltonian::evaluate_current_eobd_deriv() const
{
	return evaluate_eobd_deriv(get_variables());
}

/**
 * Compute the combinations of the physical spins on which H depends, using the stored
 * dynamical configuration
 */
std::vector<double> Hamiltonian::evaluate_current_Seffs() const
{
	return evaluate_Seffs(get_variables());
}

/**
 * Compute Hreal for the stored dynamical configuration
 */
double Hamiltonian::evaluate_current() const
{
	return evaluate(get_variables());
}

/**
 * Compute the A potential for the dynamical configuration provided in the argument
 */
double Hamiltonian::evaluate_eoba(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::evaluate_eoba: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute the D potential for the dynamical configuration provided in the argument
 */
double Hamiltonian::evaluate_eobd(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::evaluate_eobd: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute Heff potential for the dynamical configuration provided in the argument
 */
double Hamiltonian::evaluate_Heff(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::evaluate_Heff: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute the dA/dr potential for the dynamical configuration provided in the argument
 */
double Hamiltonian::evaluate_eoba_deriv(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::evaluate_eoba_deriv: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute the dD/dr potential for the dynamical configuration provided in the argument
 */
double Hamiltonian::evaluate_eobd_deriv(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::evaluate_eobd_deriv: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute the combinations of the physical spins on which H depends, using the 
 * dynamical configuration provided in the argument
 */
std::vector<double> Hamiltonian::evaluate_Seffs(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::evaluate_Seffs: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	std::vector<double> dummyvec;
	return dummyvec;
}

/**
 * Compute Hreal for the dynamical configuration provided in the argument
 */
double Hamiltonian::evaluate(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::evaluate: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute dHreal/dpr for the dynamical configuration provided in the argument
 */
double Hamiltonian::dHdpr(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::dHdpr: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute dHreal/dpphi for the dynamical configuration provided in the argument
 */
double Hamiltonian::dHdpphi(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::dHdpphi: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute dHreal/dr for the dynamical configuration provided in the argument
 */
double Hamiltonian::dHdr(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	std::cout << "Error: Hamiltonian::dHdr: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute the partial derivative of Hreal with respect to the specified coordinate, for the stored dynamical configuration
 */
double Hamiltonian::evaluate_deriv(unsigned int variable_index /**< Choose the variable for differentiation, it must be an int b/w 0 and the size of the array of dynamical variables */) const
{
	extern int dh_count;
	dh_count += 1;
	std::vector<double> variables = get_variables();
	copy_ptr<Flags>		flags_ptr = get_flags_ptr();
	
	if (variable_index >= variables.size() || variable_index < 0)
	{
		std::cout << "Error! Hamiltonian::evaluate_deriv: ";
		std::cout << "wrong variable index." << std::endl;
		abort();
	}

	if (flags_ptr -> get_flag_Hderiv() == 1) /**< Richardson-extrapolated derivative*/
	{
		double h = 0.1 * variables[variable_index];
		double err = 1.;
                //h = (h == 0 ? 1 : h);
                h = (fabs(h) < 1.0e-9 ? 1.0e-8 : h);
		double Rder = dfridr_hamiltonian(*this, variable_index, variables, h, &err);
		//return dfridr_hamiltonian(*this, variable_index, variables, h, &err);
		return Rder;
	}
	else if (flags_ptr -> get_flag_Hderiv() == 0) /**< Finite-difference derivative*/
	{
		std::vector<double> newvariables;
		double hp, hm;
		//		double hpp, hmm;
		double variable_sav = variables[variable_index];
		double h = 0.001 * variable_sav;
        // TODO Not sure which of two choices is better
//		h = (fabs(h) < 1.0e-9 ? 1 : h);
                //h = (h == 0 ? 1 : h);
                //h = (fabs(h) < 1.0e-9 ? 0.99e-8 : h);
		if ( fabs(h) < 1.0e-9){
                    //std::cout << "Stas h = " << h << std::endl; 
		    h = 0.99e-8;
		    //h = 1.e-6;
		}
		newvariables = variables;
//		for (unsigned int i = 0; i < newvariables.size(); i++)
//			std::cout << newvariables[i] << "   ";
//		std::cout << std::endl;

		newvariables[variable_index] = variable_sav + h;
//		for(int i=0; i<12; i++)
//			std::cout << newvariables[i] << "   ";
//		std::cout << std::endl;
		hp = evaluate(newvariables);

//		std::cout << "new var+: " << newvariables[variable_index] << std::endl;

		newvariables[variable_index] = variable_sav - h;
//		for(int i=0; i<12; i++)
//			std::cout << newvariables[i] << "   ";
//		std::cout << std::endl;
		hm = evaluate(newvariables);

//		std::cout << "new var-: " << newvariables[variable_index] << std::endl;
//		newvariables[variable_index] = variable_sav + 2. * h;
//		hpp = evaluate(newvariables);
//		newvariables[variable_index] = variable_sav - 2. * h;
//		hmm = evaluate(newvariables);
//		std::cout << "h: " << h << std::endl;
//		std::cout << "hp and hm: " << hp << "   " << hm << std::endl;
//		std::cout << "return index " << variable_index << " : " << (hp - hm) / (2. * h) << std::endl;
		return (hp - hm) / (2. * h);
//		return (8. * hp - 8. * hm - hpp + hmm) / (12. * h);
	}
	else
	{
		std::cout << "Error! Hamiltonian::evaluate_deriv: ";
		std::cout << "flag_Hderiv is not 0 or 1." << std::endl;
		std::cout << flags_ptr -> get_flag_Hderiv() << std::endl;
		std::cout << this->evaluate_radius() << std::endl;
		abort();
	}
}

/**
 * Compute the non-Keplerian r for a Cartesian Hamiltonian for the stored dynamical configuration
 */
double Hamiltonian::nonKeplerianRadius() const
{
	std::cout << "Error: Hamiltonian::nonKeplerianRadius: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute the non-Keplerian r for a spherical Hamiltonian for the stored dynamical configuration
 */
double Hamiltonian::nonKeplerianRadiusSph() const
{
	copy_ptr<Flags>			flags_ptr	  = get_flags_ptr();
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	//	unsigned int i;
	double pp, w, rnk;
	int warn_flag=0;
	int nppdiv = 1000, npp = 1, ppsuccess = 0;
	double pp1 = 0.01, pp2 = 40.0, pp1sol, pp2sol;
	std::vector<double> pp1list, pp2list;
	//We increase the sampling at small radii just to be sure we don't miss the root
	if (this -> evaluate_radius()>3.)
	{
		nppdiv=10;
	}
	else
	{
		nppdiv=1000;
	}
	zbrak_hamiltonian(*this, pp1, pp2, nppdiv, pp1list, pp2list, &npp);	
	ppsuccess = npp;
	if (npp > 1)
	{
		std::cout << "Warning! Hamiltonian::nonKeplerianRadius: ";
		std::cout << npp << " roots bracketed, using the last one.";
		std::cout << " Wrong root might be used at ";
		std::cout << "r = " << variables[0] << std::endl;
	}
	if (ppsuccess != 0)
	{
		pp1sol = pp1list[npp-1];
		pp2sol = pp2list[npp-1];
	}
	else
	{
		ppsuccess = zbrac_hamiltonian(*this, &pp1, &pp2);
		
		if (ppsuccess != 0)
		{
			pp1sol = pp1;
			pp2sol = pp2;
		}
		else
		{
//			std::cout << "Error! Hamiltonian::nonKeplerianRadius: ";
//			std::cout << "pp root could not be bracketed at radius ";
//			std::cout << variables[0] << ",";
//			std::cout << "(pp, dH/dr) information printed below: ";
//			for (int i = 0; i < 100; i++)
//			{
//				std::cout << i * 0.01 << "   ";
//				std::cout << this -> dHdr_pp(i * 0.01) << std::endl;
//			}
		//beware: modified
		warn_flag=0;
		}
	}
//	std::cout << "pp1sol " << pp1sol << std::endl;
//	std::cout << "pp2sol " << pp2sol << std::endl;
//	for (int i = 0; i < 100; i++)
//	{
//		std::cout << i * 0.05 << "   ";	
//		std::cout << this -> dHdr_pp(i * 0.05) << std::endl;
//	}	
	

	


	if (warn_flag==0)
	{
//		std::cout << "pp " << pp << std::endl;
		variables[3] = 0.0;
		if (Flags_ptr->get_flag_rNK_pphi() == 0) 
		{
			pp = zbrent_hamiltonian(*this, pp1sol, pp2sol, 0.00000000001);
			variables[5] = pp;
		}
		temphptr->set_variables(variables);
//		 unsigned int i;
//		std::cout << "Vars used for pphi"<< std::endl;
//		for (int i = 0; i < 12; i++)
//		 	std::cout << variables[i] << "   "<< std::endl;
		w = temphptr->evaluate_deriv(5); /**< Compute dH/dpphi for pr=0 */
//		std::cout << "w "<< w << std::endl;
		
//		std::cout << "Hsph "<< temphptr->evaluate_current()<< std::endl;

////////////////////////////////////////////////////////////////////////
//		Use mixed substitution of pphi (Mar 25, 2011)
//		copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
//		double eta;
//		eta	 = wave_params_ptr->get_eta();
//		double r = this -> evaluate_radius();
//		pp = variables[5];
//		std::vector<double> tmpCartvars(12,0.);
//		tmpCartvars[0] = r;
//		double Apot = this -> evaluate_eoba(tmpCartvars);
//		double Aprime = this -> evaluate_eoba_deriv(tmpCartvars);
//		double wpot = Apot*(1. + pp*pp/r/r);
//		std::cout << "w pporb       = " << w << std::endl;
//
//		Correct mixed substitution
//		w = 1./pow(r*pow(2.*(1. + 2*eta*(sqrt(wpot)-1.))/(r*r*Aprime),1.0/3.0),3.0/2.0);
//		std::cout << "w mixed ver 1 = " << w << std::endl;
//
//		Full expression w/out any substitution w/ pcirc
//		double ppcirc = zbrent_hamiltonian(*this, pp1sol, pp2sol, 0.00000000001);
//		w = 1./pow(r*pow((1.+ppcirc*ppcirc/r/r)*r/ppcirc/ppcirc/Apot*(1. + 2*eta*(sqrt(wpot)-1.)),1.0/3.0),3.0/2.0);
//		std::cout << "w mixed ver 2 = " << w << std::endl;
//		std::cout <<  " " << std::endl;
//
//		FILE *outH = fopen( "outputNKrad.txt", "a" );
//		if( outH != NULL )
//			fprintf(outH,"%3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n",r,pp,w,wpot,Apot,Aprime);
//		fclose(outH);	
///////////////////////////////////////////////////////////////////////		

//		std::cout << "circular orbit energy: " << temphptr->evaluate_current() << std::endl;
		if (variables[0] < 1.)
			std::cout << "circular orbit pp: " << pp << std::endl;
//		 std::cout << "circular orbit w: " << w << std::endl;
		rnk=1.0 / pow(w, 2.0/3.0);
//		std::cout << "r_NK = " << rnk << std::endl;
//		abort();
	}
	else
	{
		rnk = this -> evaluate_radius();
		std::cout << "r_NK = r = " << rnk << std::endl;
	}
	return rnk;
}


//double Hamiltonian::nonKeplerianRadiusSph() const
//{
//	copy_ptr<Hamiltonian> temphptr(this->clone());
//	double rnk = this -> evaluate_radius();
//	return rnk;
//}

/**
 * Compute the light-ring radius for a Cartesian Hamiltonian for the stored dynamical configuration
 */
double Hamiltonian::LightRingRadius() const
{
	std::cout << "Error: Hamiltonian::LightRingRadius: ";
	std::cout << "member function not implemented ";
	std::cout << "for this derived class of Hamiltonian." << std::endl;
	abort();
	return 0.;
}

/**
 * Compute the light-ring radius for a spherical Hamiltonian for the stored dynamical configuration
 */
double Hamiltonian::LightRingSph() const
{
	std::vector<double> variables = get_variables();
	std::vector<double> tmpvariables = get_variables();
	Array rvec(3);
	copy_ptr<Hamiltonian> temphptr(this->clone());
	copy_ptr<Flags>		flags_ptr       = get_flags_ptr();
	//Choose the number of points you want to sample in the search range: nrdiv
	int nrdiv = 200, nr = 1, rsuccess = 0;
	//Choose the range where to search for the LR: [r1,r2]
	double r1 = 0.5, r2 = 10., r1sol, r2sol;
	std::vector<double> r1list, r2list;	
	double rsol;
	
//	std::cout << "Looking for LR" << std::endl;
	flags_ptr->set_flag_LR(1);	
	//Actually, variables should already be spherical
	rvec[0] =  variables[0];
	rvec[1] =  variables[1];
	rvec[2] =  variables[2];
	
	//This is not really important, since we solve for r
	variables[0] = sqrt(rvec[0] * rvec[0]
						+ rvec[1] * rvec[1]
						+ rvec[2] * rvec[2]);
	variables[1] = Pi/2.;
	variables[2] = 0.;
	//Remember: for photons the result is independent of pphi
	variables[3] = 0.;
	variables[4] = 0.;
	variables[5] = 1.;

	temphptr->set_variables(variables);
	temphptr->set_flags_ptr(flags_ptr);
    
	//Note that we don't have multiple roots now: we just pick
	//the first root we find starting from r2 and moving towards r1
	//in order to avoid nan for small values of r while evaluating dH/dr
	zbrak_LR(*temphptr, r1, r2, nrdiv, r1list, r2list, &nr);
//	std::cout << "nr " << nr << std::endl;
	rsuccess = nr;
	if (nr > 1)
	{
		std::cout << "Warning! Hamiltonian::LightRingRadius: ";
		std::cout << nr << " roots bracketed, using the last one.";
		std::cout << " Wrong root might be used at ";
		std::cout << "r = " << variables[0] << std::endl;
	}
	if (rsuccess != 0)
	{
		r1sol = r1list[nr-1];
		r2sol = r2list[nr-1];
	}
	else
	{
		rsuccess = zbrac_LR(*temphptr, &r1, &r2);
		
		if (rsuccess != 0)
		{
			r1sol = r1;
			r2sol = r2;
		}
		else
		{
			std::cout << "Error! Hamiltonian::LightRingRadius: ";
			std::cout << "r root could not be bracketed at radius ";
			std::cout << variables[0] << ",";
			std::cout << "(r, dH/dr) information printed below: ";
			for (int i = 0; i < 100; i++)
			{
				std::cout << i * 0.1 << "   ";
				std::cout << temphptr -> dHdr_r(i * 0.1) << std::endl;
			}					
			abort();
		}
	}
		
//	    double delta=(5.-1.5)/100;
//		for (int i = 1; i < 100; i++)
//		{
//			variables[0] =1.5+delta*i;
//			temphptr->set_variables(variables);
//			std::cout << variables[0] << " ";
//			std::cout << temphptr->evaluate_deriv(0) << std::endl;
//		}
			
//	std::cout << "r1sol = " << r1sol << std::endl;
//	std::cout << "r2sol = " << r2sol << std::endl;
	rsol = zbrent_LR(*temphptr, r1sol, r2sol, 0.000000001);
//	std::cout << "LR = " << rsol << std::endl;

	flags_ptr->set_flag_LR(0);
	temphptr->set_variables(tmpvariables);
	temphptr->set_flags_ptr(flags_ptr);
	return rsol;
}


//**********************************************************************//
//
// EOBnospinHamiltonian member functions:
//
//**********************************************************************//
/**
 * Cartesian Hamiltonian of Phys. Rev. D 84, 124052 (2011) 
 */
EOBnospinHamiltonian::EOBnospinHamiltonian() : Hamiltonian(4, 14)
{
}

EOBnospinHamiltonian::~EOBnospinHamiltonian()
{
}

void EOBnospinHamiltonian::generate_saves()
{
	std::vector<double> saves(14);
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	double eta, lam, z1, z2;
	double eta2, a5, a4;
	eta	 = wave_params_ptr->get_eta();
	lam	 = tuning_params_ptr->get_a5();
	z1	 = tuning_params_ptr->get_z1();
	z2	 = tuning_params_ptr->get_z2();
	saves[0] = (6 * (4 - 3 * eta) * eta - 8 * z1 - 4 * z2) / 3; //z3
	saves[1] = eta2 = eta * eta;				    //eta2
	saves[2] = a5 = lam * eta;				    //a5
	saves[3] = a4 = (94./3. - 41./32.*Pi2) * eta - z1;	    //a4
	saves[4] = 32 - 24*eta -4*a4 - a5;			    //eoba_n3
	saves[5] = -16 + a4 + 8*eta;				    //eoba_n4
	saves[6] = -a4*a4 - 8*a5 - 8*a4*eta + 2*a5*eta - 16*eta2;   //eoba_d0
	saves[7] = -8*a4 - 4*a5 - 2*a4*eta - 16*eta2;		    //eoba_d1
	saves[8] = -4*a4 - 2*a5 - 16*eta;			    //eoba_d2
	saves[9] = -2*a4 - a5 - 8*eta;				    //eoba_d3
	saves[10]= -16 + a4 + 8*eta;				    //eoba_d4
	saves[11]= 36*eta2;					    //eobd_d0
	saves[12]= 2*eta*(26-3*eta);				    //eobd_d1
	saves[13]= 6*eta;					    //eobd_d2
	set_saves(saves);
	return;
}

/**
 * Compute r for the stored dynamical configuration
 */
double EOBnospinHamiltonian::evaluate_radius() const
{
	std::vector<double> variables = get_variables();
	return variables[0];
}

/**
 * Compute orbital freq for the stored dynamical configuration
 */
double EOBnospinHamiltonian::evaluate_orb_omega() const
{
	return this->evaluate_deriv(3);
}


double EOBnospinHamiltonian::evaluate_eoba(std::vector<double> variables) const
{
	copy_ptr<Flags>		flags_ptr       = get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 	saves 	   	= get_saves();
	double eoba, eta, r, r2, r3, r4, r5;
	
	eta 	= wave_params_ptr->get_eta();
	r	= variables[0];
	r2	= r * r;
	r3	= r * r2;
	r4	= r2 * r2;
	r5	= r2 * r3;
	
	switch (flags_ptr->get_flag_eoba())
	{
		case 0: 
			eoba = 1 - 2 / r + 2 * eta / r3 + saves[3] / r4 + saves[2] / r5;
			break;
		case 1:
			eoba = (r3 * saves[4] + r4 * saves[5]) /
			   (saves[6] + r * saves[7] + r2 * saves[8]
					  + r3 * saves[9] + r4 * saves[10]);
			break;
		default: 
			std::cout << "Error! EOBnospinHamiltonian::evaluate_eoba: ";
			std::cout << "flagPadeAr must be 0 or 1" << std::endl;
			abort();
	}
	return eoba;
}

double EOBnospinHamiltonian::evaluate_eoba_deriv(std::vector<double> variables) const
{
	copy_ptr<Flags>		flags_ptr       = get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double>	saves 	   	= get_saves();
	double dAdr, dAdr_denomenator, eta, r, r2, r3, r4, r5, r6;
	
	eta 	= wave_params_ptr->get_eta();
	r	= variables[0];
	r2	= r * r;
	r3	= r * r2;
	r4	= r2 * r2;
	r5	= r2 * r3;
	r6	= r3 * r3;
	
	switch (flags_ptr->get_flag_eoba())
	{
		case 0: 
			dAdr = 2 / r2 - 6 * eta / r4 - 4 * saves[3] / r5 - 5 * saves[2] / r6;
			break;
		case 1:
			dAdr_denomenator = (saves[6] + r * saves[7] + r2 * saves[8]
					  		  + r3 * saves[9] + r4 * saves[10]);
			dAdr = r2 * (3 * saves[4] * saves[6] + 
					r * (4 * saves[5] * saves[6] + 2 * saves[4] * saves[7] + 
					r * (3 * saves[5] * saves[7] + saves[4] * saves[8] + 
					r * (2 * saves[5] * saves[8] + 
					r * (saves[5] * saves[9] - saves[4] * saves[10]))))) 
					/ dAdr_denomenator / dAdr_denomenator;
			   ;
			break;
		default: 
			std::cout << "Error! EOBnospinHamiltonian::evaluate_eoba: ";
			std::cout << "flagPadeAr must be 0 or 1" << std::endl;
			abort();
	}
	return dAdr;
}

double EOBnospinHamiltonian::evaluate_eobd(std::vector<double> variables) const
{
	copy_ptr<Flags>		     flags_ptr     = get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double>	     saves 	   = get_saves();
	double eobd, eta, r, r2, r3, r4;
	
	eta 	= wave_params_ptr->get_eta();
	r	= variables[0];
	r2	= r * r;
	r3	= r * r2;
	r4	= r2 * r2;
	
	switch (flags_ptr->get_flag_eobd())
	{
		case 0:
			eobd = 1 - 6 * eta / r2 + (-52 * eta + 6 * saves[1]) / r3;
			break;
		case 1:
			eobd = r3 / (r3 + r * saves[13] + saves[12]);
			break;
		case 2:
			eobd = r4 / (r4 + r2 * saves[13] + r * saves[12] + saves[11]);
			break;
		default: 
			std::cout << "Error! EOBnospinHamiltonian::evaluate_eobd: ";
			std::cout << "flagPadeDr must be 0, 1 or 2" << std::endl;
			abort();
	}
	return eobd;
}

double EOBnospinHamiltonian::evaluate_eobd_deriv(std::vector<double> variables) const
{
	copy_ptr<Flags>			flags_ptr	= get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double>		saves 	   	= get_saves();
	double dDdr, r, r2, r3, r4, denom;

//	eta 	= wave_params_ptr->get_eta();
	r	= variables[0];
	r2	= r * r;
	r3	= r * r2;
	r4      = r2 * r2;
	denom = (saves[11] + saves[12]*r + saves[13]*r2 + r4);

	switch (flags_ptr->get_flag_eobd())
	{
		case 0:
			dDdr = 0.;
			std::cout << "Error! EOBnospinHamiltonian::evaluate_eobd_deriv: ";
			std::cout << "dDdr not implemented for this flag_eobd" << std::endl;
			abort();
			break;
		case 1:
			dDdr = r2*(3.*saves[12] + 2.*r*saves[13])
				/ (r3 + r*saves[13] + saves[12]) / (r3 + r*saves[13] + saves[12]);
			break;
		case 2:
			dDdr = r3*(4.*saves[11] + r*( 3.*saves[12] + r*( 2.*saves[13] ) ) )
				/ denom / denom;
			break;
		default: 
			std::cout << "Error! EOBnospinHamiltonian::evaluate_eobd: ";
			std::cout << "flagPadeDr must be 0, 1 or 2" << std::endl;
			abort();
	}
	return dDdr;
}

double EOBnospinHamiltonian::evaluate_Heff(std::vector<double> variables) const
{
	copy_ptr<Flags>                 flags_ptr     = get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 		saves 	      = get_saves();
	double z1, z2, z3, r, pr, pp, r2, r4, pr2, pr4, pp2, pp4, pp2pr2;
	double eoba, eobd, Heff;
	eoba = evaluate_eoba(variables);
	eobd = evaluate_eobd(variables);
	
	z1	= tuning_params_ptr->get_z1();
	z2	= tuning_params_ptr->get_z2();
	z3	= saves[0];
	r	= variables[0];
	pr	= variables[2];
	pp	= variables[3];
	r2	= r * r;
	r4	= r2 * r2;
	pr2 	= pr * pr;
	pr4 	= pr2 * pr2;
	pp2	= pp * pp;
	pp4	= pp2 * pp2;
	pp2pr2	= pp2 * pr2;
	
	switch (flags_ptr->get_flag_tortoiseR())
	{
		case 0:
			Heff = sqrt( eoba * ( 1 + pp2/r2 + (eoba / eobd) * pr2 
						+ 1/r2 * (z1 * (pp4/r4 + 2*pp2pr2/r2 + pr4)
						+ z2 * (pp2pr2/r2 + pr4) + z3 * pr4) ) );
			break;
		case 1:
			Heff = sqrt(pr2 + eoba * (1 + pp2/r2
					+ 1/r2 * (z1 * (pp4/r4 + 2*pp2pr2/r2 + pr4)
					+ z2 * (pp2pr2/r2 + pr4) + z3 * pr4)));
			break;
		
		default: 
			std::cout << "Error! EOBnospinHamiltonian::evaluate_Heff: ";
			std::cout << "flag_tortoiseR must be 0 or 1" << std::endl;
			abort();
	}
		
	return Heff;
}

double EOBnospinHamiltonian::evaluate(std::vector<double> variables) const
{
	extern int hamilton_count;
	extern time_t hamilton_time;
	time_t track_time;
	hamilton_count += 1;
	track_time = time(NULL);

	double eta, Heff, Hreal;
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();

	eta	 = wave_params_ptr->get_eta();
	Heff = evaluate_Heff(variables);

	hamilton_time += time(NULL) - track_time;

	if (eta < 1.0e-16)
		Hreal = Heff;
	else
		Hreal = sqrt(1 + 2 * eta * (Heff-1)) / eta;
	return Hreal;
}

double EOBnospinHamiltonian::dHdpr(std::vector<double> variables) const
{
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	double Hreal, Heff, r, pr, z3, onebyHHeta;
	double eta, eobA, eobD, dHdpr, pr3, r2;
	std::vector<double> saves = get_saves();
	z3 = saves[0];
	Heff = evaluate_Heff(variables);
	Hreal = evaluate(variables);
	eobA = evaluate_eoba(variables);
	eobD = evaluate_eobd(variables);
	eta = wave_params_ptr->get_eta();
	r = variables[0];
	pr = variables[2];
	r2 = r * r;
	pr3 = pr * pr * pr;
	onebyHHeta = 1./ ( eta * Hreal * Heff );

	if (eta < 1.0e-16)
		dHdpr = 1 / Heff * eobA * ( 2.*z3 * pr3 / r2 + pr * eobA / eobD );
	else 
		dHdpr = onebyHHeta * eobA * ( 2.*z3 * pr3 / r2 + pr * eobA / eobD );

	return dHdpr;
}

double EOBnospinHamiltonian::dHdpphi(std::vector<double> variables) const
{
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	double Hreal, Heff, r, pphi, onebyHHeta;
	double eta, eobA, dHdpphi, r2;
	Heff = evaluate_Heff(variables);
	Hreal = evaluate(variables);
	eobA = evaluate_eoba(variables);
	eta = wave_params_ptr->get_eta();
	r = variables[0];
	pphi = variables[3];
	r2 = r * r;
	onebyHHeta = 1./ ( eta * Hreal * Heff );
	
	if (eta < 1.0e-16)
		dHdpphi = 1 / Heff * pphi * eobA / r2;
	else
		dHdpphi = onebyHHeta * pphi * eobA / r2;

	return dHdpphi;
}

double EOBnospinHamiltonian::dHdr(std::vector<double> variables) const
{
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	double Hreal, Heff, eta, r, pr, pp, z3, onebyHHeta;
	double eobA, eobD, dAdr, dDdr, dHdr, pr2, pr4, pp2, r2, r3;
	std::vector<double> saves = get_saves();
	z3 = saves[0];
	Heff = evaluate_Heff(variables);
	Hreal = evaluate(variables);
	eobA = evaluate_eoba(variables);
	eobD = evaluate_eobd(variables);
	dAdr = evaluate_eoba_deriv(variables);
	dDdr = evaluate_eobd_deriv(variables);
	eta = wave_params_ptr->get_eta();
	r = variables[0];
	pr = variables[2];
	pp = variables[3];
	r2 = r * r;
	r3 = r * r2;
	pp2 = pp * pp;
	pr2 = pr * pr;
	pr4 = pr2 * pr2;
	onebyHHeta = 1./ ( eta * Hreal * Heff );
	
	if (eta < 1.0e-16)
		dHdr = 1 / Heff * ( ( 1. + (pp2 + z3*pr4) / r2 + pr2 * eobA / eobD ) * dAdr
			+ eobA * ( - 2.*pp2 / r3 - 2.*z3*pr4 / r3 + pr2 * dAdr / eobD
			- pr2 * eobA * dDdr / eobD / eobD ) ) / 2.;
	else
		dHdr = onebyHHeta * ( ( 1. + (pp2 + z3*pr4) / r2 + pr2 * eobA / eobD ) * dAdr
			+ eobA * ( - 2.*pp2 / r3 - 2.*z3*pr4 / r3 + pr2 * dAdr / eobD
			- pr2 * eobA * dDdr / eobD / eobD ) ) / 2.;

	return dHdr;
}

// This version replace circular orbit pp value everywhere in the Hamiltonian
// It is different from the simpler expression used in all current models.
/* double EOBnospinHamiltonian::nonKeplerianRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	double pp, w;
	pp = zbrent_hamiltonian(*this, 0.0, 5.0, 0.000001);
	variables[2] = 0.0;
	variables[3] = pp;
	temphptr->set_variables(variables);
	w = temphptr->evaluate_deriv(3);
	return 1.0 / pow(w, 2.0/3.0);
} */

double EOBnospinHamiltonian::nonKeplerianRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double>   variables 	 = get_variables();
	double eta, r, pp, eoba, dAdr, psi;
	eta	 = wave_params_ptr->get_eta();
	r	 = variables[0];
	pp	 = variables[3];
	eoba = evaluate_current_eoba();
	dAdr = evaluate_current_eoba_deriv();
	
	psi = 2. * (1. + 2. * eta * (sqrt(eoba * (1 + pp * pp / r / r)) - 1)) 
			 / (r * r * dAdr);
	return r * pow(psi, 1.0/3.0);
}

double EOBnospinHamiltonian::dHdr_r(double r) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	variables[0] = r;
	variables[2] = 0.0;
	temphptr->set_variables(variables);
	return temphptr->evaluate_deriv(0);
}

double EOBnospinHamiltonian::dHdr_pp(double pp) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	variables[2] = 0.0;
	variables[3] = pp;
	temphptr->set_variables(variables);
	return temphptr->evaluate_deriv(0);
}

double EOBnospinHamiltonian::LightRingRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	temphptr -> LightRingRadius();
	return 0.0;
}


//**********************************************************************//
//
// EOBspinHamiltonian member functions:
//
//**********************************************************************//
/**
 * Cartesian Hamiltonian of Phys. Rev. D 81, 084041 (2010)
 */
EOBspinHamiltonian::EOBspinHamiltonian() : Hamiltonian(12, 40)
{
}

EOBspinHamiltonian::~EOBspinHamiltonian()
{
}


void EOBspinHamiltonian::generate_saves()
{
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> saves(40);
	double eta, z1;
	double eta2, a6, a5, a4;

	eta	= wave_params_ptr->get_eta();
	a5	= eta * (tuning_params_ptr->get_a5());
	a6	= eta * (tuning_params_ptr->get_a6());
	z1 	= tuning_params_ptr->get_z1();
	saves[0] = eta2 = eta * eta;					//eta2
	saves[1] = a5;							//a5
	saves[2] = a4 = (94./3. - 41./32.*Pi2) * eta - z1;		//a4
	saves[3] = sqrt(1. - 4. * eta);					//dm
	saves[4] = (1. + saves[3]) / (1. - saves[3]);			//m1om2
	saves[5] = 1. / saves[4];					//m2om1
	saves[6] = 32 - 24*eta -4*a4 - a5;				//eoba_n3
	saves[7] = -16 + a4 + 8*eta;					//eoba_n4
	saves[8] = -a4*a4 - 8*a5 - 8*a4*eta + 2*a5*eta - 16*eta2;	//eoba_d0
	saves[9] = -8*a4 - 4*a5 - 2*a4*eta - 16*eta2;			//eoba_d1
	saves[10]= -4*a4 - 2*a5 - 16*eta;				//eoba_d2
	saves[11]= -2*a4 - a5 - 8*eta;					//eoba_d3
	saves[12]= -16 + a4 + 8*eta;					//eoba_d4
	saves[13]= 36*eta2;						//eobd_d0
	saves[14]= 2*eta*(26-3*eta);					//eobd_d1
	saves[15]= 6*eta;						//eobd_d2
	// below are coefficients needed by Pade(1,5) of spin A(r)
	saves[16] = -64. + 12.*a4 + 64.*eta + 4.*a5 + a6 - 4.*eta2;
	saves[17] = 32. - 4.*a4 - 24.*eta - a5;
	saves[18] = 4.*a4*a4 + 32.*a4*eta + 4.*a4*a5 + 16.*a6 - a4*a6
	+ 4.*a4*eta2 + 16.*a5*eta + a5*a5 - 8.*a6*eta + 32.*eta*eta2;
	saves[19] = 4.*a4*a4 + 32.*a4*eta + 16.*a5 + a4*a5 + 8.*a6 
	+ 32.*eta2 - 2.*a6*eta + 8.*eta*eta2;
	saves[20] = 16.*a4 + 8.*a4*eta + 8.*a5 + 4.*a6 + 32.*eta2 + 2.*a5*eta;
	saves[21] = 8.*a4 + 32.*eta + 4.*a5 + 2.*a6 - 8.*eta2;
	saves[22] = 4.*a4 + 16.*eta + 2.*a5 + a6 - 4.*eta2;
	saves[23] = 32. - 4.*a4 - 24.*eta - a5;
	
	saves[24] = 80. - 2.*a4 - 24.*eta;
	saves[25] = -32. + 4.*eta;
	saves[26] = 2.*a4*a4 + 8.*a4*eta + 16.*a5 - 12.*a6 + 48.*eta2 - 8.*a5*eta;
	saves[27] = 12.*a4 - 8.*a5 + a6 + 12.*eta2;
	saves[28] = -3.*a4 + 16.*eta;
	saves[29] = 16.*a4 - 4.*a5 - 4.*a6 + 32.*eta2;
	saves[30] = -2.*a4 + 24.*eta - a5;
	saves[31] = -4.*a4 + 32.*eta - a6 - 4.*eta2;
	saves[32] = 2.*a4 - 4.*eta;
	saves[33] = -8.*eta + a5;
	saves[34] = 8. - 4.*eta;
	saves[35] = 16. - 2.*a4 - 16.*eta;
	saves[36] = -32. + 4.*eta;
	// below are coefficients needed by Pade(0,6) of spin A(r)
	saves[37] = 64. - 12.*a4 - 64.*eta - 4.*a5*eta - a6*eta + 4.*eta2;
	saves[38] = 32. - 4.*a4 - 24.*eta - a5*eta;
	saves[39] = 16. - a4 - 8.*eta;
	set_saves(saves);
	return;
}

/**
 * Compute r for the stored dynamical configuration
 */
double EOBspinHamiltonian::evaluate_radius() const
{
	std::vector<double> variables = get_variables();	
	return sqrt(variables[0]*variables[0]+variables[1]*variables[1]+variables[2]*variables[2]);
}

/**
 * Compute orbital freq for the stored dynamical configuration
 */
double EOBspinHamiltonian::evaluate_orb_omega() const
{
	std::vector<double> variables = get_variables();	
	Array rvec(3), rdvec(3), omegavec(3);
	
	rvec[0] = variables[0];
	rvec[1] = variables[1];
	rvec[2] = variables[2];
	
	for (int i = 0; i < 3; i++)
	{
		rdvec[i] = this->evaluate_deriv(i + 3);		
	}
	
	omegavec = crossproduct(rvec, rdvec) / innerproduct(rvec, rvec);
	return sqrt(innerproduct(omegavec, omegavec));
}

double EOBspinHamiltonian::evaluate_eoba(std::vector<double> variables) const
{
	copy_ptr<Flags>		flags_ptr	= get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 	saves 	   	= get_saves();
	int i;
	double eoba, eta, eta2, e5, r, r2, r3, r4, r5, r6, aSK2, aSK4, aSK6, aSK8, aSK10;
	Array rvec(3), pvec(3), S1vec(3), S2vec(3), SKerr_vec(3);
	std::vector<double> Seffs(6);
	
	eta  = wave_params_ptr->get_eta();
	eta2 = eta*eta;
	e5   = tuning_params_ptr->get_e5();
	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
	r  = sqrt(innerproduct(rvec, rvec));
	r2 = r * r;
	r3 = r2 * r;
	r4 = r3 * r;
	r5 = r4 * r;
	r6 = r5 * r;
	
	Seffs = evaluate_Seffs(variables);
	for (i = 0; i < 3; i++)
	{
		SKerr_vec[i] = Seffs[i];
	}
//	if (flags_ptr->get_flag_SEKmodel() == 0)
//		SKerr_vec = S1vec * (1+saves[5]) + S2vec * (1+saves[4]);
//	else
//		SKerr_vec = S1vec * (1+0.75*saves[5]) + S2vec * (1+0.75*saves[4]);
	aSK2	  = innerproduct(SKerr_vec, SKerr_vec);
	aSK4	  = aSK2 * aSK2;
	aSK6	  = aSK4 * aSK2;
	aSK8	  = aSK4 * aSK4;
	aSK10	  = aSK2 * aSK8;
	/* case 2 is for the strict 3PN eoba, i.e. Pade(1,3)[A(r)] */
	switch (flags_ptr->get_flag_eoba())
	{
		case 1:
			eoba = (r3 * (saves[6]-e5*eta*aSK2 - 32*aSK2 + 6*aSK4 + 4*eta*aSK2) 
				 +  r4 * (saves[7] + 12*aSK2 - aSK4)) 
				 / (	 (saves[8]+(2*eta-8)*e5*eta*aSK2 - 8*saves[2]*aSK2 
					+ 2*saves[2]*aSK4 - aSK8 + 4*saves[1]*aSK2
					+ 4*aSK2*e5*eta*aSK2 - 12*eta*aSK4 - 8*eta2*aSK2)
		   		   + r * (saves[9]-(4-aSK2)*e5*eta*aSK2 + (saves[1]-16*eta)*aSK2 
					- 2*eta*aSK4- 2*aSK6)
		   		   + r2* (saves[10]-2*e5*eta*aSK2 - saves[2]*aSK2 - 4*aSK4 + aSK6)
		   		   + r3* (saves[11]-e5*eta*aSK2 - 8*aSK2 + 4*aSK4 + 4*eta*aSK2)
		   		   + r4* (saves[12] + 12*aSK2 - aSK4));
		   	break;
		case 10:
			eoba = (r3 * saves[6] +  r4 * saves[7]) 
				 / ( saves[8] + r * saves[9] + r2 * saves[10] + r3 * saves[11] + r4 * saves[12])
				 + aSK2 / r2;
			break;
		case 2:
		   	eoba = (r2 * (16 - saves[2] - 12 * aSK2 + aSK4 - 8*eta) 
		   		 +  r3 * (-8 + 4 * aSK2 + 2 * eta)) 
		   		 / (	 (-4*saves[2] - 4*eta2 + saves[2]*aSK2 - 8*eta*aSK2 - aSK6) 
		   		   + r * (-saves[2] - 4*eta - eta*aSK2 - aSK4) * 2. 
		   		   + r2* (-saves[2] - 4*eta - 4*aSK2 + aSK4) 
		   		   + r3* (-8 + 4 * aSK2 + 2 * eta));	
		   	break;
		case 20:
		   	eoba = (r2 * (16 - saves[2] - 8*eta) +  r3 * (-8 + 2 * eta)) 
		   		 / (	(-4*saves[2] - 4*eta2) + r * (-saves[2] - 4*eta) * 2. 
					+ r2* (-saves[2] - 4*eta) + r3* (-8 + 2 * eta)) + aSK2 / r2;	
		   	break;
		case 15:
			eoba = 	( r4 * (saves[16] + saves[24]*aSK2 - 24.*aSK4 + aSK6)
				+ r5 * (saves[17] + saves[25]*aSK2 + 6.*aSK4))
				/
				( (saves[18] + saves[26]*aSK2 + saves[27]*aSK4 + saves[28]*aSK6 + aSK10)
				+ r * (saves[19] + saves[29]*aSK2 + saves[30]*aSK4 + 2.*eta*aSK6 + 2.*aSK8)
				+ r2* (saves[20] + saves[31]*aSK2 + saves[32]*aSK4 + 4.*aSK6 - aSK8)
				+ r3* (saves[21] + saves[33]*aSK2 + saves[34]*aSK4 - 4.*aSK6)
				+ r4* (saves[22] + saves[35]*aSK2 - 12.*aSK4 + aSK6)
				+ r5* (saves[23] + saves[36]*aSK2 + 6.*aSK4));
			break;
		case 150:
			eoba = 	( r4 * saves[16] + r5 * saves[17]) /
				( saves[18] + r * saves[19] + r2* saves[20] + r3* saves[21] 
				+ r4* saves[22] + r5* saves[23]) + aSK2 / r2;
			break;
		case 6:
			eoba = r6 / ( r6 + 2.*r5 + (4.-aSK2)*r4 - (-8.+4.*aSK2+2.*eta)*r3 
				+ (saves[39]-12.*aSK2+aSK4) * r2 + (saves[38]-32.*aSK2+6.*aSK4+4.*aSK2*eta) * r 
				+ saves[37]-80.*aSK2+24.*aSK4-aSK6+2.*aSK2*saves[2]+24.*aSK2*eta);
			break;
		case 60:
			eoba = r6 / ( r6 + 2.*r5 + 4.*r4 - (-8.+2.*eta)*r3 
				+ saves[39] * r2 + saves[38] * r + saves[37]) + aSK2/r2;
			break;
		default: 
			std::cout << "Error! EOBspinHamiltonian::evaluate_eoba: ";
			std::cout << "flagPadeAr must be 1, 10, 2, 15 or 150" << std::endl;
	}
	return eoba;
}

double EOBspinHamiltonian::evaluate_eoba_deriv(std::vector<double> variables) const
{
	copy_ptr<Flags>		flags_ptr	= get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 	saves 	   	= get_saves();
	int i;
	double dAdr, eta, eta2, e5, r, r2, r3, r4, aSK2, aSK4, aSK6, aSK8;
	Array rvec(3), pvec(3), S1vec(3), S2vec(3), SKerr_vec(3);
	std::vector<double> Seffs(6);
	
	eta  = wave_params_ptr->get_eta();
	eta2 = eta*eta;
	e5   = tuning_params_ptr->get_e5();
	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
	r  = sqrt(innerproduct(rvec, rvec));
	r2 = r * r;
	r3 = r2 * r;
	r4 = r3 * r;
	
	Seffs = evaluate_Seffs(variables);
	for (i = 0; i < 3; i++)
	{
		SKerr_vec[i] = Seffs[i];
	}
//	if (flags_ptr->get_flag_SEKmodel() == 0)
//		SKerr_vec = S1vec * (1+saves[5]) + S2vec * (1+saves[4]);
//	else
//		SKerr_vec = S1vec * (1+0.75*saves[5]) + S2vec * (1+0.75*saves[4]);
	aSK2	  = innerproduct(SKerr_vec, SKerr_vec);
	aSK4	  = aSK2 * aSK2;
	aSK6	  = aSK4 * aSK2;
	aSK8	  = aSK4 * aSK4;
		
	double n3, n4;
	double d0, d1, d2, d3, d4;
	/* case 2 is for the strict 3PN eoba, i.e. Pade(1,3)[A(r)] */
	switch (flags_ptr->get_flag_eoba())
	{
		case 1:
			n3 = (saves[6]-e5*eta*aSK2 - 32*aSK2 + 6*aSK4 + 4*eta*aSK2);
			n4 = (saves[7] + 12*aSK2 - aSK4);
			d0 = (saves[8]+(2*eta-8)*e5*eta*aSK2 - 8*saves[2]*aSK2 
				+ 2*saves[2]*aSK4 - aSK8 + 4*saves[1]*aSK2
				+ 4*aSK2*e5*eta*aSK2 - 12*eta*aSK4 - 8*eta2*aSK2);
			d1 = (saves[9]-(4-aSK2)*e5*eta*aSK2 + (saves[1]-16*eta)*aSK2 
					- 2*eta*aSK4- 2*aSK6);
			d2 = (saves[10]-2*e5*eta*aSK2 - saves[2]*aSK2 - 4*aSK4 + aSK6);
			d3 = (saves[11]-e5*eta*aSK2 - 8*aSK2 + 4*aSK4 + 4*eta*aSK2);
			d4 = (saves[12] + 12*aSK2 - aSK4);
			dAdr 	= (3.*n3*r2 + 4.*n4*r3) / (d0 + d1*r + d2*r2 + d3*r3 + d4*r4)
				- (n3*r3 + n4*r4) * (d1 + 2.*d2*r + 3.*d3*r2 + 4.*d4*r3)
				/ (d0 + d1*r + d2*r2 + d3*r3 + d4*r4)
				/ (d0 + d1*r + d2*r2 + d3*r3 + d4*r4);
		   	break;
		default: 
			std::cout << "Error! EOBspinHamiltonian::evaluate_eoba_deriv: ";
			std::cout << "flagPadeAr must be 1, only" << std::endl;
	}
	return dAdr;
}

double EOBspinHamiltonian::evaluate_eobd(std::vector<double> variables) const
{
	copy_ptr<Flags>			flags_ptr	  = get_flags_ptr();
	copy_ptr<WaveformParameters> 	wave_params_ptr   = get_wave_params_ptr();
	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 		saves 	   	  = get_saves();
	int i;
	double eobd, eta, r, r2, r3, r4;
	Array rvec(3), pvec(3), S1vec(3), S2vec(3);
	
	eta	= wave_params_ptr->get_eta();
	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
	r  = sqrt(innerproduct(rvec, rvec));
	r2	= r * r;
	r3	= r * r2;
	r4	= r2 * r2;
	
	switch (flags_ptr->get_flag_eobd())
	{
		case 0:
			eobd = 1 - 6 * eta / r2 + (-52 * eta + 6 * saves[1]) / r3;
			break;
		case 1:
			eobd = r3 / (r3 + r * saves[15] + saves[14]);
			break;
		case 2:
			eobd = r4 / (r4 + r2 * saves[15] + r * saves[14] + saves[13]);
			break;
		default: 
			std::cout << "Error! EOBspinHamiltonian::evaluate_eobd: ";
			std::cout << "flagPadeDr must be 0, 1 or 2" << std::endl;
	}
	// std::cout << flags_ptr->get_flag_eobd() << "   " << eobd;
	// std::cout << "   " << r3 << "   " << saves[14] << "   " << saves[15] << std::endl;
	return eobd;
}

double EOBspinHamiltonian::evaluate_eobd_deriv(std::vector<double> variables) const
{
	copy_ptr<Flags>			flags_ptr	  = get_flags_ptr();
	copy_ptr<WaveformParameters> 	wave_params_ptr   = get_wave_params_ptr();
	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 		saves 	   	  = get_saves();
	int i;
	double dDdr, eta, r, r2, r3, r4;
	Array rvec(3), pvec(3), S1vec(3), S2vec(3);
	
	eta	= wave_params_ptr->get_eta();
	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
	r  = sqrt(innerproduct(rvec, rvec));
	r2 = r * r;
	r3 = r * r2;
	r4 = r2 * r2;
	
	switch (flags_ptr->get_flag_eobd())
	{
		case 0:
			dDdr 	= 12. * eta / r3 - 3. * (-52 * eta + 6 * saves[1]) / r4;
			break;
		case 1:
			dDdr 	= 3. * r2 / (r3 + r * saves[15] + saves[14])
				- r3 * (3. * r2 + saves[15]) / (r3 + r * saves[15] + saves[14]) 
				/ (r3 + r * saves[15] + saves[14]);
			break;
		case 2:
			dDdr 	= 4. * r3 / (r4 + r2 * saves[15] + r * saves[14] + saves[13])
				- r4 * (4. * r3 + 2. * saves[15] * r + saves[14])
				/ (r4 + r2 * saves[15] + r * saves[14] + saves[13])
				/ (r4 + r2 * saves[15] + r * saves[14] + saves[13]);
			break;
		default: 
			std::cout << "Error! EOBspinHamiltonian::evaluate_eobd: ";
			std::cout << "flagPadeDr must be 0, 1 or 2" << std::endl;
	}
	// std::cout << flags_ptr->get_flag_eobd() << "   " << eobd;
	// std::cout << "   " << r3 << "   " << saves[14] << "   " << saves[15] << std::endl;
	return dDdr;
}

std::vector<double> EOBspinHamiltonian::evaluate_Seffs(std::vector<double> variables) const
{
	copy_ptr<Flags>					flags_ptr	  	= get_flags_ptr();
	copy_ptr<WaveformParameters> 	wave_params_ptr = get_wave_params_ptr();
	std::vector<double> 			saves 	  		= get_saves();
	
	unsigned int i;
	std::vector<double> Seffs(6);	
	Array rvec(3), pvec(3), S1vec(3), S2vec(3), SKerr_vec(3), Sstar_vec(3);
	double m1om2 = saves[4];
	double m2om1 = saves[5];
	double eta = wave_params_ptr -> get_eta();

	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
	
	if (eta < 1.0e-16)
	{
		SKerr_vec = S1vec;
		Sstar_vec[0] = 0.;
		Sstar_vec[1] = 0.;
		Sstar_vec[2] = 0.;
	}
	else
	{
		if (flags_ptr->get_flag_SEKmodel() == 0)
			SKerr_vec = S1vec * (1+m2om1) + S2vec * (1+m1om2);
		else
			SKerr_vec = S1vec * (1+0.75*m2om1) + S2vec * (1+0.75*m1om2);
		Sstar_vec = S1vec * m2om1 + S2vec * m1om2;
	}

	for (i = 0; i < 3; i++)
	{
		Seffs[i]	= SKerr_vec[i];
		Seffs[i + 3]	= Sstar_vec[i];
	}
	return Seffs;
}

double EOBspinHamiltonian::evaluate_Heff(std::vector<double> variables) const
{
	copy_ptr<Flags>			flags_ptr	  = get_flags_ptr();
	copy_ptr<WaveformParameters> 	wave_params_ptr   = get_wave_params_ptr();
	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 		saves 	   	  = get_saves();
	
	int i;
	double r, r2, r3, r4, p2;
	double ndotp, ndotp2, SKhxrdotp, betadotp;
	std::vector<double> Seffs(6);
	Array rvec(3), pvec(3), S1vec(3), S2vec(3), transpvec(3);
	Array nvec(3), S_vec(3), S0_vec(3), S0_hat(3), Sstar_vec(3);
	Array Seff_vec(3), Seff_hat(3), L_vec(3), L_hat(3);
	Array sigma_vec(3), SKerr_vec(3), SKerr_hat(3), SKhxr(3);
	double gS, gS_star;
	double aSK, kSK, aSK2, kSK2, rhoSK2;
	double eoba, eobd, Deltat, Deltar, alpha_d, alpha, alpha2, SstarSkerr;
	double xi=1., xi2=1.;
	Array beta(3);
	Matrix gamma(3, 3, 0);
	double pvec_gamma_pvec = 0;
	double DHeff, Heff;
	//double M 		= wave_params_ptr->get_M();
	double eta		= wave_params_ptr->get_eta();
	//double c1		= wave_params_ptr->get_chi1();
	//double c2		= wave_params_ptr->get_chi2();
	//double lam		= tuning_params_ptr->get_a5();
	double a_eta		= tuning_params_ptr->get_aeta();
	double b_eta		= tuning_params_ptr->get_beta();
	double a_p4		= tuning_params_ptr->get_a_ph4();
	double a_r2		= tuning_params_ptr->get_a_r2();
	double b_p4		= tuning_params_ptr->get_b_ph4();
	double b_r2		= tuning_params_ptr->get_b_r2();
	double dheff_SS		= tuning_params_ptr->get_dheff_SS();
//	eta2	= saves[0];
//	a5	= saves[1];
//	a4	= saves[2];
//	dm	= saves[3];
	
//	std::cout << "beta = " << b_eta << std::endl;
//	std::cout << "dheffSS = " << dheff_SS << std::endl;
//	std::cout << "w1 =" << tuning_params_ptr->get_wfd1() << std::endl;
//	std::cout << "w2 =" << tuning_params_ptr->get_wfd2() << std::endl;
	
	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
//		std::cout << rvec[i] << "     " << pvec[i] << "     ";
//		std::cout << S1vec[i] << "     " << S2vec[i] << std::endl;
	}

	r  = sqrt(innerproduct(rvec, rvec));
	r2 = r * r;
	r3 = r2 * r;
	r4 = r3 * r;

	nvec = rvec / r;
	p2 = innerproduct(pvec, pvec);
//	ph2 = p2 / eta2;
//	ndotph = innerproduct(nvec, pvec) / eta;
//	ndotph2 = ndotph * ndotph;
	ndotp = innerproduct(nvec, pvec);
	ndotp2 = ndotp * ndotp;
//	L_vec  = crossproduct(rvec, pvec);
//	L_hat  = L_vec / sqrt(innerproduct(L_vec, L_vec));

	S_vec = S1vec + S2vec;
	Seffs = evaluate_Seffs(variables);
	for (i = 0; i < 3; i++)
	{
		SKerr_vec[i] = Seffs[i];
		Sstar_vec[i] = (1.- flags_ptr -> get_flag_LR())*Seffs[i + 3];
	}

	aSK2	   = innerproduct(SKerr_vec, SKerr_vec);
	aSK	       = sqrt(aSK2);
	if (aSK < 1e-16)
		SKerr_hat = nvec;
	else
		SKerr_hat = SKerr_vec / aSK;
//	SKerr_hat  = SKerr_vec / aSK;
	kSK	   	   = innerproduct(nvec, SKerr_hat);
	kSK2	   = kSK * kSK;
	rhoSK2	   = r2 + aSK2 * kSK2;
	SKhxr	   = crossproduct(SKerr_hat, rvec);
	SstarSkerr = innerproduct(Sstar_vec, SKerr_vec);
	
	eoba = evaluate_eoba(variables);
	eobd = evaluate_eobd(variables);

	Deltat	= r2 * eoba;
	Deltar	= Deltat / eobd;
	alpha_d = (r2 + aSK2) * (r2 + aSK2) - aSK2 * Deltat * (1-kSK2);
	alpha	= sqrt(rhoSK2 * Deltat / alpha_d);
	alpha2	= alpha * alpha;
	beta	= aSK * (r2 + aSK2 - Deltat) * SKhxr / alpha_d;
	
	if (flags_ptr->get_flag_tortoiseR() == 0)
	{
		transpvec = pvec;
		SKhxrdotp  = innerproduct(SKhxr, transpvec);

		betadotp= innerproduct(beta, transpvec);	
		pvec_gamma_pvec =(Deltar * ndotp2 + r2 * (p2 - ndotp2)
				- aSK2 / Deltat * SKhxrdotp * SKhxrdotp)
				/ rhoSK2 + betadotp * betadotp / alpha2;
		gS	= 2.0 + (0.375*eta + a_eta) * p2
			- (4.5*eta + 3.0*a_eta) * ndotp2
			- (eta + a_eta)/r + a_p4 * p2 * p2 + a_r2/r2;
		gS_star = 1.5 + (-0.625 + 0.5*eta + b_eta) * p2
			- (3.75*eta + 3.0*b_eta) * ndotp2
			- (0.5 + 1.25*eta + b_eta)/r + b_p4 * p2 * p2 + b_r2/r2;
	}
	else
	{
		//xi = alpha*sqrt(Deltar)/sqrt(rhoSK2);
		xi = Deltat / sqrt(eobd) / (r2 + aSK2);	
		xi2 = xi * xi;	
		for (i = 0; i < 3; i++)
			transpvec[i] = pvec[i] - rvec[i]/r*(xi-1.)/xi*ndotp;
		betadotp= innerproduct(beta, transpvec);	
		SKhxrdotp  = innerproduct(SKhxr, transpvec);

		pvec_gamma_pvec =(Deltar * ndotp2 / xi2 + r2 * (p2 - ndotp2)
				- aSK2 / Deltat * SKhxrdotp * SKhxrdotp)
				/ rhoSK2 + betadotp * betadotp / alpha2;
		gS	= 2.0 + (0.375*eta + a_eta) * (p2 - ndotp2)
			+ ((0.375*eta + a_eta) - (4.5*eta + 3.0*a_eta)) 
			* ndotp2 / xi2 
			- (eta + a_eta)/r + a_p4 * p2 * p2 + a_r2/r2;
		gS_star = 1.5 + (-0.625 + 0.5*eta + b_eta) * (p2 - ndotp2)
			+ ((-0.625 + 0.5*eta + b_eta) - (3.75*eta + 3.0*b_eta)) 
			* ndotp2 / xi2 
			- (0.5 + 1.25*eta + b_eta)/r + b_p4 * p2 * p2 + b_r2/r2;
	}
	sigma_vec = ((gS-2.0) * S_vec + (gS_star-2.0) * Sstar_vec) / 2.0;
	
	switch (flags_ptr->get_flag_DHeffmodel())
	{
		case 0:
			DHeff = 0.;
			break;
		case 1:
			DHeff = 2./r3 * innerproduct(transpvec, crossproduct(sigma_vec, rvec)) 
			      + dheff_SS * SstarSkerr / r4;
			break;
		case 2:
			DHeff = (1 - Deltat/r2) / r2
			      * innerproduct(transpvec, crossproduct(sigma_vec, rvec)) 
			      + dheff_SS * SstarSkerr / r4;
			break;
		case 100:
			DHeff = (r2 + aSK2 - Deltat) / alpha_d
			      * innerproduct(transpvec, crossproduct(sigma_vec, rvec))
			      + dheff_SS * SstarSkerr / r4;
			break;
	}
	
	DHeff =  (1. - flags_ptr -> get_flag_LR())*DHeff;

	Heff	= betadotp
		+ alpha * sqrt(1. - flags_ptr->get_flag_LR() + pvec_gamma_pvec
		+ 2 * (4 - 3 * eta) * eta / rhoSK2 * ndotp2 * ndotp2 / xi2 / xi2 )
		+ DHeff;
	
//	std::cout << "eoba : " << eoba << std::endl;
//	std::cout << "eobd : " << eobd << std::endl;
//	std::cout << "Deltar * ndotp2 : " << Deltar * ndotp2 << std::endl;
//	std::cout << "r2 * (ph2 - ndotph2) : " << r2 * (p2 - ndotp2) << std::endl;
//	std::cout << r2 << "   " << pvec << "   " << p2 << "   " << ndotp2 << std::endl;
//	std::cout << "betadotp : " << betadotp << std::endl;
//	std::cout << "DHeff : " << DHeff << std::endl;
//	std::cout << betadotp << "   " << eoba << "   " << eobd << "   " << Heff << std::endl;
//	abort();
	return Heff;
}

double EOBspinHamiltonian::evaluate(std::vector<double> variables) const
{
	extern int hamilton_count;
	extern time_t hamilton_time;
	time_t track_time;
	hamilton_count += 1;
	track_time = time(NULL);
	copy_ptr<Flags>			flags_ptr	  = get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();

	double eta, Heff, Hreal;
	eta	= wave_params_ptr->get_eta();
	Heff 	= evaluate_Heff(variables);
	hamilton_time += time(NULL) - track_time;
	
	if (eta < 1.0e-16 || flags_ptr->get_flag_LR()==1)
		Hreal = Heff;
	else
		Hreal = sqrt(1 + 2 * eta * (Heff-1)) / eta;
	return Hreal;
}

double EOBspinHamiltonian::nonKeplerianRadius() const
{
	SphEOBspinHamiltonian SphspinH(*this);
	return SphspinH.nonKeplerianRadius();
}

double EOBspinHamiltonian::dHdr_r(double r) const
{
	std::cout << "Error! EOBspinHamiltonian::dHdr_r: ";
	std::cout << "member function not implemented for EOBspinHamiltonian, ";
	std::cout << "use SphEOBspinHamiltonian instead." << std::endl;
	abort();
	return 0.0;
}

double EOBspinHamiltonian::dHdr_pp(double pp) const
{
	std::cout << "Error! EOBspinHamiltonian::dHdr_pp: ";
	std::cout << "member function not implemented for EOBspinHamiltonian, ";
	std::cout << "use SphEOBspinHamiltonian instead." << std::endl;
	abort();
	return 0.0;
}

double EOBspinHamiltonian::LightRingRadius() const
{
	SphEOBspinHamiltonian SphspinH(*this);
	return SphspinH.LightRingRadius();
}

/**
 * Compute the unit vector orthogonal to the orbital plane for the dynamical configuration provided in the argument
 */
Array EOBspinHamiltonian::evaluate_LNhat(std::vector<double> variables /**< Array of the dynamical variables, in general rvec, pvec, and spin vectors */) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	
	unsigned int i;
	Array LNhat(3);	
	Array rvec(3), pvec(3), S1(3), S2(3), rdotvec(3);
	
	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1[i] = variables[i + 6];
		S2[i] = variables[i + 9];
	}
	
	//	for(i=0; i<12; i++)
	//		std::cout << variables[i] << std::endl;
	//Compute rdot vector
	rdotvec[0] = temphptr->evaluate_deriv(3);
	rdotvec[1] = temphptr->evaluate_deriv(4);
	rdotvec[2] = temphptr->evaluate_deriv(5);
	//	std::cout << rdotvec[0] << std::endl;
	//	std::cout << rdotvec[1] << std::endl;
	//	std::cout << rdotvec[2] << std::endl;
	
	//Build LNhat
	LNhat =  crossproduct(rvec,rdotvec);
	LNhat = LNhat / sqrt(innerproduct(LNhat,LNhat));	
	return LNhat;
}


//**********************************************************************//
//
// SphEOBspinHamiltonian member functions:
//
//**********************************************************************//

SphEOBspinHamiltonian::SphEOBspinHamiltonian() : EOBspinHamiltonian(),LNhat(3)
{
}
/**
 * Constructor of the spherical Hamiltonian corresponding to the Cartesian Hamiltonian
 * of Phys. Rev. D 81, 084041 (2010). It takes the stored dynamical configuration and rearranges it
 * in spherical components
 */
SphEOBspinHamiltonian::SphEOBspinHamiltonian(const EOBspinHamiltonian & targetH /**< Cartesian spin Hamiltonian object */)
: EOBspinHamiltonian(),LNhat(3)
{	
	int i,j;
	Array rvec(3), pvec(3), S1(3), S2(3), Xhat(3);
	//	Array Yhat(3,0.), Zhat(3,0.);
	Array rvecprime(3), pvecprime(3), S1prime(3), S2prime(3);
	Array Xprime(3);
	Array Yprime(3), Zprime(3);
	std::vector<double> variables = targetH.get_variables();
	
	// Fixed Cartesian frame
	Xhat[0]=1.;
	Xhat[1]=0.;
	Xhat[2]=0.;
	
	//	Yhat[0]=0.;
	//	Yhat[1]=1.;
	//	Yhat[2]=0.;
	//	
	//	Zhat[0]=0.;
	//	Zhat[1]=0.;
	//	Zhat[2]=1.;
	
	//	std::cout << "Vars wrt fixed Cartesian frame" << std::endl;
	//	for(int i=0; i<12;i++)
	//		std::cout << variables[i] << std::endl;
	//	std::cout << " " << std::endl;
	
	//	Organize the vars
	rvec[0]    = variables[0];
	rvec[1]    = variables[1];
	rvec[2]    = variables[2];
	
	pvec[0]    = variables[3];
	pvec[1]    = variables[4];
	pvec[2]    = variables[5];
	
	S1[0]      = variables[6];
	S1[1]      = variables[7];
	S1[2]      = variables[8];
	
	S2[0]      = variables[9];
	S2[1]      = variables[10];
	S2[2]      = variables[11];
	
	//	std::cout << "S2 in constructor before rotation" << variables[9]<< std::endl;
	
	//LNhat is XIhat
	LNhat = targetH.evaluate_LNhat(targetH.get_variables());
	//	std::cout << "LNhat " << LNhat << std::endl;
	
	Xprime = LNhat;
	//	std::cout << "LNhat " << Xprime << std::endl;
	for (i=0; i<3; i++)
		if (fabs(Xprime[i])<1.e-8) 
		{
			Xprime[i]=0.;
		}
	for (i=0; i<3; i++)
		Xprime = Xprime/sqrt(innerproduct(Xprime,Xprime));
	//	std::cout << "LNhat " << Xprime << std::endl;
	
	/////////////////////////////////////////////////////////////////
	//	NOT USED
	//	Compute polar angles for LNhat wrt to a polar axis along zhat	
	//	double thetaN = acos(innerproduct(Xprime,zvec));
	//	double phiN = atan2(innerproduct(Xprime,yvec),innerproduct(Xprime,xvec));
	//	std::cout << "thetaN " << thetaN << std::endl;
	//	std::cout << "phiN " << phiN << std::endl;
	
	//	RotationMatrix[0][0] = sin(thetaN)*cos(phiN);
	//	RotationMatrix[0][1] = -sin(phiN);
	//	RotationMatrix[0][2] = -cos(thetaN)*cos(phiN);
	//	RotationMatrix[1][0] = sin(thetaN)*sin(phiN);
	//	RotationMatrix[1][1] = cos(phiN);
	//	RotationMatrix[1][2] = -cos(thetaN)*sin(phiN);
	//	RotationMatrix[2][0] = cos(thetaN);
	//	RotationMatrix[2][1] = 0.;
	//	RotationMatrix[2][2] = sin(thetaN);
	
	//	InvRotationMatrix[0][0] = sin(thetaN)*cos(phiN);
	//	InvRotationMatrix[0][1] = sin(thetaN)*sin(phiN);
	//	InvRotationMatrix[0][2] = cos(thetaN);
	//	InvRotationMatrix[1][0] = -sin(phiN);
	//	InvRotationMatrix[1][1] = cos(phiN);
	//	InvRotationMatrix[1][2] = 0.;
	//	InvRotationMatrix[2][0] = -cos(thetaN)*cos(phiN);
	//	InvRotationMatrix[2][1] = -cos(thetaN)*sin(phiN);
	//	InvRotationMatrix[2][2] = sin(thetaN);
	////////////////////////////////////////////////////////////////
	
	//	Define the basis-transofrmation matrix from the fixed Cartesian
	//	frame to the inst. orb. plane Cartesian frame
	//	Xprime = RotationMatrix.Xhat
	//	Yprime = RotationMatrix.Yhat
	//	Zprime = RotationMatrix.Zhat
	double a=Xprime[0], b=Xprime[1], c=Xprime[2];
	//	std::vector< std::vector<double> > RotationMatrix(3,3), InvRotationMatrix(3,3);
	//	
	//	if(Xprime[0]==1. && Xprime[1]==0. && Xprime[2]==0.)
	//	{
	//		RotationMatrix[0][0] = 1.;
	//		RotationMatrix[0][1] = 0.;
	//		RotationMatrix[0][2] = 0.;
	//		RotationMatrix[1][0] = 0.;
	//		RotationMatrix[1][1] = 1.;
	//		RotationMatrix[1][2] = 0.;
	//		RotationMatrix[2][0] = 0.;
	//		RotationMatrix[2][1] = 0.;
	//		RotationMatrix[2][2] = 1.;
	//	}
	//	else
	//	{
	//		RotationMatrix[0][0] = a;
	//		RotationMatrix[0][1] = sqrt(1-a*a);
	//		RotationMatrix[0][2] = 0.;
	//		RotationMatrix[1][0] = b;
	//		RotationMatrix[1][1] = -a*b/RotationMatrix[0][1];
	//		RotationMatrix[1][2] = c/RotationMatrix[0][1];
	//		RotationMatrix[2][0] = c;
	//		RotationMatrix[2][1] = -a*c/RotationMatrix[0][1];
	//		RotationMatrix[2][2] = -b/RotationMatrix[0][1];
	//	}
	
	std::vector< std::vector<double> > RotationMatrix(3), InvRotationMatrix;
	
	if(Xprime[0]==1. && Xprime[1]==0. && Xprime[2]==0.)
	{
		(RotationMatrix[0]).push_back(1.);
		(RotationMatrix[0]).push_back(0.);
		(RotationMatrix[0]).push_back(0.);
		(RotationMatrix[1]).push_back(0.);
		(RotationMatrix[1]).push_back(1.);
		(RotationMatrix[1]).push_back(0.);
		(RotationMatrix[2]).push_back(0.);
		(RotationMatrix[2]).push_back(0.);
		(RotationMatrix[2]).push_back(1.);
	}
	else
	{
		(RotationMatrix[0]).push_back(a);
		(RotationMatrix[0]).push_back(sqrt(1-a*a));
		(RotationMatrix[0]).push_back(0.);
		(RotationMatrix[1]).push_back(b);
		(RotationMatrix[1]).push_back(-a*b/RotationMatrix[0][1]);
		(RotationMatrix[1]).push_back(c/RotationMatrix[0][1]);
		(RotationMatrix[2]).push_back(c);
		(RotationMatrix[2]).push_back(-a*c/RotationMatrix[0][1]);
		(RotationMatrix[2]).push_back(-b/RotationMatrix[0][1]);
	}
	
	//	std::cout << "Rot matrix" << std::endl;
	//	for (i=0; i<3; i++)
	//		for(j=0; j<3; j++)
	//			std::cout << RotationMatrix[i][j] << std::endl;
	
	
	
	//	Chop rotation matrix
	//	for (i=0; i<3; i++)
	//		for(j=0; j<3; j++)
	//			if (abs(RotationMatrix[i][j])<1.e-15) RotationMatrix[i][j]=0.;
	
	//  RotationMatrix is orthogonal so the inverse is just its transpose	
	for(int i=0; i<3; i++)
		InvRotationMatrix.push_back(std::vector<double>(3));
	for (i=0; i<3; i++)
		for(j=0; j<3; j++)
			InvRotationMatrix[i][j] = RotationMatrix[j][i];
	
	//  Rotate all the vectors and bring them in the inst. orb. plane frame	
	for (i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			//			std::cout << "InvRotationMatrix " << InvRotationMatrix[i][j] << std::endl;
			rvecprime[i] += InvRotationMatrix[i][j]*rvec[j];
			pvecprime[i] += InvRotationMatrix[i][j]*pvec[j];
			S1prime[i] += InvRotationMatrix[i][j]*S1[j];
			S2prime[i] += InvRotationMatrix[i][j]*S2[j];
		}
	variables[0] = rvecprime[0];
	variables[1] = rvecprime[1];
	variables[2] = rvecprime[2];
	
	variables[3] = pvecprime[0];
	variables[4] = pvecprime[1];
	variables[5] = pvecprime[2];
	
	variables[6]  = S1prime[0];
	variables[7]  = S1prime[1];
	variables[8]  = S1prime[2];
	
	variables[9]  = S2prime[0];
	variables[10] = S2prime[1];
	variables[11] = S2prime[2];
	
	//	std::cout << "S2 in constructor after rotation" << variables[9]<< std::endl;
	
	//	std::cout << "Vars after rotation" << std::endl;
	//	for(int i=0; i<12;i++)
	//		std::cout << variables[i] << std::endl;
	//	std::cout << " " << std::endl;
	
	//	Move to spherical coordinates wrt to a polar axis along {1,0,0} which is LNhat
	variables[0] = sqrt(innerproduct(rvecprime,rvecprime));
	variables[1] = acos(rvecprime[0] / variables[0]);
	variables[2] = atan2(-rvecprime[1], rvecprime[2]);
	
	variables[3] = innerproduct(rvecprime, pvecprime) / variables[0];
	variables[4] = -innerproduct(crossproduct(crossproduct(rvecprime, Xhat), rvecprime), pvecprime) 
	/ variables[0] / sin(variables[1]);
	variables[5] = -innerproduct(crossproduct(rvecprime, Xhat), pvecprime);
	
	//	std::cout << "List of SphH variables:" << std::endl;
	//	for (int i = 0; i < 12; i++)
	//		std::cout << variables[i] << "   " << std::endl;
	//	std::cout << " " << std::endl;
	
	/// OLD IMPLEMENTATION /////////////////
	//	//Define instantaneous orbital plane ref frame
	//	double r = sqrt(innerproduct(rvec,rvec));
	//	
	//	if(Xprime[0]==1. && Xprime[1]==0. && Xprime[2]==0.)
	//	{
	//		zvec[0]=0.;
	//		zvec[1]=1.;
	//		zvec[2]=0.;
	//	}
	//	else
	//	{
	//		zvec[0]=1.;
	//		zvec[1]=0.;
	//		zvec[2]=0.;
	//	}
	//	
	//	Zprime = crossproduct(Xprime,zvec)+0.*rvec/r;
	//	Yprime = crossproduct(Zprime,Xprime);
	//	
	////	std::cout << "Instantaneous orbital plane frame" << std::endl;
	//	std::cout << "Xprime = " << Xprime << std::endl;
	//	std::cout << "Yprime = " << Yprime << std::endl;
	//	std::cout << "Zprime = " << Zprime << std::endl;
	////	abort();
	//	
	//	//Build the variables in the instantaneous orbital plane frame
	//	rvecprime[0] = innerproduct(rvec,Xprime);
	//	rvecprime[1] = innerproduct(rvec,Yprime);
	//	rvecprime[2] = innerproduct(rvec,Zprime);
	//	
	//	pvecprime[0] = innerproduct(pvec,Xprime);
	//	pvecprime[1] = innerproduct(pvec,Yprime);
	//	pvecprime[2] = innerproduct(pvec,Zprime);
	//	
	//	S1prime[0] = innerproduct(S1,Xprime);
	//	S1prime[1] = innerproduct(S1,Yprime);
	//	S1prime[2] = innerproduct(S1,Zprime);
	//	
	//	S2prime[0] = innerproduct(S2,Xprime);
	//	S2prime[1] = innerproduct(S2,Yprime);
	//	S2prime[2] = innerproduct(S2,Zprime);
	//	
	//	//Arrange everything in the variables array
	//	variables[0] = rvecprime[0];
	//	variables[1] = rvecprime[1];
	//	variables[2] = rvecprime[2];
	//	
	//	variables[3] = pvecprime[0];
	//	variables[4] = pvecprime[1];
	//	variables[5] = pvecprime[2];
	//	
	//	variables[6]  = S1prime[0];
	//	variables[7]  = S1prime[1];
	//	variables[8]  = S1prime[2];
	//	
	//	variables[9]  = S2prime[0];
	//	variables[10] = S2prime[1];
	//	variables[11] = S2prime[2];
	//	
	//	std::cout << "New vars" << std::endl;
	//	for(int i=0; i<12;i++)
	//		std::cout << variables[i] << std::endl;
	//	
	//	rvec[0] =  variables[2];
	//	rvec[1] = -variables[1];
	//	rvec[2] =  variables[0];
	//	pvec[0] =  variables[5];
	//	pvec[1] = -variables[4];
	//	pvec[2] =  variables[3];
	//	
	//	variables[0] = r;
	//	variables[1] = acos(rvec[2] / variables[0]);
	//	variables[2] = atan2(rvec[1], rvec[0]);
	//	variables[3] = innerproduct(rvec, pvec) / variables[0];
	//	variables[4] = -innerproduct(crossproduct(crossproduct(rvec, zvec), rvec), pvec) 
	//	/ variables[0] / sin(variables[1]);
	//	variables[5] = -innerproduct(crossproduct(rvec, zvec), pvec);
	//	
	//	std::cout << "List of SphH variables:" << std::endl;
	//	for (int i = 0; i < 6; i++)
	//		std::cout << variables[i] << "   " << std::endl;
	//	abort();
	//std::cout << crossproduct(rvec, zvec) << std::endl;
	//std::cout << std::endl;
	//	abort();
	/////////////////////////////////////////////////////////
	//	std::cout << "S2 in constructor" << variables[9]<< std::endl;
	
	set_flags_ptr(targetH.get_flags_ptr());
	set_wave_params_ptr(targetH.get_wave_params_ptr());
	set_tuning_params_ptr(targetH.get_tuning_params_ptr());
	set_variables(variables);
	set_saves(targetH.get_saves());
}


SphEOBspinHamiltonian::~SphEOBspinHamiltonian()
{
}

/**
 * Compute r for the stored dynamical configuration
 */
double SphEOBspinHamiltonian::evaluate_radius() const
{
	std::vector<double> variables = get_variables();
	return variables[0];
}

double SphEOBspinHamiltonian::evaluate(std::vector<double> variables) const
{
	int i,j;
	double r, th, ph, pr, pt, pp;
	std::vector<double> cartvar(variables.size());
	
	
	r	= variables[0];
	th	= variables[1];
	ph	= variables[2];
	pr	= variables[3];
	pt	= variables[4];
	pp	= variables[5];
	for (i = 6; i < 12; i++)
		cartvar[i] = variables[i];
	//	std::cout << "S2 in eval " << variables[9]<< std::endl;
	
	//	Assume the polar axis along {1,0,0} and move to Cartesian coordinates:
	//	these are Cart coords of the instantaneous orb. plane frame
	cartvar[0] = r*cos(th);
	cartvar[1] =-r*sin(th)*sin(ph);
	cartvar[2] = r*sin(th)*cos(ph);
	
	if(th==0. || th==Pi)
	{
		if(th==0.)
		{
			th = Pi/2;
			ph = 0.;
		}
		else
		{
			th = Pi/2;
			ph = Pi;
		}
		cartvar[3] = pr*sin(th)*cos(ph) + pt/r*cos(th)*cos(ph) - pp/r/sin(th)*sin(ph);
		cartvar[4] = pr*sin(th)*sin(ph) + pt/r*cos(th)*sin(ph) + pp/r/sin(th)*cos(ph);
		cartvar[5] = pr*cos(th)         - pt/r*sin(th);		
	}
	else
	{
		cartvar[3] = pr*cos(th) 		-pt/r*sin(th);
		cartvar[4] =-pr*sin(th)*sin(ph) -pt/r*cos(th)*sin(ph) -pp/r/sin(th)*cos(ph);
		cartvar[5] = pr*sin(th)*cos(ph) +pt/r*cos(th)*cos(ph) -pp/r/sin(th)*sin(ph);
	}
	
	//	Vars in the Cart inst orb frame
	//	std::cout << "Vars in the Cart inst orb frame" << std::endl;
	//	for(int i=0; i<12;i++)
	//		std::cout << cartvar[i] << std::endl;
	//	abort();
	
	//  Rotate to the fixed Cartesian frame
	double a=LNhat[0], b=LNhat[1], c=LNhat[2];
	//	std::cout << "LNhat " << LNhat << std::endl;
	
	//	std::vector< std::vector<double> > RotationMatrix(3,3);
	//	if(LNhat[0]==1. && LNhat[1]==0. && LNhat[2]==0.)
	//	{
	//		RotationMatrix[0][0] = 1.;
	//		RotationMatrix[0][1] = 0.;
	//		RotationMatrix[0][2] = 0.;
	//		RotationMatrix[1][0] = 0.;
	//		RotationMatrix[1][1] = 1.;
	//		RotationMatrix[1][2] = 0.;
	//		RotationMatrix[2][0] = 0.;
	//		RotationMatrix[2][1] = 0.;
	//		RotationMatrix[2][2] = 1.;
	//	}
	//	else
	//	{
	//		RotationMatrix[0][0] = a;
	//		RotationMatrix[0][1] = sqrt(1-a*a);
	//		RotationMatrix[0][2] = 0.;
	//		RotationMatrix[1][0] = b;
	//		RotationMatrix[1][1] = -a*b/RotationMatrix[0][1];
	//		RotationMatrix[1][2] = c/RotationMatrix[0][1];
	//		RotationMatrix[2][0] = c;
	//		RotationMatrix[2][1] = -a*c/RotationMatrix[0][1];
	//		RotationMatrix[2][2] = -b/RotationMatrix[0][1];
	//	}
	
	std::vector< std::vector<double> > RotationMatrix(3);
	//	RotationMatrix.push_back(std::vector<double>());
	
	if(LNhat[0]==1. && LNhat[1]==0. && LNhat[2]==0.)
	{
		(RotationMatrix[0]).push_back(1.);
		(RotationMatrix[0]).push_back(0.);
		(RotationMatrix[0]).push_back(0.);
		(RotationMatrix[1]).push_back(0.);
		(RotationMatrix[1]).push_back(1.);
		(RotationMatrix[1]).push_back(0.);
		(RotationMatrix[2]).push_back(0.);
		(RotationMatrix[2]).push_back(0.);
		(RotationMatrix[2]).push_back(1.);
	}
	else
	{
		(RotationMatrix[0]).push_back(a);
		(RotationMatrix[0]).push_back(sqrt(1-a*a));
		(RotationMatrix[0]).push_back(0.);
		(RotationMatrix[1]).push_back(b);
		(RotationMatrix[1]).push_back(-a*b/RotationMatrix[0][1]);
		(RotationMatrix[1]).push_back(c/RotationMatrix[0][1]);
		(RotationMatrix[2]).push_back(c);
		(RotationMatrix[2]).push_back(-a*c/RotationMatrix[0][1]);
		(RotationMatrix[2]).push_back(-b/RotationMatrix[0][1]);
	}
	
	
	//	//	Chop rotation matrix
	//	for (i=0; i<3; i++)
	//		for(j=0; j<3; j++)
	//			if (abs(RotationMatrix[i][j])<1.e-15) RotationMatrix[i][j]=0.;
	
	Array rvec(3), rvecprime(3);
	Array pvec(3), pvecprime(3);
	Array S1(3), S1prime(3);
	Array S2(3), S2prime(3);
	
	rvecprime[0] = cartvar[0];
	rvecprime[1] = cartvar[1];
	rvecprime[2] = cartvar[2];
	
	pvecprime[0] = cartvar[3];
	pvecprime[1] = cartvar[4];
	pvecprime[2] = cartvar[5];
	
	S1prime[0] = cartvar[6];
	S1prime[1] = cartvar[7];
	S1prime[2] = cartvar[8];
	
	S2prime[0] = cartvar[9];
	S2prime[1] = cartvar[10];
	S2prime[2] = cartvar[11];
	
	
	for (i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			rvec[i] += RotationMatrix[i][j]*rvecprime[j];
			pvec[i] += RotationMatrix[i][j]*pvecprime[j];
			S1[i] += RotationMatrix[i][j]*S1prime[j];
			S2[i] += RotationMatrix[i][j]*S2prime[j];
		}
	
	cartvar[0] = rvec[0];
	cartvar[1] = rvec[1];
	cartvar[2] = rvec[2];
	cartvar[3] = pvec[0];
	cartvar[4] = pvec[1];
	cartvar[5] = pvec[2];		
	cartvar[6] = S1[0];
	cartvar[7] = S1[1];
	cartvar[8] = S1[2];	
	cartvar[9] =  S2[0];
	cartvar[10] = S2[1];
	cartvar[11] = S2[2];	
	
	//	for (i = 0; i < 12; i++)
	//		std::cout << cartvar[i] << std::endl;
	//	abort();
	return EOBspinHamiltonian::evaluate(cartvar);
}

//double SphEOBspinHamiltonian::evaluate(std::vector<double> variables) const
//{
//	int i;
//	double r, th, ph, pr, pt, pp;
//	std::vector<double> cartvar(variables.size());
//
//	r	= variables[0];
//	th	= variables[1];
//	ph	= variables[2];
//	pr	= variables[3];
//	pt	= variables[4];
//	pp	= variables[5];
//	for (i = 6; i < 12; i++)
//		cartvar[i] = variables[i];
//
//	cartvar[2] = r*sin(th)*cos(ph);
//	cartvar[1] =-r*sin(th)*sin(ph);
//	cartvar[0] = r*cos(th);
//	cartvar[5] = pr*sin(th)*cos(ph) +pt/r*cos(th)*cos(ph) -pp/r/sin(th)*sin(ph);
//	cartvar[4] =-pr*sin(th)*sin(ph) -pt/r*cos(th)*sin(ph) -pp/r/sin(th)*cos(ph);
//	cartvar[3] = pr*cos(th) 		-pt/r*sin(th);
//
//	//for (i = 0; i < 12; i++)
//	//	std::cout << cartvar[i] << std::endl;
//
//	return EOBspinHamiltonian::evaluate(cartvar);
//}


double SphEOBspinHamiltonian::nonKeplerianRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	return temphptr -> nonKeplerianRadiusSph();	
}

double SphEOBspinHamiltonian::dHdr_r(double r) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	variables[0] = r;	
	variables[3] = 0.0;
	temphptr->set_variables(variables);
	return temphptr->evaluate_deriv(0);
}

double SphEOBspinHamiltonian::dHdr_pp(double pp) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	variables[3] = 0.0;
	variables[5] = pp;
	temphptr->set_variables(variables);
	return temphptr->evaluate_deriv(0);
}


double SphEOBspinHamiltonian::LightRingRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	return temphptr -> LightRingSph();
}

//**********************************************************************//
//
// BBHamiltonian member functions:
//
//**********************************************************************//
/**
 * Cartesian Hamiltonian of Phys.Rev. D81 (2010) 084024
 */
BBHamiltonian::BBHamiltonian() : Hamiltonian(12, 8)
{
}

BBHamiltonian::~BBHamiltonian()
{
}


void BBHamiltonian::generate_saves()
{
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> saves(8);
	double eta;
	eta	= wave_params_ptr->get_eta();

	//Mass saves
	saves[0] = sqrt(1. - 4. * eta);					//dm
	saves[1] = (1. + saves[0]) / (1. - saves[0]);			//m1om2
	saves[2] = 1. / saves[1]; //m2om2
	
	//K(eta) polynomial saves
	//saves[3]  = 1.4467 - 2.7868 *eta; 
	saves[3]  = tuning_params_ptr->get_K();
	//	saves[3] = 0.75;	
	
	//Tuning saves
	saves[4]  = tuning_params_ptr->get_Ka();//spin term in K(eta)
    saves[5]  = tuning_params_ptr->get_wfd1();//wfd1 in BB paper
	saves[6]  = tuning_params_ptr->get_wfd2();//wfd2 in BB paper
	saves[7]  = tuning_params_ptr->get_d1();//wfd2 in BB paper
	set_saves(saves);
	return;
}

/**
 * Compute r for the stored dynamical configuration
 */
double BBHamiltonian::evaluate_radius() const
{
	std::vector<double> variables = get_variables();	
	return sqrt(variables[0]*variables[0]+variables[1]*variables[1]+variables[2]*variables[2]);
}

/**
 * Compute the orbital freq for the stored dynamical configuration
 */
double BBHamiltonian::evaluate_orb_omega() const
{
	std::vector<double> variables = get_variables();	
	Array rvec(3), rdvec(3), omegavec(3);
	
	rvec[0] = variables[0];
	rvec[1] = variables[1];
	rvec[2] = variables[2];
	
	for (int i = 0; i < 3; i++)
	{
		rdvec[i] = this->evaluate_deriv(i + 3);		
	}
	
	omegavec = crossproduct(rvec, rdvec) / innerproduct(rvec, rvec);
	return sqrt(innerproduct(omegavec, omegavec));
}

double BBHamiltonian::evaluate_eoba(std::vector<double> variables) const
{
	std::vector<double> saves = get_saves();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	int i;
	double eoba, eta, r, bulk, LogTerms, u;
	double KK, a, deltaU;
//	double KK0, c0, c1, c2;
	double k0, k1, k2, k3, k4;
	Array rvec(3);
	std::vector<double> SKerr(3);
	
	eta  = wave_params_ptr->get_eta();
	
   	for (i = 0; i < 3; i++)
		rvec[i] = variables[i];
	r  = sqrt(innerproduct(rvec, rvec)); 
	u = 1./r;
	
	SKerr = evaluate_sigmaKerr(variables);
	a = sqrt(SKerr[0]*SKerr[0]+SKerr[1]*SKerr[1]+SKerr[2]*SKerr[2]);
	
	//Tuning saves
//	KK0  = saves[0];
//	Ka   = saves[4];
//	wfd1 = saves[5];
//	wfd2 = saves[6];
	
	//K(eta) polynomial
//	c0  = saves[4];
//	c1  = -8.*(c0 - KK0 + Ka*a*a);//c1 or K_1 in BB paper Eq. 6.11 
//	c2  = 16.*(c0 - KK0 + Ka*a*a);//c2 or K_2 in BB paper Eq. 6.11 
//	KK  = c0 + c1*eta + c2*eta*eta;//KK or K in BB paper Eq. 6.9
	KK  = saves[3] + saves[4]*a*a;
	
	//\Delta_i's
	k0  = KK*(-2. + eta*KK);//k0 or \Delta_1 in BB paper Eq. 5.77
	k1  = -2.*(k0 + KK)*(-1. + eta*KK);//k1 or \Delta_2 in BB paper Eq. 5.78
    k2  = (k1*(4. + k1 - 4.*eta*KK))/2. - intPower(a,2)*k0*intPower(-1. + eta*KK,2);
    k3  = -intPower(k1,3)/3. + k1*k2 + intPower(k1,2)*(-1. + eta*KK) - 2.*(1. + k2 - eta*KK)*(-1. + eta*KK) - intPower(a,2)*k1*intPower(-1. + eta*KK,2);
    k4  = (24.*intPower(k1,4) - 96.*intPower(k1,2)*k2 + 48.*intPower(k2,2) - 64.*intPower(k1,3)*(-1. + eta*KK) 
		   + 48.*intPower(a,2)*(intPower(k1,2) - 2.*k2)*intPower(-1. + eta*KK,2) + 
		   96.*k1*(k3 + 2.*k2*(-1. + eta*KK)) - (-1. + eta*KK)*(192.*k3 + (-1. + eta*KK)*(-3008. + 123.*intPower(Pi,2))))/96.;  
	bulk = intPower(-1. + eta*KK,-2) + (2.*u)/(-1. + eta*KK) + intPower(a,2)*intPower(u,2);	
    LogTerms = ( 1. + eta*k0 + eta*log(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4)) );
    deltaU = bulk*LogTerms;
	eoba = deltaU;
//	std::cout << k0 << " " << k1 << " " << k2 << " " << k3 << " " << k4 << std::endl;
//	std::cout << "bulk "<< bulk << std::endl;
//	std::cout << "LogTerms " << LogTerms << std::endl;
//	std::cout << "KK " << KK << std::endl;
//	std::cout << "a " << a << std::endl;
//	std::cout << "u " << u << std::endl;
//	std::cout << "eta " << eta << std::endl;
//	std::cout << "EOBA " << eoba << std::endl;

/////////////////////////////////////////////	
//  Tidal correction from Tanja May 11, 2012
//	double Lambda1;
//	double Lambda2;
//	double m1 = wave_params_ptr->get_m1();
//	double m2 = wave_params_ptr->get_m2();
////	NSBH
////  Lambda1 = 0.;	
////	Lambda2 = 0.784965; //NSBH q=5, chi=0.5
////	Lambda2 = 0.260786; //NSBH q=7, chi=0.5
//
////	NSNS
//	double lambda21;
//	double lambda22;
//	lambda21 = 24.3;  // compactness = 0.145
////	lambda21 = 8.216; // compactness = 0.17
////	lambda21 = 2.426; // compactness = 0.2
//	lambda22 = lambda21;
//	Lambda1 = 3.*lambda21*m2/2./m1;
//	Lambda2 = 3.*lambda22*m1/2./m2;
//	
////	Adiabatic quadrupolar effects
//	eoba = eoba - 2.*intPower(u,6)*(Lambda2*(1. + 2.5*u*m2 + u*u*(337./28.*m2*m2 + m2/8. + 3.))
//	+ Lambda1*(1. + 2.5*u*m1 + u*u*(337./28.*m1*m1 + m1/8. + 3.)));
//	
////	Octupole
//	double lambda31, lambda32;
//	lambda31 = 16.91;  // compactness = 0.145
////	lambda31 = 3.9; // compactness = 0.17
////	lambda31 = 0.77; // compactness = 0.2
//	lambda32 = lambda31;
//	eoba = eoba - 15.*(m1/m2*lambda31 + m2/m1*lambda32)*intPower(u,7);
////
////  f-mode
//	double omega0 = 0.175; // for equal masses it should be in [0.146,0.2]
//	double W = 4.*u*u*u/omega0/omega0;
//	eoba = eoba + 2.*intPower(u,6)/(1. - W)*(Lambda2*((4. - W)/4. + u/(1. - W)*(2.5*m2 + 
//		   W/8.*(-31. + W*(11. - 5.*m1) + 43.*m1))) + Lambda1*((4. - W)/4. 
//		   + u/(1. - W)*(2.5*m2 + W/8.*(-31. + W*(11. - 5.*m1) + 43.*m1))));
/////////////////////////////////////////////	
	return eoba;
}


//////////////////////////////////////////////////////////////////////////
// Non-spinning A(r) Pade (1,5) potential, the same as Delta_u since there
// is no spin
//
//double BBHamiltonian::evaluate_eoba(std::vector<double> variables) const
//{
//	std::vector<double> saves = get_saves();
//	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
//	int i;
//	double eoba, eta, r, u;
//	double a5, a6;
//	Array rvec(3);
//	
//	eta  = wave_params_ptr->get_eta();
//	
//   	for (i = 0; i < 3; i++)
//		rvec[i] = variables[i];
//	r  = sqrt(innerproduct(rvec, rvec)); 
//	u = 1./r;
//	
//	a5 = -82.5384 + 508.681*eta - 787.826*eta*eta;
//	a6 = 500 - 1800*eta;
//
//	eoba =  (96*(768 + (-3584 - 24*a5 + 123*Pi*Pi)*eta + 
//				 3*u*(-512 + (3520 + 32*a5 + 8*a6 - 123*Pi*Pi - 32*eta)*eta)))/
//	(73728 + eta*(2304*a5*a5*u*u*u*u*u*eta + 15129*Pi*Pi*Pi*Pi*u*u*u*u*(1 + u)*eta + 
//				24*a5*(-96 + 192*u*(1 + 2*u)*(1 + 4*u*u) + 
//			u*u*u*(192 + u*(-123*Pi*Pi*(1 + 4*u) + 64*(47 + 212*u)))*eta) + 
//				  256*(-1344 + 3*u*(8*(53 + 2*u*(53 + 94*u)) + 
//				3*a6*(1 + 2*u*(1 + 2*u)*(1 + 4*u*u))) - 
//				2*u*(18 + u*(36 + u*(-1272 + u*(-22328 + 9*a6 + 59*
//				(-376 + 3*a6)*u))))*eta + 24*u*u*u*u*(3 + 59*u)*
//				eta*eta) - 984*Pi*Pi*(-12 + 
//	u*(12 + u*(24 + u*(48 + eta*(24 + u*(848 + u*(848 - 3*a6 + 12*eta)))))))));
//	
//	return eoba;
//}
//////////////////////////////////////////////////////////////////////////


double BBHamiltonian::evaluate_eoba_deriv(std::vector<double> variables) const
{
	std::vector<double> saves = get_saves();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	int i;
	double eta, r, bulk, LogTerms, u;
	double a;
	//double KK0, c0, c1, c2;
	double KK, k0, k1, k2, k3, k4;
	Array rvec(3);
	std::vector<double> SKerr(3);
	
	eta  = wave_params_ptr->get_eta();
	
   	for (i = 0; i < 3; i++)
		rvec[i] = variables[i];
	r  = sqrt(innerproduct(rvec, rvec)); 
	u = 1./r;
	
	SKerr = evaluate_sigmaKerr(variables);
	a = sqrt(SKerr[0]*SKerr[0]+SKerr[1]*SKerr[1]+SKerr[2]*SKerr[2]);
	
	//Tuning saves
//	KK0  = saves[0];
//	Ka   = saves[4];
//	wfd1 = saves[5];
//	wfd2 = saves[6];
	
	//K(eta) polynomial
//	c0  = saves[4];
//	c1  = -8.*(c0 - KK0 + Ka*a*a);//c1 or K_1 in BB paper Eq. 6.11 
//	c2  = 16.*(c0 - KK0 + Ka*a*a);//c2 or K_2 in BB paper Eq. 6.11 
//	KK  = c0 + c1*eta + c2*eta*eta;//KK or K in BB paper Eq. 6.9
	KK = saves[3] + saves[4]*a*a;
	
	//\Delta_i's
	k0  = KK*(-2. + eta*KK);//k0 or \Delta_1 in BB paper Eq. 5.77
	k1  = -2.*(k0 + KK)*(-1. + eta*KK);//k1 or \Delta_2 in BB paper Eq. 5.78
    k2  = (k1*(4. + k1 - 4.*eta*KK))/2. - intPower(a,2)*k0*intPower(-1. + eta*KK,2);
    k3  = -intPower(k1,3)/3. + k1*k2 + intPower(k1,2)*(-1. + eta*KK) - 2.*(1. + k2 - eta*KK)*(-1. + eta*KK) - intPower(a,2)*k1*intPower(-1. + eta*KK,2);
    k4  = (24.*intPower(k1,4) - 96.*intPower(k1,2)*k2 + 48.*intPower(k2,2) - 64.*intPower(k1,3)*(-1. + eta*KK) 
		   + 48.*intPower(a,2)*(intPower(k1,2) - 2.*k2)*intPower(-1. + eta*KK,2) + 
		   96.*k1*(k3 + 2.*k2*(-1. + eta*KK)) - (-1. + eta*KK)*(192.*k3 + (-1. + eta*KK)*(-3008. + 123.*intPower(Pi,2))))/96.;  
	
	
	bulk = intPower(-1. + eta*KK,-2) + (2.*u)/(-1. + eta*KK) + intPower(a,2)*intPower(u,2);	
    LogTerms = ( 1. + eta*k0 + eta*log(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4)) );
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Implementation where eoba_deriv is deoba/dr
	double dbulkdr, dLogTermsdr, ddeltaUdr, dAdr;
	dbulkdr =  - (2.*u*u)/(-1. + eta*KK) - 2.*intPower(a,2)*intPower(u,3);
	dLogTermsdr = - eta*u*u*(k1 + u*(2.*k2 + u*(3.*k3 + 4.*k4*u)))/(1. + k1*u + k2*u*u + k3*u*u*u + k4*u*u*u*u);
	ddeltaUdr = dbulkdr*LogTerms + bulk*dLogTermsdr;
	dAdr = ddeltaUdr;
	
/////////////////////////////////////////////	
//  Tidal correction from Tanja Dec 8, 2011
//	double Lambda1 = 0.;
//	double Lambda2;
//	Lambda2 = 0.784965; //NSBH q=5, chi=0.5
//	Lambda2 = 0.260786; //NSBH q=7, chi=0.5
//	dAdr = dAdr + intPower(u,7)*(Lambda1*(12. + 35.*wave_params_ptr->get_m1()*u) + Lambda2*(12. + 35.*wave_params_ptr->get_m2()*u));
/////////////////////////////////////////////	

/////////////////////////////////////////////	
//  Tidal correction from Tanja May 11, 2012
//	double Lambda1;
//	double Lambda2;
//	double m1 = wave_params_ptr->get_m1();
//	double m2 = wave_params_ptr->get_m2();
////	NSBH
////  Lambda1 = 0.;	
////	Lambda2 = 0.784965; //NSBH q=5, chi=0.5
////	Lambda2 = 0.260786; //NSBH q=7, chi=0.5
//	
////	NSNS
//	double lambda21;
//	double lambda22;
//	lambda21 = 24.3;  // compactness = 0.145
////	lambda21 = 8.216; // compactness = 0.17
////	lambda21 = 2.426; // compactness = 0.2
//	lambda22 = lambda21;
//	Lambda1 = 3.*lambda21*m2/2./m1;
//	Lambda2 = 3.*lambda22*m1/2./m2;
//	
//	//	Adiabatic quadrupolar effects
//	dAdr = dAdr + intPower(u,7)*(Lambda2*(12. + 35.*u*m2 + 16.*u*u*(337./28.*m2*m2 + m2/8. + 3.))
//				+ Lambda1*(12. + 35.*u*m1 + 16.*u*u*(337./28.*m1*m1 + m1/8. + 3.)));
//	
////	//	Octupole
//	double lambda31, lambda32;
//	lambda31 = 16.91;  // compactness = 0.145
////	lambda31 = 3.9; // compactness = 0.17
////	lambda31 = 0.77; // compactness = 0.2
//	lambda32 = lambda31;
//	dAdr = dAdr + 105.*(m1/m2*lambda31 + m2/m1*lambda32)*intPower(u,8);
////	
////	//  f-mode
//	double omega0 = 0.175; // for equal masses it should be in [0.146,0.2]
//	double W = 4.*u*u*u/omega0/omega0;
//	dAdr = dAdr - intPower(u,7)/(1. - W)/(1. - W)/(1. - W)*
//	(Lambda1*(-12. + (22.5 + 77.5*u - 107.5*m2*u)*W + 
//			  (-13.5 - 66.75*u + 59.25*m2*u)*W*W + 
//			  (3. + 19.25*u - 8.75*m2*u)*W*W*W + m1*u*(-35. + 5.*W)) + 
//	 Lambda2*(-12. + (22.5 + 77.5*u - 107.5*m1*u)*W + 
//			  (-13.5 - 66.75*u + 59.25*m1*u)*W*W + 
//			  (3. + 19.25*u - 8.75*m1*u)*W*W*W + m2*u*(-35. + 5.*W))); 
/////////////////////////////////////////////		
	
	return dAdr;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//  Implementation where eoba_dervi is -u*u*deoba/du, matchinig Enrico's original code for deoba/du
//	double dbulkdu, dLogTermsdu, deobadu;
//	dbulkdu = 2.*a*a*u + 2./(-1. + KK*eta);
//	dLogTermsdu = eta*(k1 + u*(2.*k2 + u*(3.*k3 + 4.*k4*u)))/
//					(1. + k1*u + k2*u*u + k3*u*u*u + k4*u*u*u*u);
//	deobadu = bulk*dLogTermsdu + dbulkdu*LogTerms;
////	deobadu = 2.*(1./(-1. + eta*KK) + intPower(a,2)*u)*LogTerms + 
////	bulk* (eta*(k1 + u*(2.*k2 + u*(3.*k3 + 4.*k4*u))))/(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4));
//	return -deobadu/r/r;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////		
}


//////////////////////////////////////////////////////////////////////////
// Non-spinning dA/dr Pade (1,5) potential
//
//double BBHamiltonian::evaluate_eoba_deriv(std::vector<double> variables) const
//{
//	std::vector<double> saves = get_saves();
//	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
//	int i;
//	double dAdr, eta, r, bulk, LogTerms, u, u2, u3, u4, u5;
//	double a5, a52, a53, a6, a62, a63, Pi2, Pi4, Pi6;
//	double eta2, eta3, eta4, eta5;
//	Array rvec(3);
//	
//	eta  = wave_params_ptr->get_eta();
//	eta2 = eta*eta;
//	eta3 = eta2*eta; 
//	eta4 = eta3*eta;
//	eta5 = eta4*eta;
//	
//   	for (i = 0; i < 3; i++)
//		rvec[i] = variables[i];
//	r  = sqrt(innerproduct(rvec, rvec)); 
//	u = 1./r;
//	u2 = u*u;
//	u3 = u2*u;
//	u4 = u3*u;
//	u5 = u4*u;
//	
//	Pi2 = Pi*Pi;
//	Pi4 = Pi*Pi*Pi*Pi;
//	Pi6 = Pi4*Pi*Pi;
//	
//	a5 = -82.5384 + 508.681*eta - 787.826*eta*eta;
//	a6 = 500 - 1800*eta;
//
//	a52 = a5*a5;
//	a53 = a52*a5;
//	a62 = a6*a6;
//	
//	dAdr =(384*u2*(28311552 - 73728*(-123*Pi2*(1 - u - 2*u2 + 4*u3) + 
//									 24*a5*(1 - 2*u - 4*u2 - 8*u3 + 24*u4) + 
//									 8*(448 - (424 + 3*a6)*u - 2*(352 + 3*a6)*u2 - 
//										4*(-376 + 3*a6)*u3 - 24*a6*u4 + 96*a6*u5))*eta - 
//				   48*(15129*Pi4*(-1 + 2*u + 3*u2 - 28*u3 + 4*u4 + 32*u5) + 
//					   576*a52*(-1 + 4*u + 4*u2 - 176*u4 + 128*u5) + 
//					   1968*Pi2*(448 - (872 + 3*a6)*u - 3*(344 + a6)*u2 + 11488*u3 - 
//								 4*(424 + 9*a6)*u4 + 64*(-212 + 3*a6)*u5) - 
//					   64*(200704 - 128*(2986 + 21*a6)*u + 3*(-108864 - 944*a6 + 3*a62)*
//						   u2 + 4*(1170112 - 576*a6 + 9*a62)*u3 + 
//						   4*(-184384 - 8208*a6 + 27*a62)*u4 + 
//						   32*(-177472 + 5376*a6 + 9*a62)*u5) - 
//					   48*a5*(123*Pi2*(-1 + 3*u + 4*u2 - 12*u3 - 64*u4 + 64*u5) + 
//							  8*(448 - 3*(440 + a6)*u - 1504*u2 + 12*(328 + a6)*u3 + 
//								 16*(1768 + 3*a6)*u4 + 64*(-424 + 3*a6)*u5)))*eta2 + 
//				   u*(13824*a53*u3*(-5 + 16*u) - 1860867*Pi6*u2*(-1 + u + 3*u2) - 
//					  726192*Pi4*u*(3 + 210*u - (213 + a6)*u2 + (-644 + a6)*u3) - 
//					  576*a52*u*(144 - 41*(-64 + 3*Pi2)*u + (25856 - 861*Pi2)*u2 - 
//								 4*(24128 + 24*a6 - 861*Pi2)*u3) + 
//					  23616*Pi2*(-48 + 5328*u - 4*(-43976 + 9*a6)*u2 - 
//								 4*(45512 + 433*a6)*u3 + (-549824 + 1888*a6 + 3*a62)*u4) - 
//					  1024*(-32256 + 72*(24616 + 3*a6)*u - 32*(-1146688 + 729*a6)*u2 + 
//							(-39045824 - 560496*a6 + 81*a62)*u3 + 
//							12*(-9715648 + 55984*a6 + 177*a62)*u4) + 
//					  24*a5*(9216 + 288*(-3584 + 123*Pi2)*u + 
//							 2*(-9966592 + 3456*a6 + 781296*Pi2 - 15129*Pi4)*u2 + 
//							 (-16973824 + 153600*a6 + 1023360*Pi2 - 5904*a6*Pi2 - 
//							  15129*Pi4)*u3 + 16*(11794432 - 2304*a6 - 846240*Pi2 + 
//												  15129*Pi4)*u4))*eta3 - 
//				   192*u2*(-2304 + 144*(3520 + 24*a5 - 123*Pi2)*u + 
//						   (11984896 - 3456*a6 - 852144*Pi2 + 15129*Pi4 + 
//							a5*(76800 - 2952*Pi2))*u2 + (-18432*a5 + 1152*a52 + 
//														 (3776 + 24*a6 - 123*Pi2)*(-3776 + 123*Pi2))*u3)*eta4 + 
//				   9216*u4*(-144 + (-3776 + 123*Pi2)*u)*eta5))/
//	(-73728 - 96*(-123*Pi2*(-1 + u + 2*u2 + 4*u3) + 
//				  24*a5*(-1 + 2*u + 4*u2 + 8*u3 + 16*u4) + 
//				  8*(-448 + (424 + 3*a6)*u + (848 + 6*a6)*u2 + 4*(376 + 3*a6)*u3 + 
//					 24*a6*u4 + 48*a6*u5))*eta + 
//	 u*(9216 + 18432*u - 192*(3392 + 24*a5 - 123*Pi2)*u2 + 
//		(-11431936 + 4608*a6 + 834432*Pi2 - 15129*Pi4 + 
//		 24*a5*(-3008 + 123*Pi2))*u3 - 
//		(2304*a52 - 96*a5*(-3392 + 123*Pi2) + (-3776 + 123*Pi2)*
//         (-3008 + 24*a6 + 123*Pi2))*u4)*eta2 + 
//	 96*u4*(-192 + (-3776 + 123*Pi2)*u)*eta3)/
//	(-73728 - 96*(-123*Pi2*(-1 + u + 2*u2 + 4*u3) + 
//	24*a5*(-1 + 2*u + 4*u2 + 8*u3 + 16*u4) + 
//	8*(-448 + (424 + 3*a6)*u + (848 + 6*a6)*u2 + 4*(376 + 3*a6)*u3 + 
//	24*a6*u4 + 48*a6*u5))*eta + 
//	u*(9216 + 18432*u - 192*(3392 + 24*a5 - 123*Pi2)*u2 + 
//	(-11431936 + 4608*a6 + 834432*Pi2 - 15129*Pi4 + 
//	24*a5*(-3008 + 123*Pi2))*u3 - 
//	(2304*a52 - 96*a5*(-3392 + 123*Pi2) + (-3776 + 123*Pi2)*
//	(-3008 + 24*a6 + 123*Pi2))*u4)*eta2 + 
//	 96*u4*(-192 + (-3776 + 123*Pi2)*u)*eta3);
//	return dAdr;
//}
//////////////////////////////////////////////////////////////////////////


double BBHamiltonian::evaluate_eobd(std::vector<double> variables) const
{
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
	const double Dpot4PN = tuning_params_ptr->get_Dpot4PN();

	int i;
	double eobd, eta, r, u;
	Array rvec(3);
	
	eta  = wave_params_ptr->get_eta();
	
   	for (i = 0; i < 3; i++)
		rvec[i] = variables[i];
	r  = sqrt(innerproduct(rvec, rvec)); 
	u = 1./r;
	
	eobd = 1./(1.+log(1. + 6.*eta*u*u + 2.*(26. - 3.*eta)*eta*u*u*u+ (36. - Dpot4PN)*eta*eta*u*u*u*u));
	return eobd;
}


//////////////////////////////////////////////////////////////////////////
// Non-spinning D(r) Pade (1,3) potential
//
//double BBHamiltonian::evaluate_eobd(std::vector<double> variables) const
//{
//	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
//	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
//	const double Dpot4PN = tuning_params_ptr->get_Dpot4PN();
//	
//	int i;
//	double eobd, eta, r;
//	Array rvec(3);
//	
//	eta  = wave_params_ptr->get_eta();
//	
//   	for (i = 0; i < 3; i++)
//		rvec[i] = variables[i];
//	r  = sqrt(innerproduct(rvec, rvec)); 
//		
//	eobd = r*r*r/(52.*eta-6.*eta*eta+6.*eta*r+r*r*r);
//	return eobd;
//}
//////////////////////////////////////////////////////////////////////////


double BBHamiltonian::evaluate_eobd_deriv(std::vector<double> variables) const
{
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();	
	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
	const double Dpot4PN = tuning_params_ptr->get_Dpot4PN();
	int i;
	double dDdr, eta, r, u, D;
	Array rvec(3);
	
	eta  = wave_params_ptr->get_eta();
	
   	for (i = 0; i < 3; i++)
		rvec[i] = variables[i];
	r  = sqrt(innerproduct(rvec, rvec)); 
	u = 1./r;
	
	D = 1./(1.+log(1. + 6.*eta*u*u + 2.*(26. - 3.*eta)*eta*u*u*u + (36. - Dpot4PN)*eta*eta*u*u*u*u));
	//dDdr = - u*u*D*D*6*eta*u*(2. + u*(-3.*eta+26.))/(-1. + 6.*eta*eta*u*u*u -2.*eta*u*u*(3. + 26.*u));
	dDdr =  u*u*D*D*(12.*eta*u + 6.*(26. - 3.*eta)*eta*u*u + 4.*(36. - Dpot4PN)*eta*eta*u*u*u)
	/(1. + 6.*eta*u*u + 2.*(26. - 3.*eta)*eta*u*u*u + (36. - Dpot4PN)*eta*eta*u*u*u*u);
	return dDdr;
}


//////////////////////////////////////////////////////////////////////////
// Non-spinning dD/dr Pade (1,3) potential
//
//double BBHamiltonian::evaluate_eobd_deriv(std::vector<double> variables) const
//{
//	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();	
//	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
//	const double Dpot4PN = tuning_params_ptr->get_Dpot4PN();
//	int i;
//	double dDdr, eta, r;
//	Array rvec(3);
//	
//	eta  = wave_params_ptr->get_eta();
//	
//   	for (i = 0; i < 3; i++)
//		rvec[i] = variables[i];
//	r  = sqrt(innerproduct(rvec, rvec)); 
//	
//	dDdr = (6.*eta*r*r*(26.-3.*eta+2.*r))/(-6.*eta*eta+r*r*r+eta*(52.+6.*r))/(-6.*eta*eta+r*r*r+eta*(52.+6.*r));
//	
//	return dDdr;
//}
//////////////////////////////////////////////////////////////////////////


std::vector<double> BBHamiltonian::evaluate_Seffs(std::vector<double> variables) const
{
	copy_ptr<WaveformParameters> 	wave_params_ptr   = get_wave_params_ptr();

	unsigned int i;
	std::vector<double> Seffs(6);	
	Array rvec(3), pvec(3), S1vec(3), S2vec(3), SKerr_vec(3), Sstar_vec(3);
	double M = wave_params_ptr->get_M();	

	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
	
	for (i = 0; i < 3; i++)
		SKerr_vec[i]	= (S1vec[i]+S2vec[i])/(M*M);
	
	for (i = 0; i < 3; i++)
	{
		Seffs[i]	= SKerr_vec[i];
		Seffs[i + 3]	= 0.;
	}
	return Seffs;
}

/**
 * Compute sigmaKerr as defined in Eq. (5.2) of Phys.Rev. D81 (2010) 084024
 */
std::vector<double> BBHamiltonian::evaluate_sigmaKerr(std::vector<double> variables) const
{
	copy_ptr<WaveformParameters> 	wave_params_ptr   = get_wave_params_ptr();
	
	unsigned int i;
	std::vector<double> sigmaKerr(6);	
	Array S1vec(3), S2vec(3);
	double M = wave_params_ptr->get_M();	
	for (i = 0; i < 3; i++)
	{
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
//	std::cout << "S1vec: "<< " "<< S1vec << std::endl;
//	std::cout << "S2vec: "<< " "<< S2vec << std::endl;

	//BB paper Eq. 5.2
	for (i = 0; i < 3; i++)
	{
		sigmaKerr[i]	= (S1vec[i]+S2vec[i])/(M*M);
//		std::cout << "sigmaKerr: "<< " "<< sigmaKerr[i] << std::endl;
	}
	return sigmaKerr;	
}

//std::vector<double> BBHamiltonian::evaluate_sigmaStar(std::vector<double> variables) const
//{
//	copy_ptr<WaveformParameters> 	wave_params_ptr   = get_wave_params_ptr();
//	std::vector<double> 		saves 	   	  = get_saves();
//	double eta = wave_params_ptr->get_eta();	
//
//	unsigned int i;
//	std::vector<double> sigmaStar(6);	
//	Array S1vec(3), S2vec(3);
//	double m1om2 = saves[1];
//	double m2om1 = saves[2];
//	double M = wave_params_ptr->get_M(); 	
//	
//	for (i = 0; i < 3; i++)
//	{
//		S1vec[i] = variables[i + 6];
//		S2vec[i] = variables[i + 9];
//	}
//	
//	//BB paper Eq. 5.3
//	if (eta > 1.0e-16)
//	{
//		for (i = 0; i < 3; i++)
//			sigmaStar[i] = (m2om1*S1vec[i] + m1om2*S2vec[i])/(M*M);
//	}
//	else
//	{
//		for (i = 0; i < 3; i++)
//			sigmaStar[i] = 0.;
//	}
//	return sigmaStar;	
//}

/**
 * Compute sigmaStar as defined in Eq. (5.3) of Phys.Rev. D81 (2010) 084024
 */
std::vector<double> BBHamiltonian::evaluate_sigmaStar(std::vector<double> variables) const
{
	copy_ptr<Flags>    flags_ptr   = get_flags_ptr();
	copy_ptr<WaveformParameters>    wave_params_ptr   = get_wave_params_ptr();
	std::vector<double>             saves             = get_saves();
	double eta = wave_params_ptr->get_eta();
	
	unsigned int i;
	std::vector<double> sigmaStar(6);
	Array S1vec(3), S2vec(3);
//	double m1om2 = saves[1];
//	double m2om1 = saves[2];
	double M = wave_params_ptr->get_M();
	double etaF = wave_params_ptr->get_etaF();
	double m1 =0.5*(1.+sqrt(1.-4.*etaF));
	double m2 =0.5*(1.-sqrt(1.-4.*etaF));
	
	for (i = 0; i < 3; i++)
	{
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
	
	
	if (eta > 1.0e-16 || flags_ptr->get_flag_ParticleSpin() == 1)		
	{
		for (i = 0; i < 3; i++)
			//                      sigmaStar[i] = (m2om1*S1vec[i] + m1om2*S2vec[i])/(M*M);
			sigmaStar[i] = (m2/m1*S1vec[i] + m1/m2*S2vec[i])/(M*M);
	}
	else
	{
		for (i = 0; i < 3; i++)
			sigmaStar[i] = 0.;
	}
//	std::cout << "BBHamiltonian::evaluate_sigmaStar" << std::endl;
//	std::cout << "S1vec = " << S1vec << std::endl;
//	std::cout << "S2vec = " << S2vec << std::endl;
////	abort();
	return sigmaStar;
}

double BBHamiltonian::evaluate_Heff(std::vector<double> variables) const 
{	
	
	copy_ptr<Flags>			flags_ptr	  = get_flags_ptr();
	copy_ptr<WaveformParameters> 	wave_params_ptr   = get_wave_params_ptr();
	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 		saves 	   	  = get_saves();
	
	Array rvec(3), pvec(3), S1vec(3), S2vec(3), nvec(3);
	std::vector<double> sigmaKerr(3), sigmaStar(3);
	int i;
	double eta = wave_params_ptr->get_eta();
	double a;
	double wfd2, wfd1, d1;
	
	double r, px, py, pz, nx, ny, nz, Skerr_x, Skerr_y, Skerr_z, costheta, xi2, xi_x, xi_y, xi_z, vx, vy, vz, pxir, pvr, pn, pr, pf, ptheta2,
	w2, rho2, u, deltaT, deltaR, Lambda, D, qq, ww, B, w, MU, nu, BR, wr, nur, mur, wcos, nucos, mucos, ww_r, Lambda_r,
    deltaU, deltaU_u, Q, deltaT_r, pn2, pp, deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z, sx, sy, sz, sxi, sv, sn, s3,
	H, Hns, Hs, Hss, Hwcos, Hwr, HSOL, HSONL, e3_x, e3_y, e3_z;
	
	double aSK2, r2, Deltat, eobd, eoba, xi, ndotph;
	
    const double aa=tuning_params_ptr->get_aeta(), bb=tuning_params_ptr->get_beta(), dheff_SS=tuning_params_ptr->get_dheff_SS(); //spin gauge parameters
	
   	for (i = 0; i < 3; i++)
	{
		rvec[i] = variables[i];
		pvec[i] = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
//	std::cout << "Heff variables" << std::endl;
//	for(int i=0; i<12; i++)
//		std::cout << variables[i] << std::endl;
	
	sigmaKerr = evaluate_sigmaKerr(variables);
	sigmaStar = evaluate_sigmaStar(variables);
	a = sqrt(sigmaKerr[0]*sigmaKerr[0]+sigmaKerr[1]*sigmaKerr[1]+sigmaKerr[2]*sigmaKerr[2]);
	aSK2 = a * a;
	
	r  = sqrt(innerproduct(rvec, rvec));
	r2 = r * r;
	nvec = rvec / r;

////////////////////////////////////////////////////////////////
//	double prold;
//	prold = innerproduct(nvec,pvec);
//	if (flags_ptr->get_flag_tortoiseR() == 1)
//	{
//		ndotph = innerproduct(nvec, pvec);
//		eoba = evaluate_eoba(variables);
//		eobd = evaluate_eobd(variables);
//		Deltat	= r2 * eoba;
////		std::cout << "Deltat HAM " << Deltat <<std::endl;
////		std::cout << "eoba HAM " << eoba <<std::endl;
////		std::cout << "eobd HAM " << eobd <<std::endl;
////		std::cout << "a HAM " << a <<std::endl;
//		//Tortoise coord tranformation	
//		xi = Deltat / sqrt(eobd) / (r2 + aSK2);		
////		std::cout << "xi in Ham= "<< xi <<std::endl;
////		std::cout << "ndotp Ham " << ndotph <<std::endl;
//		for (i = 0; i < 3; i++)
//			pvec[i] = pvec[i] - rvec[i]/r*(xi-1.)/xi*ndotph;
//	}
//	if(flags_ptr->get_flag_pr4tortoise() == 1)
//		prold = innerproduct(nvec,pvec);
////////////////////////////////////////////////////////////////
//	double prold;
//	prold = innerproduct(nvec,pvec);
//	//	if (flags_ptr->get_flag_tortoiseR() == 1)
//	//	{
//	ndotph = innerproduct(nvec, pvec);
//	eoba = evaluate_eoba(variables);
//	eobd = evaluate_eobd(variables);
//	Deltat	= r2 * eoba;
//	//		std::cout << "Deltat HAM " << Deltat <<std::endl;
//	//		std::cout << "eoba HAM " << eoba <<std::endl;
//	//		std::cout << "eobd HAM " << eobd <<std::endl;
//	//		std::cout << "a HAM " << a <<std::endl;
//	//Tortoise coord tranformation	
//	xi = Deltat / sqrt(eobd) / (r2 + aSK2);		
//	//		std::cout << "xi in Ham= "<< xi <<std::endl;
//	//		std::cout << "ndotp Ham " << ndotph <<std::endl;
//	if (flags_ptr->get_flag_tortoiseR() == 1)
//	{	
//		for (i = 0; i < 3; i++)
//			pvec[i] = pvec[i] - rvec[i]/r*(xi-1.)/xi*ndotph;
//	}
//	else
//	{
//		prold = 1.*prold;
//	}
//	if(flags_ptr->get_flag_pr4tortoise() == 1)
//		prold = innerproduct(nvec,pvec);
////////////////////////////////////////////////////////////////
//	std::cout << "flags_ptr->get_flag_pr4tortoise() "<< flags_ptr->get_flag_pr4tortoise() << std::endl;
//	std::cout << "flags_ptr->get_flag_tortoiseR() "<< flags_ptr->get_flag_tortoiseR() << std::endl;
	double prold;
	if (flags_ptr->get_flag_pr4tortoise() == 1)
	{
		prold = innerproduct(nvec,pvec);
	}
	else if (flags_ptr->get_flag_pr4tortoise() == 0)
	{
		ndotph = innerproduct(nvec, pvec);
		eoba = evaluate_eoba(variables);
		eobd = evaluate_eobd(variables);
		Deltat	= r2 * eoba;
//		std::cout << "Deltat HAM " << Deltat <<std::endl;
//		std::cout << "eoba HAM " << eoba <<std::endl;
//		std::cout << "eobd HAM " << eobd <<std::endl;
//		std::cout << "a HAM " << a <<std::endl;
		//Tortoise coord tranformation	
		xi = Deltat / sqrt(eobd) / (r2 + aSK2);		
//		std::cout << "xi in Ham= "<< xi <<std::endl;
		//		std::cout << "ndotp Ham " << ndotph <<std::endl;
		if (flags_ptr->get_flag_tortoiseR() == 1)
		{		
			prold = innerproduct(nvec,pvec);
			for (i = 0; i < 3; i++)
				pvec[i] = pvec[i] - rvec[i]/r*(xi-1.)/xi*ndotph;
		}
		else if (flags_ptr->get_flag_tortoiseR() == 0)
		{
			prold = xi*innerproduct(nvec,pvec);
		}
	}
	else
	{
		std::cout << "In BBHamiltonian::evaluate_Heff inconsistent flag_tortoiseR" << std::endl;
		abort();
	}
	
//	std::cout << "xi " << xi << std::endl;
//	std::cout << "prold " << prold << std::endl;
//	
//	px = pvec[0]/eta;
//	py = pvec[1]/eta;
//	pz = pvec[2]/eta;
	px = pvec[0];
	py = pvec[1];
	pz = pvec[2];
	
//	std::cout << "px = " << px << std::endl;
//	std::cout << "py = " << py << std::endl;
//	std::cout << "pz = " << pz << std::endl;
//	abort();	
    nx = nvec[0];
    ny = nvec[1];
    nz = nvec[2];   	
	
	Array Lhat(3);
	if(innerproduct(pvec,pvec)!=0.){
		Lhat = crossproduct(rvec,pvec);
		Lhat = Lhat/sqrt(innerproduct(Lhat,Lhat));
	}
	else {
		Lhat[0]=0.;
		Lhat[1]=0.;
		Lhat[2]=0.;
	}
    //Tuning saves
//	KK0  = saves[0];
//	Ka   = saves[4];
	wfd1 = saves[5];
	wfd2 = saves[6];
	d1   = saves[7];
	

//	std::cout << "a " << a <<std::endl;
//	std::cout << "eta " << eta <<std::endl;
//	std::cout << "KK " << KK <<std::endl;
//	std::cout << "k0 " << k0 <<std::endl;
//	std::cout << "k1 " << k1 <<std::endl;
//	std::cout << "k2 " << k2 <<std::endl;
//	std::cout << "k3 " << k3 <<std::endl;
//	std::cout << "k4 " << k4 <<std::endl;
//	std::cout << " " <<std::endl;
//	abort();

	//BB paper Eq. 5.67 so that Skerr = sigmaKerr
    Skerr_x = sigmaKerr[0];
    Skerr_y = sigmaKerr[1];
    Skerr_z = sigmaKerr[2];

//	std::cout << "nhat = " << nx << " " << ny << " " << nz << std::endl;
	double S1Norm=innerproduct(S1vec,S1vec);
        double S2Norm=innerproduct(S2vec,S2vec);
	double LhatNorm=innerproduct(Lhat,Lhat);
//        std::cout << "S1Norm " << S1Norm << std::endl;
//        std::cout << "S2Norm " << S2Norm << std::endl;
//        std::cout << "LhatNorm " << LhatNorm << std::endl;
//	std::cout << "Lhat " << Lhat << std::endl;

	if(S1Norm==0. && S2Norm==0. && LhatNorm==0.) {
		e3_x = 0.;
		e3_y = 1.;
		e3_z = 0.;
	}
        if(S1Norm==0. && S2Norm==0. && LhatNorm!=0.) {
                e3_x = Lhat[0];
                e3_y = Lhat[1];
                e3_z = Lhat[2];
        }
        if((S1Norm!=0. || S2Norm!=0.) && a!=0.) {
		e3_x = Skerr_x/a;
		e3_y = Skerr_y/a;
		e3_z = Skerr_z/a;
	}
        if((S1Norm!=0. || S2Norm!=0.) && a==0. && LhatNorm!=0) {
                e3_x = Lhat[0];
                e3_y = Lhat[1];
                e3_z = Lhat[2];
        }
        if((S1Norm!=0. || S2Norm!=0.) && a==0. && LhatNorm==0) {
                e3_x = 0.;
                e3_y = 1.;
                e3_z = 0.;
	}
	
    costheta = e3_x*nx+e3_y*ny+e3_z*nz; 
	
    xi2 = 1.-costheta*costheta; 
//	std::cout << "xi2 " << xi2 << std::endl;
	
    xi_x =  -e3_z*ny + e3_y*nz;
    xi_y =   e3_z*nx - e3_x*nz;
    xi_z =  -e3_y*nx + e3_x*ny;
	
    vx =  -nz*xi_y + ny*xi_z;
    vy =   nz*xi_x - nx*xi_z;
    vz =  -ny*xi_x + nx*xi_y;
	
    pxir = (px*xi_x+py*xi_y+pz*xi_z)*r;
    pvr  = (px*vx+py*vy+pz*vz)*r;
    pn   = px*nx+py*ny+pz*nz;
	
//	std::cout << "flags_ptr->get_flag_pr4tortoise()" << flags_ptr->get_flag_pr4tortoise() << std::endl;
//	std::cout << "flags_ptr->get_flag_tortoiseR()" << flags_ptr->get_flag_tortoiseR() << std::endl;
//	std::cout << px << std::endl;
//	std::cout << py << std::endl;
//	std::cout << pz << std::endl;
//	std::cout << xi_x << std::endl;
//	std::cout << xi_y << std::endl;
//	std::cout << xi_z << std::endl;
	
    pr      = pn;
    pf      = pxir;
    ptheta2 = pvr*pvr/xi2;
	
//	std::cout << "pr "<< pr << std::endl;
//	std::cout << "pf "<< pf << std::endl;
//	std::cout << "ptheta2 "<< ptheta2 << std::endl;
//	std::cout <<  "a " << a << std::endl; 
//	abort();	
    w2   = r*r + a*a;
    rho2 = r*r + a*a*costheta*costheta;
    u    = 1./r;
	
///////////////////////////////////////////////////////////////////////////////////////////////	
//  Original Enrico's implementdation of Delta_u
//  double bulk, LogTerms, KK, k0, k1, k2, k3, k4;
////
//////K(eta) polynomial
////////	c0  = saves[4];
////////	c1  = -8.*(c0 - KK0 + Ka*a*a);//c1 or K_1 in BB paper Eq. 6.11 
////////	c2  = 16.*(c0 - KK0 + Ka*a*a);//c2 or K_2 in BB paper Eq. 6.11 
////////	KK  = c0 + c1*eta + c2*eta*eta;//KK or K in BB paper Eq. 6.9
//	KK = saves[3] + saves[4]*a*a;
////	
////	//\Delta_i's
//	k0  = KK*(-2. + eta*KK);//k0 or \Delta_1 in BB paper Eq. 5.77
//	k1  = -2.*(k0 + KK)*(-1. + eta*KK);//k1 or \Delta_2 in BB paper Eq. 5.78
//	k2  = (k1*(4. + k1 - 4.*eta*KK))/2. - intPower(a,2)*k0*intPower(-1. + eta*KK,2);
//	k3  = -intPower(k1,3)/3. + k1*k2 + intPower(k1,2)*(-1. + eta*KK) - 2.*(1. + k2 - eta*KK)*(-1. + eta*KK) - intPower(a,2)*k1*intPower(-1. + eta*KK,2);
//	k4  = (24.*intPower(k1,4) - 96.*intPower(k1,2)*k2 + 48.*intPower(k2,2) - 64.*intPower(k1,3)*(-1. + eta*KK) 
//	   + 48.*intPower(a,2)*(intPower(k1,2) - 2.*k2)*intPower(-1. + eta*KK,2) + 
//	   96.*k1*(k3 + 2.*k2*(-1. + eta*KK)) - (-1. + eta*KK)*(192.*k3 + (-1. + eta*KK)*(-3008. + 123.*intPower(Pi,2))))/96.;  
//	
//	bulk = intPower(-1. + eta*KK,-2) + (2.*u)/(-1. + eta*KK) + intPower(a,2)*intPower(u,2);	
//    LogTerms = ( 1. + eta*k0 + eta*log(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4)) );
//	deltaU = bulk*LogTerms;
/////////////////////////////////////////////////////////////////////////////////////////////		
	deltaU = evaluate_eoba(variables);
    deltaT = r*r*deltaU;
//    std::cout << "bulk = " << bulk << std::endl;
//    std::cout << "LogTerms = " << LogTerms << std::endl;
//	std::cout << "DE = " << deltaT << std::endl;
//	std::cout << "deltaU = " << deltaU << std::endl;
//	std::cout << "eoba = " << evaluate_eoba(variables) << std::endl;

//    deltaU_u = 2.*(1./(-1. + eta*KK) + intPower(a,2)*u)*LogTerms + 
//	bulk* (eta*(k1 + u*(2.*k2 + u*(3.*k3 + 4.*k4*u))))/(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4));
	deltaU_u = -r*r*evaluate_eoba_deriv(variables);
//	deltaU_u = evaluate_eoba_deriv(variables);

//	if(r<10){
//	std::cout << "deltaU_u new = " << deltaU_u << std::endl;
//		std::cout << "deltaU_u Andrea = " <<  2.*(1./(-1. + eta*KK) + intPower(a,2)*u)*LogTerms + 
//			bulk* (eta*(k1 + u*(2.*k2 + u*(3.*k3 + 4.*k4*u))))/(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4))<< std::endl;
//	std::cout << " " << std::endl;
//	}

    deltaT_r = 2.*r*deltaU-deltaU_u;
	
    Lambda = w2*w2 - a*a*deltaT*xi2;
//	std::cout << "Lambda "<< Lambda << std::endl;
//	abort();
	
///////////////////////////////////////////////////////////////////////////////////////////////	
//  Original Enrico's implementdation of 1/eobd
//  D = 1.+log(1. + 6.*eta*u*u + 2.*(26. - 3.*eta)*eta*u*u*u + (36. - Dpot4PN)*eta*eta*u*u*u*u);
///////////////////////////////////////////////////////////////////////////////////////////////	
	D = 1./evaluate_eobd(variables);
//	std::cout << "D "<< D << std::endl;
//	abort();
	
    deltaR = deltaT*D;
	
    qq = 2.*eta*(4. - 3.*eta);
	
    ww = 2.*a*r + wfd2*eta*a*a*a*u + wfd1*eta*a*u;
//	std::cout << "ww " << a << " " << eta <<" " <<  wfd1 << " " << wfd2 << " " << ww << std::endl;
//	std::cout << "prold " << prold << std::endl;
////////////////////////////////////////////////////////////////////	
	//I remove the 1 under the sqrt for the LR
//    Hns=sqrt(1. - flags_ptr -> get_flag_LR() + intPower(pr,4)*qq/intPower(r,2) 
//			 + ptheta2/rho2 + intPower(pf,2)*rho2/(Lambda*xi2) + intPower(pr,2)*deltaR/rho2)/
//	sqrt(Lambda/(rho2*deltaT)) + pf*ww/Lambda;
	
//	Mar 30, 2011: introduced prold = prstar to avoid the D^2/A^4 term when doing the
//	inverse tortoise transformation
    Hns=sqrt(1. - flags_ptr -> get_flag_LR() + intPower(prold,4)*qq/intPower(r,2) 
			 + ptheta2/rho2 + intPower(pf,2)*rho2/(Lambda*xi2) + intPower(pr,2)*deltaR/rho2)/
	sqrt(Lambda/(rho2*deltaT)) + pf*ww/Lambda;
//////////////////////////////////////////////////////////////////	
//	std::cout << "pr " << pr << std::endl;
//	std::cout << "prold " << prold << std::endl;
//	std::cout << "qq " << qq << std::endl;
//	std::cout << "rho2 " << rho2 << std::endl;
//	std::cout << "pf " << pf << std::endl;
//	std::cout << "Lambda " << Lambda << std::endl;
//	std::cout << "xi2 " << xi2 << std::endl;
//	std::cout << "deltaR " << deltaR << std::endl;
//	std::cout << "intPower(prold,4)*qq/intPower(r,2) " << intPower(prold,4)*qq/intPower(r,2) << std::endl;
//	std::cout << "intPower(pf,2)*rho2/(Lambda*xi2) " << intPower(pf,2)*rho2/(Lambda*xi2) << std::endl;
//	std::cout << "intPower(pr,2)*deltaR/rho2 " << intPower(pr,2)*deltaR/rho2 << std::endl;
//	std::cout << "sqrt(Lambda/(rho2*deltaT)) " << sqrt(Lambda/(rho2*deltaT)) << std::endl;
    B = sqrt(deltaT);
    w = ww/Lambda;
    nu = 0.5*log(deltaT*rho2/Lambda);
    MU = 0.5*log(rho2);  
	
    Lambda_r=4.*r*w2 - a*a*deltaT_r*xi2;
	
    ww_r=2.*a - (intPower(a,3)*wfd2*eta)/intPower(r,2) - wfd1*eta*a/(r*r);
	
    BR = (-2.*deltaT + sqrt(deltaR)*deltaT_r)/(2.*sqrt(deltaR*deltaT));
    wr = (-Lambda_r*ww + Lambda*ww_r)/(Lambda*Lambda);
    nur = (r/rho2 + (w2 * (-4.*r*deltaT + w2*deltaT_r) ) / (2.*deltaT*Lambda) );
    mur = (r/rho2 - 1./sqrt(deltaR));
	
    wcos = -2.*a*a*costheta*deltaT*ww/(Lambda*Lambda);  
    nucos= a*a*costheta*w2*(w2-deltaT)/(rho2*Lambda);  
    mucos=a*a*costheta/rho2;
	
    Q=1. + intPower(pvr,2)/(exp(2.*MU)*xi2) + exp(2.*nu)*intPower(pxir,2)/(intPower(B,2)*xi2)
	+ intPower(pn,2)*deltaR/exp(2.*MU);
	
    pn2=pr*pr*deltaR/rho2;
    pp=Q-1.;
	
	deltaSigmaStar_x = 0.;
	deltaSigmaStar_y = 0.;
	deltaSigmaStar_z = 0.;
	
	switch (flags_ptr->get_flag_DeltaSigmaStar_25PN())
	{
		case 1:
		{
			deltaSigmaStar_x +=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sigmaKerr[0] - 8.*bb*(1. + 3.*pn2*r - pp*r)*sigmaStar[0] + 
							  eta*(-8.*sigmaKerr[0] - 36.*pn2*r*sigmaKerr[0] + 3.*pp*r*sigmaKerr[0] + 14.*sigmaStar[0] 
								   - 30.*pn2*r*sigmaStar[0] + 4.*pp*r*sigmaStar[0]))/(12.*r) + d1*eta*sigmaStar[0]/r/r/r;
			
			deltaSigmaStar_y +=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sigmaKerr[1] - 8.*bb*(1. + 3.*pn2*r - pp*r)*sigmaStar[1] + 
							  eta*(-8.*sigmaKerr[1] - 36.*pn2*r*sigmaKerr[1] + 3.*pp*r*sigmaKerr[1] + 14.*sigmaStar[1] 
								   - 30.*pn2*r*sigmaStar[1] + 4.*pp*r*sigmaStar[1]))/(12.*r) + d1*eta*sigmaStar[1]/r/r/r;
			
			deltaSigmaStar_z +=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sigmaKerr[2] - 8.*bb*(1. + 3.*pn2*r - pp*r)*sigmaStar[2] + 
							  eta*(-8.*sigmaKerr[2] - 36.*pn2*r*sigmaKerr[2] + 3.*pp*r*sigmaKerr[2] + 14.*sigmaStar[2] 
								   - 30.*pn2*r*sigmaStar[2] + 4.*pp*r*sigmaStar[2]))/(12.*r) + d1*eta*sigmaStar[2]/r/r/r;
			break;
		}
		case 0:
			break;
		default:
		{
			std::cout << "flag_DeltaSigmaStar_15PN must be 0 or 1" << std::endl;
			abort();
			break;
		}
	}

	switch (flags_ptr->get_flag_DeltaSigmaStar_35PN())
	{
		case 1:
		{
			double u1, u2;
			double a1 = 0.;
			double a2 = 0.;
			double a3 = 0.;
			double b1 = 0.;
			double b2 = 0.;
			double b3 = 0.;
//			Wrong mapping
//			u1=(-48*b2 + 706*eta - 54*(eta*eta) + (-48*b1 + 48*b2 + eta*(-368 +
//				39*eta))*pp*r - 120*(2*b3 - 3*(eta*eta))*(pn2*pn2)*(r*r) +
//				(48*b1 + eta*(-23 + 78*eta))*(pp*pp)*(r*r) - 6*pn2*r*(16*b1 + 32*b2 +
//				24*b3 - 47*eta + 54*(eta*eta) +	(24*b1 - 24*b3 - eta*(16 + 33*eta))*pp*r))/(72.*(r*r));
//
//			u2=-(16*(6*a2 + 7*eta*(8 + 3*eta)) + 4*(24*a1 - 24*a2 + eta*(217 +
//				3*eta))*pp*r + 30*(16*a3 - 27*(eta*eta))*(pn2*pn2)*(r*r)
//				 - 3*(32*a1 + 3*eta*(-5 + 24*eta))*(pp*pp)*(r*r) + 6*pn2*r*(32*a1 +	64*a2 + 48*a3 + 16*eta + 147*(eta*eta) + 
//			    3*(16*a1 - 16*a3- eta*(2 + 35*eta))*pp*r))/(144.*(r*r));

//			New mapping
//			u1 =  -(1/(72*r*r))*(48*b2 - 706*eta + 54*eta*eta +
//					2*(24*b1 - 24*b2 + (103 - 60*eta)*eta)*pp*r +
//					120*(2*b3 - 3*eta*eta)*pn2*pn2*r*r + (-48*b1 + eta*(23 + 3*eta))*
//					pp*pp*r*r + 6*pn2*r*(16*b1 + 32*b2 + 24*b3 - 47*eta +
//					54*eta*eta + (24*b1 - 24*b3 + eta*(-16 + 21*eta))*pp*r));
//			
//			u2 = 1/(144*r*r)*(-16*(6*a2 + 7*eta*(8 + 3*eta)) +
//					4*(-24*a1 + 24*a2 + eta*(-109 + 51*eta))*pp*r +
//					30*(-16*a3 + 27*eta*eta)*pn2*pn2*r*r +
//					3*(32*a1 - 15*eta)*pp*pp*r*r -
//					6*pn2*r*(32*a1 + 64*a2 + 48*a3 + 16*eta + 147*eta*eta +
//					3*(16*a1 - 16*a3 + eta*(-2 + 13*eta))*pp*r));

//			New mapping	with gauge params
			double aaa = 0.;
			double bbb = 0.;
			
//			double aaa = -3/2*eta;
//			double bbb = -5/4*eta;
//			b1 = 1/16.*eta*(19 + 10*eta);
//			b2 = -(1/16.)*eta*(-39 + 20*eta);
//			b3 = 3./2*eta*eta;
//			a1 = 1/8.*eta*(6 + 7*eta);
//			a2 = -(1/4.)*eta*(-5 + 7*eta);
//			a3 = 27./16*eta*eta;
			
//			u1 = -(16*bbb + 48*b2 - 706*eta + 96*bbb*eta + 54*eta*eta +
//				   2*(24*b1 - 24*b2 + bbb*(2 - 48*eta) + 103*eta - 60*eta*eta)*pp*
//				   r + 120*(2*b3 - 3*eta*eta)*pn2*pn2*r*
//				   r + (-20*bbb - 48*b1 + eta*(23 + 3*eta))*pp*pp*r*r +
//				   6*pn2*r*(16*b1 + 32*b2 + 24*b3 - 47*eta + 54*eta*eta +
//							8*bbb*(5 + 6*eta) + (10*bbb + 24*b1 - 24*b3 - 16*eta +
//												21*eta*eta)*pp*r))/(72*r*r);
//			
//			u2 = (-16*(6*a2 + 7*eta*(8 + 3*eta) + 2*aaa*(1 + 6*eta)) +
//				  4*(-24*a1 + 24*a2 - 109*eta + 51*eta*eta + aaa*(-2 + 48*eta))*pp*
//				  r - 30*(16*a3 - 27*eta*eta)*pn2*pn2*r*
//				  r + (40*aaa + 96*a1 - 45*eta)*pp*pp*r*r -
//				  6*pn2*r*(32*a1 + 64*a2 + 48*a3 + 16*eta + 147*eta*eta +
//						   16*aaa*(5 + 6*eta) + (20*aaa + 48*a1 - 48*a3 - 6*eta +
//												39*eta*eta)*pp*r))/(144*r*r);

//			Final mapping with gauge params
			aaa = -3/2*eta;
			bbb = -5/4*eta;
			a1 = eta*eta/2;
			a2 = -(1/8)*eta*(-7 + 8*eta);
			a3 = -((9*eta*eta)/16);
			b1 = 1/16*eta*(9 + 5*eta);
			b2 = -(1/8)*eta*(-17 + 5*eta);
			b3 = -3/8*eta*eta;
		        
                        aaa = 0.;
                        bbb = 0.;
                        a1 = 0.;
                        a2 = 0.;
                        a3 = 0.;
                        b1 = 0.;
                        b2 = 0.;
                        b3 = 0.;
	
			u1 =-(2*(24*b2 + eta*(-353 + 27*eta) + bbb*(56 + 60*eta)) + 
				  2*(24*b1 - 24*b2 + bbb*(14 - 66*eta) + 103*eta - 60*intPower(eta,2))*pp*
				  r + 120*(2*b3 - 3*eta*(bbb + eta))*intPower(pn2,2)*intPower(r,2) + 
				  (-48*b1 + 4*bbb*(1 + 3*eta) + eta*(23 + 3*eta))*intPower(pp,2)*
				  intPower(r,2) + 6*pn2*r*(16*b1 + 32*b2 + 24*b3 - 47*eta + 
										54*intPower(eta,2) + 24*bbb*(1 + eta) + 
										(24*b1 - 24*b3 - 16*eta + 21*intPower(eta,2) + bbb*(-2 + 30*eta))*pp*
										r))/(72.*intPower(r,2));			
			
			u2 = (-16*(6*a2 + 7*eta*(8 + 3*eta) + aaa*(14 + 15*eta)) + 
				  4*(-24*a1 + 24*a2 - 109*eta + 51*intPower(eta,2) + 2*aaa*(-7 + 33*eta))*
				  pp*r + 30*(-16*a3 + 3*eta*(8*aaa + 9*eta))*intPower(pn2,2)*intPower(r,2) + 
				  (96*a1 - 45*eta - 8*aaa*(1 + 3*eta))*intPower(pp,2)*intPower(r,2) - 
				  6*pn2*r*(32*a1 + 64*a2 + 48*a3 + 16*eta + 147*intPower(eta,2) + 
						   48*aaa*(1 + eta) + (48*a1 - 48*a3 - 6*eta + 39*intPower(eta,2) + 
											 aaa*(-4 + 60*eta))*pp*r))/(144.*intPower(r,2));
			
			deltaSigmaStar_x += u1*sigmaStar[0] + u2*sigmaKerr[0];
			deltaSigmaStar_y += u1*sigmaStar[1] + u2*sigmaKerr[1];
			deltaSigmaStar_z += u1*sigmaStar[2] + u2*sigmaKerr[2];
			break;
		}
		case 0:
			break;
		default:
		{
			std::cout << "flag_DeltaSigmaStar_35PN must be 0 or 1" << std::endl;
			abort();
			break;
		}
	}
	
//	double etaF=wave_params_ptr->get_etaF();
//	//I set Sstar to 0 for the LR
//    sx= (1. - flags_ptr->get_flag_LR())*(sigmaStar[0]+deltaSigmaStar_x)/etaF;
//    sy= (1. - flags_ptr->get_flag_LR())*(sigmaStar[1]+deltaSigmaStar_y)/etaF;
//    sz= (1. - flags_ptr->get_flag_LR())*(sigmaStar[2]+deltaSigmaStar_z)/etaF;     
	//I set Sstar to 0 for the LR
    sx= (1. - flags_ptr->get_flag_LR())*(sigmaStar[0]+deltaSigmaStar_x);
    sy= (1. - flags_ptr->get_flag_LR())*(sigmaStar[1]+deltaSigmaStar_y);
    sz= (1. - flags_ptr->get_flag_LR())*(sigmaStar[2]+deltaSigmaStar_z);
   	
    sxi=sx*xi_x+sy*xi_y+sz*xi_z;
    sv=sx*vx+sy*vy+sz*vz;
    sn=sx*nx+sy*ny+sz*nz; 
	
    s3=sx*e3_x+sy*e3_y+sz*e3_z;  
	
    Hwr=(exp(-3.*MU - nu)*sqrt(deltaR)*(exp(2.*(MU + nu))*intPower(pxir,2)*sv - B*exp(MU + nu)*pvr*pxir*sxi + 
										intPower(B,2)*xi2*(exp(2.*MU)*(sqrt(Q) + Q)*sv + pn*pvr*sn*sqrt(deltaR) - intPower(pn,2)*sv*deltaR)))/(2.*B*(1. + sqrt(Q))*sqrt(Q)*xi2);
	
    Hwcos=(exp(-3.*MU - nu)*(sn*(-(exp(2.*(MU + nu))*intPower(pxir,2)) + intPower(B,2)*(intPower(pvr,2) - exp(2.*MU)*(sqrt(Q) + Q)*xi2)) - 
							 B*pn*(B*pvr*sv - exp(MU + nu)*pxir*sxi)*sqrt(deltaR)))/(2.*B*(1. + sqrt(Q))*sqrt(Q));
	
    HSOL=(exp(-MU + 2.*nu)*(-B + exp(MU + nu))*pxir*s3)/(intPower(B,2)*sqrt(Q)*xi2);
	
    HSONL=(exp(-2.*MU + nu)*(-(B*exp(MU + nu)*nucos*pxir*(1. + 2.*sqrt(Q))*sn*xi2) + 
							 (-(BR*exp(MU + nu)*pxir*(1. + sqrt(Q))*sv) + B*(exp(MU + nu)*nur*pxir*(1. + 2.*sqrt(Q))*sv + B*mur*pvr*sxi + 
																			 B*sxi*(-(mucos*pn*xi2) + sqrt(Q)*(mur*pvr - nur*pvr + (-mucos + nucos)*pn*xi2))))*sqrt(deltaR)))/(intPower(B,2)*(sqrt(Q) + Q)*xi2);   
	
    Hs= w*s3 + Hwr*wr + Hwcos*wcos + HSOL + HSONL;
	
    Hss=-0.5/(r*r*r)*(sx*sx+sy*sy+sz*sz-3.*sn*sn);
//	std::cout <<"sx " <<  sx << std::endl;
//	std::cout <<"sy " <<  sy << std::endl;
//	std::cout <<"sz " <<  sz << std::endl;
//	std::cout <<"sn " <<  sn << std::endl;


    H=Hns+Hs+Hss + (1. - flags_ptr->get_flag_LR())*dheff_SS * eta * (sigmaStar[0]*sigmaKerr[0]+sigmaStar[1]*sigmaKerr[1]+sigmaStar[2]*sigmaKerr[2]) / (r*r*r*r);
	
//	for(int i=0; i<12; i++)
//	{
//	    std::cout << variables[i] << std::endl;
//	}
//	std::cout <<"BBH " <<  H << std::endl;
//	std::cout <<"eta " <<  eta << std::endl;
//	std::cout <<"r " <<  r << std::endl;
//	std::cout <<"Hns " <<  Hns << std::endl;
//	std::cout <<"Hs " <<  Hs << std::endl;
//	std::cout <<"Hss " <<  Hss << std::endl;
//	std::cout <<"dSS term " <<   (1. - flags_ptr->get_flag_LR())*dheff_SS * eta * (sigmaStar[0]*sigmaKerr[0]+sigmaStar[1]*sigmaKerr[1]+sigmaStar[2]*sigmaKerr[2]) / (r*r*r*r) << std::endl;
//	abort();
	
//	if(r<15.)
//		std::cout << r<< " " << H << " " << deltaU << " " << 1./D << std::endl;

    /*if(H!=H)
    {
  	    std::cout << "Heff is a nan! H" << std::endl;
	    std::cout << "Heff variables" << std::endl;
	    for(int i=0; i<12; i++)
	    {
	        std::cout << variables[i] << std::endl;
	    }
	    std::cout <<"H " <<  H << std::endl;
	    std::cout <<"eta " <<  eta << std::endl;
	    std::cout <<"r " <<  r << std::endl;
	    std::cout <<"Hns " <<  Hns << std::endl;
	    std::cout <<"Hs " <<  Hs << std::endl;
	    std::cout <<"Hss " <<  Hss << std::endl;
	    std::cout << "S1Norm " << S1Norm << std::endl;
	    std::cout << "S2Norm " << S2Norm << std::endl;
        std::cout << "LhatNorm " << LhatNorm << std::endl;
	    for(int i=0; i<12; i++)
            std::cout << variables[i] << std::endl;
		abort();
    }*/
	return H;
	
}


double BBHamiltonian::evaluate(std::vector<double> variables) const
{
	extern int hamilton_count;
	extern time_t hamilton_time;
	time_t track_time;
	hamilton_count += 1;
	track_time = time(NULL);
	copy_ptr<Flags>			flags_ptr	  = get_flags_ptr();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	copy_ptr<TuningParameters> tuning_params_ptr = get_tuning_params_ptr();
	
	double eta, Heff,Hreal;
	//std::cout << variables[0] << "   " << variables[5] << std::endl;
	eta	= wave_params_ptr->get_eta();
	Heff 	= evaluate_Heff(variables);

	hamilton_time += time(NULL) - track_time;
	//For the LR we have eta=0
	if(eta > 1.0e-16 && flags_ptr->get_flag_LR()==0)
	{
		Hreal=sqrt(1 + 2 * eta * (Heff - 1.) )/eta;
	}
	else
	{
		Hreal=Heff;
	}
				   
	return Hreal;
}


double BBHamiltonian::nonKeplerianRadius() const
{
	SphBBHamiltonian SphspinH(*this);
	return SphspinH.nonKeplerianRadius();
}

double BBHamiltonian::dHdr_r(double r) const
{
	std::cout << "Error! BBHamiltonian::dHdr_r: ";
	std::cout << "member function not implemented for BBHamiltonian." << std::endl;

	abort();
	return 0.0;
}

double BBHamiltonian::dHdr_pp(double pp) const
{
	std::cout << "Error! BBHamiltonian::dHdr_pp: ";
	std::cout << "member function not implemented for BBHamiltonian."<< std::endl;
	abort();
	return 0.0;
}

double BBHamiltonian::LightRingRadius() const
{
	SphBBHamiltonian SphspinH(*this);
	return SphspinH.LightRingRadius();
}

/**
 * Compute the direction orthogonal to the instantaneous orb plane for the dynamical configuration prvided in the argument
 */
Array BBHamiltonian::evaluate_LNhat(std::vector<double> variables /**< Array of the dynamical configuration*/) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	
	unsigned int i;
	Array LNhat(3);	
	Array rvec(3), pvec(3), S1(3), S2(3), rdotvec(3);
	
	for (i = 0; i < 3; i++)
	{
		rvec[i]	 = variables[i];
		pvec[i]  = variables[i + 3];
		S1[i] = variables[i + 6];
		S2[i] = variables[i + 9];
	}
	
//	std::cout << "Vars in LNhat" << std::endl;
//	for(i=0; i<12; i++)
//		std::cout << variables[i] << std::endl;
//	std::cout << " " << std::endl;

	//Compute rdot vector
	rdotvec[0] = temphptr->evaluate_deriv(3);
	rdotvec[1] = temphptr->evaluate_deriv(4);
	rdotvec[2] = temphptr->evaluate_deriv(5);
//	std::cout << "rdotvec" << std::endl;
//	std::cout << rdotvec[0] << std::endl;
//	std::cout << rdotvec[1] << std::endl;
//	std::cout << rdotvec[2] << std::endl;

	//Build LNhat
	LNhat =  crossproduct(rvec,rdotvec);
//	std::cout << "LNhat: " << LNhat[0] << "   " << LNhat[1] << "   " << LNhat[2] << std::endl;
	LNhat = LNhat / sqrt(innerproduct(LNhat,LNhat));
	for (i = 0; i < 3; i++)
		if (fabs(LNhat[i]) < 1.0e-6) LNhat[i] = 0.0;
	LNhat = LNhat / sqrt(innerproduct(LNhat,LNhat));
//	std::cout << "LNhat: " << LNhat[0] << "   " << LNhat[1] << "   " << LNhat[2] << std::endl;
	return LNhat;
}





//**********************************************************************//
//
// SphBBHamiltonian member functions:
//
//**********************************************************************//
/**
 * Spherical implementation of the Hamiltonian of Phys.Rev. D81 (2010) 084024
 */
SphBBHamiltonian::SphBBHamiltonian() : BBHamiltonian(), LNhat(3)
{
}

SphBBHamiltonian::SphBBHamiltonian(const BBHamiltonian & targetH)
: BBHamiltonian(), LNhat(3)
{	
	int i,j;
	Array rvec(3), pvec(3), S1(3), S2(3), Xhat(3);
	Array Yhat(3,0.), Zhat(3,0.);
	Array rvecprime(3), pvecprime(3), S1prime(3), S2prime(3);
	Array Xprime(3);
	Array Yprime(3), Zprime(3);
	std::vector<double> variables = targetH.get_variables();
	
// Fixed Cartesian frame
	Xhat[0]=1.;
	Xhat[1]=0.;
	Xhat[2]=0.;
	
	Yhat[0]=0.;
	Yhat[1]=1.;
	Yhat[2]=0.;
	
	Zhat[0]=0.;
	Zhat[1]=0.;
	Zhat[2]=1.;
	
//	std::cout << "Vars wrt fixed Cartesian frame" << std::endl;
//	for(int i=0; i<12;i++)
//		std::cout << variables[i] << std::endl;
//	std::cout << " " << std::endl;
	
//	Organize the vars
	rvec[0]    = variables[0];
	rvec[1]    = variables[1];
	rvec[2]    = variables[2];
	
	pvec[0]    = variables[3];
	pvec[1]    = variables[4];
	pvec[2]    = variables[5];
	
	S1[0]      = variables[6];
	S1[1]      = variables[7];
	S1[2]      = variables[8];
	
	S2[0]      = variables[9];
	S2[1]      = variables[10];
	S2[2]      = variables[11];
	
	LNhat = targetH.evaluate_LNhat(targetH.get_variables());
//	LNhat[0] = 1;
//	LNhat[1] = 1;
//	LNhat[2] = 0;
//	LNhat = LNhat/sqrt(innerproduct(LNhat,LNhat));
//	std::cout << "H before rot " << targetH.evaluate(variables) << std::endl;
//	std::cout << "LNhat " << LNhat << std::endl;
	std::vector< std::vector<double> > Rot1(3), Rot2(3);
	Array LNhatTmp(3);
	if (innerproduct(LNhat,Xhat) < 0.9) {
		Xprime = LNhat;
		Yprime = crossproduct(Xprime,Xhat);
		Yprime = Yprime/sqrt(innerproduct(Yprime,Yprime));
		Zprime = crossproduct(Xprime,Yprime);
		Zprime = Zprime/sqrt(innerproduct(Zprime,Zprime));
		Rot1[0].push_back(1.);
		Rot1[0].push_back(0.);
		Rot1[0].push_back(0.);
		Rot1[1].push_back(0.);
		Rot1[1].push_back(1.);
		Rot1[1].push_back(0.);
		Rot1[2].push_back(0.);
		Rot1[2].push_back(0.);
		Rot1[2].push_back(1.);
	}
	else {
		Rot1[0].push_back(1./sqrt(2.));
		Rot1[0].push_back(-1./sqrt(2.));
		Rot1[0].push_back(0.);
		Rot1[1].push_back(1./sqrt(2.));
		Rot1[1].push_back(1./sqrt(2.));
		Rot1[1].push_back(0.);
		Rot1[2].push_back(0.);
		Rot1[2].push_back(0.);
		Rot1[2].push_back(1.);
		for (i=0; i<3; i++)
			for(j=0; j<3; j++)
				LNhatTmp[i] += Rot1[i][j]*LNhat[j];
		Xprime = LNhatTmp;
		Yprime = crossproduct(Xprime,Xhat);
		Yprime = Yprime/sqrt(innerproduct(Yprime,Yprime));
		Zprime = crossproduct(Xprime,Yprime);
		Zprime = Zprime/sqrt(innerproduct(Zprime,Zprime));		
	}
	Rot2[0].push_back(innerproduct(Xprime,Xhat));
	Rot2[0].push_back(innerproduct(Xprime,Yhat));
	Rot2[0].push_back(innerproduct(Xprime,Zhat));
	Rot2[1].push_back(innerproduct(Yprime,Xhat));
	Rot2[1].push_back(innerproduct(Yprime,Yhat));
	Rot2[1].push_back(innerproduct(Yprime,Zhat));
	Rot2[2].push_back(innerproduct(Zprime,Xhat));
	Rot2[2].push_back(innerproduct(Zprime,Yhat));
	Rot2[2].push_back(innerproduct(Zprime,Zhat));
	
//	std::cout << "Rot1" << std::endl;
//	for (i=0; i<3; i++) {
//		for(j=0; j<3; j++)
//			std::cout << Rot1[i][j] << " ";
//		std::cout << std::endl;
//	}
//	std::cout << "Rot2" << std::endl;
//	for (i=0; i<3; i++) {
//		for(j=0; j<3; j++)
//			std::cout << Rot2[i][j] << " ";
//		std::cout << std::endl;
//	}
//	abort();
		
	Array rvectmp(3), pvectmp(3), S1tmp(3), S2tmp(3);
	Array LNhatprime(3);
	LNhatTmp[0] = 0.;
	LNhatTmp[1] = 0.;
	LNhatTmp[2] = 0.;
	
	for (i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			rvectmp[i] += Rot1[i][j]*rvec[j];
			pvectmp[i] += Rot1[i][j]*pvec[j];
			S1tmp[i] += Rot1[i][j]*S1[j];
			S2tmp[i] += Rot1[i][j]*S2[j];
			LNhatTmp[i] += Rot1[i][j]*LNhat[j];
		}
	for (i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			rvecprime[i] += Rot2[i][j]*rvectmp[j];
			pvecprime[i] += Rot2[i][j]*pvectmp[j];
			S1prime[i] += Rot2[i][j]*S1tmp[j];
			S2prime[i] += Rot2[i][j]*S2tmp[j];
			LNhatprime[i] += Rot2[i][j]*LNhatTmp[j];
		}
//	std::cout << "LNhatprime " << LNhatprime << std::endl;
	variables[0] = rvecprime[0];
	variables[1] = rvecprime[1];
	variables[2] = rvecprime[2];
	
	variables[3] = pvecprime[0];
	variables[4] = pvecprime[1];
	variables[5] = pvecprime[2];
	
	variables[6]  = S1prime[0];
	variables[7]  = S1prime[1];
	variables[8]  = S1prime[2];
	
	variables[9]  = S2prime[0];
	variables[10] = S2prime[1];
	variables[11] = S2prime[2];
	
//	std::cout << "Vars after rotation" << std::endl;
//	for(int i=0; i<12;i++)
//		std::cout << variables[i] << std::endl;
//	std::cout << g" " << std::endl;
//	std::cout << "H after rot " << targetH.evaluate(variables) << std::endl;
	//abort();

//	Move to spherical coordinates wrt to a polar axis along {1,0,0} which is LNhat
	variables[0] = sqrt(innerproduct(rvecprime,rvecprime));
	variables[1] = acos(rvecprime[0] / variables[0]);
	variables[2] = atan2(-rvecprime[1], rvecprime[2]);
	
	variables[3] = innerproduct(rvecprime, pvecprime) / variables[0];
	variables[4] = -innerproduct(crossproduct(crossproduct(rvecprime, Xhat), rvecprime), pvecprime) 
					/ variables[0] / sin(variables[1]);
	variables[5] = -innerproduct(crossproduct(rvecprime, Xhat), pvecprime);
	
//	std::cout << "List of SphH variables:" << std::endl;
//	for (int i = 0; i < 12; i++)
//		std::cout << variables[i] << "   " << std::endl;
//	std::cout << " " << std::endl;
	
	set_flags_ptr(targetH.get_flags_ptr());
	set_wave_params_ptr(targetH.get_wave_params_ptr());
	set_tuning_params_ptr(targetH.get_tuning_params_ptr());
	set_variables(variables);
	set_saves(targetH.get_saves());
}

//SphBBHamiltonian::SphBBHamiltonian(const BBHamiltonian & targetH)
//: BBHamiltonian(), LNhat(3)
//{
//	double r;
//	Array rvec(3), pvec(3), S1(3), S2(3), zvec(3, 0.0), rdotvec(3);
//	Array rvecprime(3), pvecprime(3), S1prime(3), S2prime(3), LNhat(3);
//	Array Xprime(3), Yprime(3), Zprime(3);
//	std::vector<double> variables = targetH.get_variables();
//	
//	
//	rvec[0]    = variables[0];
//	rvec[1]    = variables[1];
//	rvec[2]    = variables[2];
//	
//	pvec[0]    = variables[3];
//	pvec[1]    = variables[4];
//	pvec[2]    = variables[5];
//	
//	S1[0]      = variables[6];
//	S1[1]      = variables[7];
//	S1[2]      = variables[8];
//	
//	S2[0]      = variables[9];
//	S2[1]      = variables[10];
//	S2[2]      = variables[11];
//	//	std::cout << "Passed vars" << std::endl;
//	//	for(int i=0; i<12;i++)
//	//		std::cout << variables[i] << std::endl;
//	//	std::cout << " " << std::endl;
//	//	targetH.set_variables(variables);
//	//	std::cout << "Hpassed = " << temphptr->evaluate(variables) << std::endl;
//	//	std::cout << "rvec = " << rvec << std::endl;
//	//	std::cout << "pvec = " << pvec << std::endl;
//	
//	//Compute rdot vector
//	rdotvec[0] = targetH.evaluate_deriv(3);
//	rdotvec[1] = targetH.evaluate_deriv(4);
//	rdotvec[2] = targetH.evaluate_deriv(5);
//	//	std::cout << "rdotvec = " << rdotvec << std::endl;
//	
//	//Build LNhat
//	LNhat =  crossproduct(rvec,rdotvec);
//	LNhat = LNhat / sqrt(innerproduct(LNhat,LNhat));	
//	//	std::cout << "LNhat = " << LNhat << std::endl;
//	
//	//Define instantaneous orbital plane ref frame
//	r = sqrt(innerproduct(rvec,rvec));
//	Xprime = LNhat;
//	Zprime = rvec/r;
//	Yprime = crossproduct(Zprime,Xprime);
//	
//	//	std::cout << "Instantaneous orbital plane frame" << std::endl;	
////	std::cout << "Xprime = " << Xprime << std::endl;
////	std::cout << "Yprime = " << Yprime << std::endl;
////	std::cout << "Zprime = " << Zprime << std::endl;
//	
//	//Build the variables in the instantaneous orbital plane frame
//	rvecprime[0] = 0.0;
//	rvecprime[1] = 0.0;
//	rvecprime[2] = r;
//	
//	pvecprime[0] = innerproduct(pvec,Xprime);
//	pvecprime[1] = innerproduct(pvec,Yprime);
//	pvecprime[2] = innerproduct(pvec,Zprime);
//	
//	S1prime[0] = innerproduct(S1,Xprime);
//	S1prime[1] = innerproduct(S1,Yprime);
//	S1prime[2] = innerproduct(S1,Zprime);
//	
//	S2prime[0] = innerproduct(S2,Xprime);
//	S2prime[1] = innerproduct(S2,Yprime);
//	S2prime[2] = innerproduct(S2,Zprime);
//	
//	//Arrange everything in the variables array
//	variables[0] = rvecprime[0];
//	variables[1] = rvecprime[1];
//	variables[2] = rvecprime[2];
//	
//	variables[3] = pvecprime[0];
//	variables[4] = pvecprime[1];
//	variables[5] = pvecprime[2];
//	
//	variables[6]  = S1prime[0];
//	variables[7]  = S1prime[1];
//	variables[8]  = S1prime[2];
//	
//	variables[9]  = S2prime[0];
//	variables[10] = S2prime[1];
//	variables[11] = S2prime[2];
//	
//	//	std::cout << "New vars" << std::endl;
//	//	for(int i=0; i<12;i++)
//	//		std::cout << variables[i] << std::endl;
//	
//	
//	rvec[0] =  variables[2];
//	rvec[1] = -variables[1];
//	rvec[2] =  variables[0];
//	pvec[0] =  variables[5];
//	pvec[1] = -variables[4];
//	pvec[2] =  variables[3];
//	zvec[2] = 1.0;
//	
//	unsigned int i;
//	std::cout << "List of CartH variables:" << std::endl;
//	for (i = 0; i < 6; i++)
//	std::cout << variables[i] << "   " << std::endl;
//	std::cout << std::endl;
//	
//	variables[0] = r;
//	variables[1] = acos(rvec[2] / variables[0]);
//	variables[2] = atan2(rvec[1], rvec[0]);
//	variables[3] = innerproduct(rvec, pvec) / variables[0];
//	variables[4] = -innerproduct(crossproduct(crossproduct(rvec, zvec), rvec), pvec) 
//	/ variables[0] / sin(variables[1]);
//	variables[5] = -innerproduct(crossproduct(rvec, zvec), pvec);
//	
//	std::cout << "List of SphH variables:" << std::endl;
//	for (i = 0; i < 6; i++)
//	std::cout << variables[i] << "   " << std::endl;
//	std::cout << std::endl;
//	
//	set_flags_ptr(targetH.get_flags_ptr());
//	set_wave_params_ptr(targetH.get_wave_params_ptr());
//	set_tuning_params_ptr(targetH.get_tuning_params_ptr());
//	set_variables(variables);
//	set_saves(targetH.get_saves());
//}




SphBBHamiltonian::~SphBBHamiltonian()
{
}

/**
 * Compute r for the stored dynamical configuration
 */
double SphBBHamiltonian::evaluate_radius() const
{
	std::vector<double> variables = get_variables();
	return variables[0];
}

double SphBBHamiltonian::evaluate(std::vector<double> variables) const
{
	int i;
	double r, th, ph, pr, pt, pp;
	std::vector<double> cartvar(variables.size());


	r	= variables[0];
	th	= variables[1];
	ph	= variables[2];
	pr	= variables[3];
	pt	= variables[4];
	pp	= variables[5];
	for (i = 6; i < 12; i++)
		cartvar[i] = variables[i];

//	Assume the polar axis along {1,0,0} and move to Cartesian coordinates:
//	these are Cart coords of the instantaneous orb. plane frame
	cartvar[0] = r*cos(th);
	cartvar[1] =-r*sin(th)*sin(ph);
	cartvar[2] = r*sin(th)*cos(ph);
	
	if(th==0. || th==Pi)
	{
		if(th==0.)
		{
			th = Pi/2;
			ph = 0.;
		}
		else
		{
			th = Pi/2;
			ph = Pi;
		}
		cartvar[3] = pr*sin(th)*cos(ph) + pt/r*cos(th)*cos(ph) - pp/r/sin(th)*sin(ph);
		cartvar[4] = pr*sin(th)*sin(ph) + pt/r*cos(th)*sin(ph) + pp/r/sin(th)*cos(ph);
		cartvar[5] = pr*cos(th)         - pt/r*sin(th);		
	}
	else
	{
		cartvar[3] = pr*cos(th) 		-pt/r*sin(th);
		cartvar[4] =-pr*sin(th)*sin(ph) -pt/r*cos(th)*sin(ph) -pp/r/sin(th)*cos(ph);
		cartvar[5] = pr*sin(th)*cos(ph) +pt/r*cos(th)*cos(ph) -pp/r/sin(th)*sin(ph);
	}
	
//	Vars in the Cart inst orb frame
//	std::cout << "Vars in the Cart inst orb frame" << std::endl;
//	for(int i=0; i<12;i++)
//		std::cout << cartvar[i] << std::endl;
//	abort();
	return BBHamiltonian::evaluate(cartvar);
}

//double SphBBHamiltonian::evaluate(std::vector<double> variables) const
//{
//	int i;
//	double r, th, ph, pr, pt, pp;
//	std::vector<double> cartvar(variables.size());
//	
//	r	= variables[0];
//	th	= variables[1];
//	ph	= variables[2];
//	pr	= variables[3];
//	pt	= variables[4];
//	pp	= variables[5];
//	for (i = 6; i < 12; i++)
//		cartvar[i] = variables[i];
//	
//	cartvar[2] = r*sin(th)*cos(ph);
//	cartvar[1] =-r*sin(th)*sin(ph);
//	cartvar[0] = r*cos(th);
//	cartvar[5] = pr*sin(th)*cos(ph) +pt/r*cos(th)*cos(ph) -pp/r/sin(th)*sin(ph);
//	cartvar[4] =-pr*sin(th)*sin(ph) -pt/r*cos(th)*sin(ph) -pp/r/sin(th)*cos(ph);
//	cartvar[3] = pr*cos(th) 		-pt/r*sin(th);
//	
//	//for (i = 0; i < 12; i++)
//	//	std::cout << cartvar[i] << std::endl;
//	
//	return BBHamiltonian::evaluate(cartvar);
//}

double SphBBHamiltonian::nonKeplerianRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	return temphptr -> nonKeplerianRadiusSph();
}

/**
 * Compute dHreal/dr at a given radial separation
 */
double SphBBHamiltonian::dHdr_r(double r /**< Radius */) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	variables[0] = r;	
	variables[3] = 0.0;
	temphptr->set_variables(variables);
	return temphptr->evaluate_deriv(0);
}

/**
 * Compute dHreal/dpphi at a given angular momentum 
 */
double SphBBHamiltonian::dHdr_pp(double pp /**< Angular momentum */) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	variables[3] = 0.0;
	variables[5] = pp;
	temphptr->set_variables(variables);
	return temphptr->evaluate_deriv(0);
}


double SphBBHamiltonian::LightRingRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	return temphptr -> LightRingSph();
}



//**********************************************************************//
//
// SEv2Hamiltonian member functions:
//
//**********************************************************************//


SEv2Hamiltonian::SEv2Hamiltonian() : BBHamiltonian()
{
}

SEv2Hamiltonian::~SEv2Hamiltonian()
{
}

double SEv2Hamiltonian::evaluate_eoba(std::vector<double> variables) const
{
	std::vector<double> saves = get_saves();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	int i;
	double eoba, eta, r, bulk, LogTerms, u;
	double KK, a, deltaU;
//	double KK0, c0, c1, c2;
	double k0, k1, k2, k3, k4, k5, k1p2, k1p3;
	Array rvec(3);
	std::vector<double> SKerr(3);

	eta  = wave_params_ptr->get_eta();

   	for (i = 0; i < 3; i++)
		rvec[i] = variables[i];
	r  = sqrt(innerproduct(rvec, rvec));
	u = 1./r;

	SKerr = evaluate_sigmaKerr(variables);
	a = sqrt(SKerr[0]*SKerr[0]+SKerr[1]*SKerr[1]+SKerr[2]*SKerr[2]);

	//Tuning saves
//	KK0  = saves[0];
//	Ka   = saves[4];
//	wfd1 = saves[5];
//	wfd2 = saves[6];

	//K(eta) polynomial
	KK  = saves[3] + saves[4]*a*a;

	//\Delta_i's
	k0  = KK*(-2. + eta*KK);//k0 or \Delta_1 in BB paper Eq. 5.77
	k1  = -2.*(k0 + KK)*(-1. + eta*KK);//k1 or \Delta_2 in BB paper Eq. 5.78
	k1p2= k1*k1;
	k1p3= k1*k1p2;
    k2  = (k1*(4. + k1 - 4.*eta*KK))/2. - intPower(a,2)*k0*intPower(-1. + eta*KK,2);
    k3  = -intPower(k1,3)/3. + k1*k2 + intPower(k1,2)*(-1. + eta*KK) - 2.*(1. + k2 - eta*KK)*(-1. + eta*KK) - intPower(a,2)*k1*intPower(-1. + eta*KK,2);
    k4  = (24.*intPower(k1,4) - 96.*intPower(k1,2)*k2 + 48.*intPower(k2,2) - 64.*intPower(k1,3)*(-1. + eta*KK)
		   + 48.*intPower(a,2)*(intPower(k1,2) - 2.*k2)*intPower(-1. + eta*KK,2) +
		   96.*k1*(k3 + 2.*k2*(-1. + eta*KK)) - (-1. + eta*KK)*(192.*k3 + (-1. + eta*KK)*(-3008. + 123.*intPower(Pi,2))))/96.;
    k5  = (eta*KK-1)*(eta*KK-1)/eta * (64./5.*eta*log(u) +
    	   eta*(-4237./60.+128./5.*EulerGamma+2275.*Pi*Pi/512.
    			-1./3.*a*a*(k1*k1*k1-3.*k1*k2+3.*k3)
    			-(k1p3*k1p2-5.*k1p3*k2+5.*k1*k2*k2+5.*k1p2*k3-5.*k2*k3-5.*k1*k4)/5./(KK*eta-1.)/(KK*eta-1.)
    			+(k1p2*k1p2-4.*k1p2*k2+2.*k2*k2+4.*k1*k3-4.*k4)/2./(-1.+eta*KK)+256./5.*log(2.)));
	bulk = intPower(-1. + eta*KK,-2) + (2.*u)/(-1. + eta*KK) + intPower(a,2)*intPower(u,2);
    LogTerms = ( 1. + eta*k0 + eta*log(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4) + k5*intPower(u,5)) );

    deltaU = bulk*LogTerms;
	eoba = deltaU;
//	std::cout << k0 << " " << k1 << " " << k2 << " " << k3 << " " << k4 << std::endl;
//	std::cout << "bulk "<< bulk << std::endl;
//	std::cout << "LogTerms " << LogTerms << std::endl;
//	std::cout << "KK " << KK << std::endl;
//	std::cout << "a " << a << std::endl;
//	std::cout << "u " << u << std::endl;
//	std::cout << "eta " << eta << std::endl;
//	std::cout << "EOBA " << eoba << std::endl;
//	abort();
	return eoba;
}

double SEv2Hamiltonian::evaluate_eoba_deriv(std::vector<double> variables) const
{
	std::vector<double> saves = get_saves();
	copy_ptr<WaveformParameters> wave_params_ptr = get_wave_params_ptr();
	int i;
	double eta, r, bulk, LogTerms, u;
	double a;
	//double KK0, c0, c1, c2;
	double KK, k0, k1, k2, k3, k4, k5, k1p2, k1p3;
	Array rvec(3);
	std::vector<double> SKerr(3);

	eta  = wave_params_ptr->get_eta();

   	for (i = 0; i < 3; i++)
		rvec[i] = variables[i];
	r  = sqrt(innerproduct(rvec, rvec));
	u = 1./r;

	SKerr = evaluate_sigmaKerr(variables);
	a = sqrt(SKerr[0]*SKerr[0]+SKerr[1]*SKerr[1]+SKerr[2]*SKerr[2]);

	//Tuning saves
//	KK0  = saves[0];
//	Ka   = saves[4];
//	wfd1 = saves[5];
//	wfd2 = saves[6];

	//K(eta) polynomial
//	c0  = saves[4];
//	c1  = -8.*(c0 - KK0 + Ka*a*a);//c1 or K_1 in BB paper Eq. 6.11
//	c2  = 16.*(c0 - KK0 + Ka*a*a);//c2 or K_2 in BB paper Eq. 6.11
//	KK  = c0 + c1*eta + c2*eta*eta;//KK or K in BB paper Eq. 6.9
	KK = saves[3] + saves[4]*a*a;

	//\Delta_i's
	k0  = KK*(-2. + eta*KK);//k0 or \Delta_1 in BB paper Eq. 5.77
	k1  = -2.*(k0 + KK)*(-1. + eta*KK);//k1 or \Delta_2 in BB paper Eq. 5.78
    k2  = (k1*(4. + k1 - 4.*eta*KK))/2. - intPower(a,2)*k0*intPower(-1. + eta*KK,2);
    k3  = -intPower(k1,3)/3. + k1*k2 + intPower(k1,2)*(-1. + eta*KK) - 2.*(1. + k2 - eta*KK)*(-1. + eta*KK) - intPower(a,2)*k1*intPower(-1. + eta*KK,2);
    k4  = (24.*intPower(k1,4) - 96.*intPower(k1,2)*k2 + 48.*intPower(k2,2) - 64.*intPower(k1,3)*(-1. + eta*KK)
		   + 48.*intPower(a,2)*(intPower(k1,2) - 2.*k2)*intPower(-1. + eta*KK,2) +
		   96.*k1*(k3 + 2.*k2*(-1. + eta*KK)) - (-1. + eta*KK)*(192.*k3 + (-1. + eta*KK)*(-3008. + 123.*intPower(Pi,2))))/96.;
    k1p2= k1*k1;
    k1p3= k1*k1p2;
    k5  = (eta*KK-1)*(eta*KK-1)/eta * (64./5.*eta*log(u) +
    	   eta*(-4237./60.+128./5.*EulerGamma+2275.*Pi*Pi/512.
    			-1./3.*a*a*(k1*k1*k1-3.*k1*k2+3.*k3)
    			-(k1p3*k1p2-5.*k1p3*k2+5.*k1*k2*k2+5.*k1p2*k3-5.*k2*k3-5.*k1*k4)/5./(KK*eta-1.)/(KK*eta-1.)
    			+(k1p2*k1p2-4.*k1p2*k2+2.*k2*k2+4.*k1*k3-4.*k4)/2./(-1.+eta*KK)+256./5.*log(2.)));


	bulk = intPower(-1. + eta*KK,-2) + (2.*u)/(-1. + eta*KK) + intPower(a,2)*intPower(u,2);
    LogTerms = ( 1. + eta*k0 + eta*log(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4) + k5*intPower(u,5)) );

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Implementation where eoba_deriv is deoba/dr
	double dbulkdr, dLogTermsdr, ddeltaUdr, dAdr;
	dbulkdr =  - (2.*u*u)/(-1. + eta*KK) - 2.*intPower(a,2)*intPower(u,3);
	dLogTermsdr = - eta*u*u*(k1 + u*(2.*k2 + u*(3.*k3 + u*(4.*k4+5.*k5*u))))/(1. + k1*u + k2*u*u + k3*u*u*u + k4*u*u*u*u + k5*u*u*u*u*u);

	ddeltaUdr = dbulkdr*LogTerms + bulk*dLogTermsdr;
	dAdr = ddeltaUdr;

	return dAdr;
}

double SEv2Hamiltonian::evaluate_Heff(std::vector<double> variables) const
{

	copy_ptr<Flags>			flags_ptr	  = get_flags_ptr();
	copy_ptr<WaveformParameters> 	wave_params_ptr   = get_wave_params_ptr();
	copy_ptr<TuningParameters> 	tuning_params_ptr = get_tuning_params_ptr();
	std::vector<double> 		saves 	   	  = get_saves();

	Array rvec(3), pvec(3), S1vec(3), S2vec(3), nvec(3);
	std::vector<double> sigmaKerr(3), sigmaStar(3);
	int i;
	double eta = wave_params_ptr->get_eta();
	double a;
	double wfd2, wfd1, d1;

	double r, px, py, pz, nx, ny, nz, Skerr_x, Skerr_y, Skerr_z, costheta, xi2, xi_x, xi_y, xi_z, vx, vy, vz, pxir, pvr, pn, pr, pf, ptheta2,
	w2, rho2, u, deltaT, deltaR, Lambda, D, qq, ww, B, w, MU, nu, BR, wr, nur, mur, wcos, nucos, mucos, ww_r, Lambda_r,
    deltaU, deltaU_u, Q, deltaT_r, pn2, pp, deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z, sx, sy, sz, sxi, sv, sn, s3,
	H, Hns, Hs, Hss, Hwcos, Hwr, HSOL, HSONL, e3_x, e3_y, e3_z;

	double aSK2, r2, Deltat, eobd, eoba, xi, ndotph;

    const double aa=tuning_params_ptr->get_aeta(), bb=tuning_params_ptr->get_beta(), dheff_SS=tuning_params_ptr->get_dheff_SS(); //spin gauge parameters

   	for (i = 0; i < 3; i++)
	{
		rvec[i] = variables[i];
		pvec[i] = variables[i + 3];
		S1vec[i] = variables[i + 6];
		S2vec[i] = variables[i + 9];
	}
//	std::cout << "Heff variables" << std::endl;
//	for(int i=0; i<12; i++)
//		std::cout << variables[i] << std::endl;

	sigmaKerr = evaluate_sigmaKerr(variables);
	sigmaStar = evaluate_sigmaStar(variables);
	a = sqrt(sigmaKerr[0]*sigmaKerr[0]+sigmaKerr[1]*sigmaKerr[1]+sigmaKerr[2]*sigmaKerr[2]);
	aSK2 = a * a;

	r  = sqrt(innerproduct(rvec, rvec));
	r2 = r * r;
	nvec = rvec / r;

//	std::cout << "flags_ptr->get_flag_pr4tortoise() "<< flags_ptr->get_flag_pr4tortoise() << std::endl;
//	std::cout << "flags_ptr->get_flag_tortoiseR() "<< flags_ptr->get_flag_tortoiseR() << std::endl;
	double prold;
	if (flags_ptr->get_flag_pr4tortoise() == 1)
	{
		prold = innerproduct(nvec,pvec);
	}
	else if (flags_ptr->get_flag_pr4tortoise() == 0)
	{
		ndotph = innerproduct(nvec, pvec);
		eoba = evaluate_eoba(variables);
		eobd = evaluate_eobd(variables);
		Deltat	= r2 * eoba;
		//Tortoise coord tranformation
		xi = Deltat / sqrt(eobd) / (r2 + aSK2);
//		std::cout << "xi in Ham= "<< xi <<std::endl;
		//		std::cout << "ndotp Ham " << ndotph <<std::endl;
		if (flags_ptr->get_flag_tortoiseR() == 1)
		{
			prold = innerproduct(nvec,pvec);
			for (i = 0; i < 3; i++)
				pvec[i] = pvec[i] - rvec[i]/r*(xi-1.)/xi*ndotph;
		}
		else if (flags_ptr->get_flag_tortoiseR() == 0)
		{
			prold = xi*innerproduct(nvec,pvec);
		}
	}
	else
	{
		std::cout << "In BBHamiltonian::evaluate_Heff inconsistent flag_tortoiseR" << std::endl;
		abort();
	}
	
//	std::cout << "xi " << xi << std::endl;
//	std::cout << "prold " << prold << std::endl;
//
//	px = pvec[0]/eta;
//	py = pvec[1]/eta;
//	pz = pvec[2]/eta;
	px = pvec[0];
	py = pvec[1];
	pz = pvec[2];
	
//	std::cout << "px = " << px << std::endl;
//	std::cout << "py = " << py << std::endl;
//	std::cout << "pz = " << pz << std::endl;
//	abort();
    nx = nvec[0];
    ny = nvec[1];
    nz = nvec[2];
	
	Array Lhat(3);
	if(innerproduct(pvec,pvec)!=0.){
		Lhat = crossproduct(rvec,pvec);
		Lhat = Lhat/sqrt(innerproduct(Lhat,Lhat));
	}
	else {
		Lhat[0]=0.;
		Lhat[1]=0.;
		Lhat[2]=0.;
	}
    //Tuning saves
	wfd1 = saves[5];
	wfd2 = saves[6];
	d1   = saves[7];

	//BB paper Eq. 5.67 so that Skerr = sigmaKerr
    Skerr_x = sigmaKerr[0];
    Skerr_y = sigmaKerr[1];
    Skerr_z = sigmaKerr[2];

//	std::cout << "nhat = " << nx << " " << ny << " " << nz << std::endl;
	double S1Norm=innerproduct(S1vec,S1vec);
        double S2Norm=innerproduct(S2vec,S2vec);
	double LhatNorm=innerproduct(Lhat,Lhat);

	if(S1Norm==0. && S2Norm==0. && LhatNorm==0.) {
		e3_x = 0.;
		e3_y = 1.;
		e3_z = 0.;
	}
        if(S1Norm==0. && S2Norm==0. && LhatNorm!=0.) {
                e3_x = Lhat[0];
                e3_y = Lhat[1];
                e3_z = Lhat[2];
        }
        if((S1Norm!=0. || S2Norm!=0.) && a!=0.) {
		e3_x = Skerr_x/a;
		e3_y = Skerr_y/a;
		e3_z = Skerr_z/a;
	}
        if((S1Norm!=0. || S2Norm!=0.) && a==0. && LhatNorm!=0) {
                e3_x = Lhat[0];
                e3_y = Lhat[1];
                e3_z = Lhat[2];
        }
        if((S1Norm!=0. || S2Norm!=0.) && a==0. && LhatNorm==0) {
                e3_x = 0.;
                e3_y = 1.;
                e3_z = 0.;
	}

    costheta = e3_x*nx+e3_y*ny+e3_z*nz;
	
    xi2 = 1.-costheta*costheta;
//	std::cout << "xi2 " << xi2 << std::endl;
	
    xi_x =  -e3_z*ny + e3_y*nz;
    xi_y =   e3_z*nx - e3_x*nz;
    xi_z =  -e3_y*nx + e3_x*ny;
	
    vx =  -nz*xi_y + ny*xi_z;
    vy =   nz*xi_x - nx*xi_z;
    vz =  -ny*xi_x + nx*xi_y;
	
    pxir = (px*xi_x+py*xi_y+pz*xi_z)*r;
    pvr  = (px*vx+py*vy+pz*vz)*r;
    pn   = px*nx+py*ny+pz*nz;
	
    pr      = pn;
    pf      = pxir;
    ptheta2 = pvr*pvr/xi2;
	
//	std::cout << "pr "<< pr << std::endl;
//	std::cout << "pf "<< pf << std::endl;
//	std::cout << "ptheta2 "<< ptheta2 << std::endl;
//	std::cout <<  "a " << a << std::endl;
//	abort();
    w2   = r*r + a*a;
    rho2 = r*r + a*a*costheta*costheta;
    u    = 1./r;
	
	deltaU = evaluate_eoba(variables);
    deltaT = r*r*deltaU;
//    std::cout << "bulk = " << bulk << std::endl;
//    std::cout << "LogTerms = " << LogTerms << std::endl;
//	std::cout << "DE = " << deltaT << std::endl;
//	std::cout << "deltaU = " << deltaU << std::endl;
//	std::cout << "eoba = " << evaluate_eoba(variables) << std::endl;

//    deltaU_u = 2.*(1./(-1. + eta*KK) + intPower(a,2)*u)*LogTerms +
//	bulk* (eta*(k1 + u*(2.*k2 + u*(3.*k3 + 4.*k4*u))))/(1. + k1*u + k2*intPower(u,2) + k3*intPower(u,3) + k4*intPower(u,4));
	deltaU_u = -r*r*evaluate_eoba_deriv(variables);
//	deltaU_u = evaluate_eoba_deriv(variables);

    deltaT_r = 2.*r*deltaU-deltaU_u;
	
    Lambda = w2*w2 - a*a*deltaT*xi2;
//	std::cout << "Lambda "<< Lambda << std::endl;
//	abort();
	
///////////////////////////////////////////////////////////////////////////////////////////////
//  Original Enrico's implementdation of 1/eobd
//  D = 1.+log(1. + 6.*eta*u*u + 2.*(26. - 3.*eta)*eta*u*u*u + (36. - Dpot4PN)*eta*eta*u*u*u*u);
///////////////////////////////////////////////////////////////////////////////////////////////
	D = 1./evaluate_eobd(variables);
//	std::cout << "D "<< D << std::endl;
//	abort();
	
    deltaR = deltaT*D;
	
    qq = 2.*eta*(4. - 3.*eta);
	
    ww = 2.*a*r + wfd2*eta*a*a*a*u + wfd1*eta*a*u;
//	std::cout << "ww " << a << " " << eta <<" " <<  wfd1 << " " << wfd2 << " " << ww << std::endl;
//	std::cout << "prold " << prold << std::endl;

    Hns=sqrt(1. - flags_ptr -> get_flag_LR() + intPower(prold,4)*qq/intPower(r,2)
			 + ptheta2/rho2 + intPower(pf,2)*rho2/(Lambda*xi2) + intPower(pr,2)*deltaR/rho2)/
	sqrt(Lambda/(rho2*deltaT)) + pf*ww/Lambda;
//////////////////////////////////////////////////////////////////
//	std::cout << "pr " << pr << std::endl;
//	std::cout << "prold " << prold << std::endl;
//	std::cout << "qq " << qq << std::endl;
//	std::cout << "rho2 " << rho2 << std::endl;
//	std::cout << "pf " << pf << std::endl;
//	std::cout << "Lambda " << Lambda << std::endl;
//	std::cout << "xi2 " << xi2 << std::endl;
//	std::cout << "deltaR " << deltaR << std::endl;
//	std::cout << "intPower(prold,4)*qq/intPower(r,2) " << intPower(prold,4)*qq/intPower(r,2) << std::endl;
//	std::cout << "intPower(pf,2)*rho2/(Lambda*xi2) " << intPower(pf,2)*rho2/(Lambda*xi2) << std::endl;
//	std::cout << "intPower(pr,2)*deltaR/rho2 " << intPower(pr,2)*deltaR/rho2 << std::endl;
//	std::cout << "sqrt(Lambda/(rho2*deltaT)) " << sqrt(Lambda/(rho2*deltaT)) << std::endl;
    B = sqrt(deltaT);
    w = ww/Lambda;
    nu = 0.5*log(deltaT*rho2/Lambda);
    MU = 0.5*log(rho2);

    Lambda_r=4.*r*w2 - a*a*deltaT_r*xi2;

    ww_r=2.*a - (intPower(a,3)*wfd2*eta)/intPower(r,2) - wfd1*eta*a/(r*r);
	
    BR = (-2.*deltaT + sqrt(deltaR)*deltaT_r)/(2.*sqrt(deltaR*deltaT));
    wr = (-Lambda_r*ww + Lambda*ww_r)/(Lambda*Lambda);
    nur = (r/rho2 + (w2 * (-4.*r*deltaT + w2*deltaT_r) ) / (2.*deltaT*Lambda) );
    mur = (r/rho2 - 1./sqrt(deltaR));
	
    wcos = -2.*a*a*costheta*deltaT*ww/(Lambda*Lambda);
    nucos= a*a*costheta*w2*(w2-deltaT)/(rho2*Lambda);
    mucos=a*a*costheta/rho2;
	
    Q=1. + intPower(pvr,2)/(exp(2.*MU)*xi2) + exp(2.*nu)*intPower(pxir,2)/(intPower(B,2)*xi2)
	+ intPower(pn,2)*deltaR/exp(2.*MU);
	
    pn2=pr*pr*deltaR/rho2;
    pp=Q-1.;
	
	deltaSigmaStar_x = 0.;
	deltaSigmaStar_y = 0.;
	deltaSigmaStar_z = 0.;
	
	switch (flags_ptr->get_flag_DeltaSigmaStar_25PN())
	{
		case 1:
		{
			deltaSigmaStar_x +=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sigmaKerr[0] - 8.*bb*(1. + 3.*pn2*r - pp*r)*sigmaStar[0] +
							  eta*(-8.*sigmaKerr[0] - 36.*pn2*r*sigmaKerr[0] + 3.*pp*r*sigmaKerr[0] + 14.*sigmaStar[0]
								   - 30.*pn2*r*sigmaStar[0] + 4.*pp*r*sigmaStar[0]))/(12.*r) + d1*eta*sigmaKerr[0]/r/r/r;

			deltaSigmaStar_y +=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sigmaKerr[1] - 8.*bb*(1. + 3.*pn2*r - pp*r)*sigmaStar[1] +
							  eta*(-8.*sigmaKerr[1] - 36.*pn2*r*sigmaKerr[1] + 3.*pp*r*sigmaKerr[1] + 14.*sigmaStar[1]
								   - 30.*pn2*r*sigmaStar[1] + 4.*pp*r*sigmaStar[1]))/(12.*r) + d1*eta*sigmaKerr[1]/r/r/r;

			deltaSigmaStar_z +=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sigmaKerr[2] - 8.*bb*(1. + 3.*pn2*r - pp*r)*sigmaStar[2] +
							  eta*(-8.*sigmaKerr[2] - 36.*pn2*r*sigmaKerr[2] + 3.*pp*r*sigmaKerr[2] + 14.*sigmaStar[2]
								   - 30.*pn2*r*sigmaStar[2] + 4.*pp*r*sigmaStar[2]))/(12.*r) + d1*eta*sigmaKerr[2]/r/r/r;
			break;
		}
		case 0:
			break;
		default:
		{
			std::cout << "flag_DeltaSigmaStar_15PN must be 0 or 1" << std::endl;
			abort();
			break;
		}
	}

	switch (flags_ptr->get_flag_DeltaSigmaStar_35PN())
	{
		case 1:
		{
			double u1, u2;
			double a1 = 0.;
			double a2 = 0.;
			double a3 = 0.;
			double b1 = 0.;
			double b2 = 0.;
			double b3 = 0.;

//			New mapping	with gauge params
			double aaa = 0.;
			double bbb = 0.;

//			Final mapping with gauge params
			aaa = -3/2*eta;
			bbb = -5/4*eta;
			a1 = eta*eta/2;
			a2 = -(1/8)*eta*(-7 + 8*eta);
			a3 = -((9*eta*eta)/16);
			b1 = 1/16*eta*(9 + 5*eta);
			b2 = -(1/8)*eta*(-17 + 5*eta);
			b3 = -3/8*eta*eta;

                        aaa = 0.;
                        bbb = 0.;
                        a1 = 0.;
                        a2 = 0.;
                        a3 = 0.;
                        b1 = 0.;
                        b2 = 0.;
                        b3 = 0.;
	
			u1 =-(2*(24*b2 + eta*(-353 + 27*eta) + bbb*(56 + 60*eta)) +
				  2*(24*b1 - 24*b2 + bbb*(14 - 66*eta) + 103*eta - 60*intPower(eta,2))*pp*
				  r + 120*(2*b3 - 3*eta*(bbb + eta))*intPower(pn2,2)*intPower(r,2) +
				  (-48*b1 + 4*bbb*(1 + 3*eta) + eta*(23 + 3*eta))*intPower(pp,2)*
				  intPower(r,2) + 6*pn2*r*(16*b1 + 32*b2 + 24*b3 - 47*eta +
										54*intPower(eta,2) + 24*bbb*(1 + eta) +
										(24*b1 - 24*b3 - 16*eta + 21*intPower(eta,2) + bbb*(-2 + 30*eta))*pp*
										r))/(72.*intPower(r,2));

			u2 = (-16*(6*a2 + 7*eta*(8 + 3*eta) + aaa*(14 + 15*eta)) +
				  4*(-24*a1 + 24*a2 - 109*eta + 51*intPower(eta,2) + 2*aaa*(-7 + 33*eta))*
				  pp*r + 30*(-16*a3 + 3*eta*(8*aaa + 9*eta))*intPower(pn2,2)*intPower(r,2) +
				  (96*a1 - 45*eta - 8*aaa*(1 + 3*eta))*intPower(pp,2)*intPower(r,2) -
				  6*pn2*r*(32*a1 + 64*a2 + 48*a3 + 16*eta + 147*intPower(eta,2) +
						   48*aaa*(1 + eta) + (48*a1 - 48*a3 - 6*eta + 39*intPower(eta,2) +
											 aaa*(-4 + 60*eta))*pp*r))/(144.*intPower(r,2));

			deltaSigmaStar_x += u1*sigmaStar[0] + u2*sigmaKerr[0];
			deltaSigmaStar_y += u1*sigmaStar[1] + u2*sigmaKerr[1];
			deltaSigmaStar_z += u1*sigmaStar[2] + u2*sigmaKerr[2];
			break;
		}
		case 0:
			break;
		default:
		{
			std::cout << "flag_DeltaSigmaStar_35PN must be 0 or 1" << std::endl;
			abort();
			break;
		}
	}
	
	//I set Sstar to 0 for the LR
    sx= (1. - flags_ptr->get_flag_LR())*(sigmaStar[0]+deltaSigmaStar_x);
    sy= (1. - flags_ptr->get_flag_LR())*(sigmaStar[1]+deltaSigmaStar_y);
    sz= (1. - flags_ptr->get_flag_LR())*(sigmaStar[2]+deltaSigmaStar_z);

    sxi=sx*xi_x+sy*xi_y+sz*xi_z;
    sv=sx*vx+sy*vy+sz*vz;
    sn=sx*nx+sy*ny+sz*nz;
	
    s3=sx*e3_x+sy*e3_y+sz*e3_z;
	
    Hwr=(exp(-3.*MU - nu)*sqrt(deltaR)*(exp(2.*(MU + nu))*intPower(pxir,2)*sv - B*exp(MU + nu)*pvr*pxir*sxi +
										intPower(B,2)*xi2*(exp(2.*MU)*(sqrt(Q) + Q)*sv + pn*pvr*sn*sqrt(deltaR) - intPower(pn,2)*sv*deltaR)))/(2.*B*(1. + sqrt(Q))*sqrt(Q)*xi2);
	
    Hwcos=(exp(-3.*MU - nu)*(sn*(-(exp(2.*(MU + nu))*intPower(pxir,2)) + intPower(B,2)*(intPower(pvr,2) - exp(2.*MU)*(sqrt(Q) + Q)*xi2)) -
							 B*pn*(B*pvr*sv - exp(MU + nu)*pxir*sxi)*sqrt(deltaR)))/(2.*B*(1. + sqrt(Q))*sqrt(Q));
	
    HSOL=(exp(-MU + 2.*nu)*(-B + exp(MU + nu))*pxir*s3)/(intPower(B,2)*sqrt(Q)*xi2);
	
    HSONL=(exp(-2.*MU + nu)*(-(B*exp(MU + nu)*nucos*pxir*(1. + 2.*sqrt(Q))*sn*xi2) +
							 (-(BR*exp(MU + nu)*pxir*(1. + sqrt(Q))*sv) + B*(exp(MU + nu)*nur*pxir*(1. + 2.*sqrt(Q))*sv + B*mur*pvr*sxi +
																			 B*sxi*(-(mucos*pn*xi2) + sqrt(Q)*(mur*pvr - nur*pvr + (-mucos + nucos)*pn*xi2))))*sqrt(deltaR)))/(intPower(B,2)*(sqrt(Q) + Q)*xi2);
	
    Hs= w*s3 + Hwr*wr + Hwcos*wcos + HSOL + HSONL;
	
    Hss=-0.5/(r*r*r)*(sx*sx+sy*sy+sz*sz-3.*sn*sn);

    H=Hns+Hs+Hss + (1. - flags_ptr->get_flag_LR())*dheff_SS * eta * (innerproduct(S1vec,S1vec)+innerproduct(S2vec,S2vec)) / (r*r*r*r);

//	for(int i=0; i<12; i++)
//	{
//	    std::cout << variables[i] << std::endl;
//	}
//	std::cout <<"SEv2H " <<  H << std::endl;
//	std::cout <<"eta " <<  eta << std::endl;
//	std::cout <<"r " <<  r << std::endl;
//	std::cout <<"Hns " <<  Hns << std::endl;
//	std::cout <<"Hs " <<  Hs << std::endl;
//	std::cout <<"Hss " <<  Hss << std::endl;
//	std::cout <<"dSS term " <<   (1. - flags_ptr->get_flag_LR())*dheff_SS * eta * (sigmaStar[0]*sigmaKerr[0]+sigmaStar[1]*sigmaKerr[1]+sigmaStar[2]*sigmaKerr[2]) / (r*r*r*r) << std::endl;
//	abort();

//	if(r<15.)
//		std::cout << r<< " " << H << " " << deltaU << " " << 1./D << std::endl;

    if(H!=H)
    {
  	    std::cout << "Heff is a nan! H" << std::endl;
	    std::cout << "Heff variables" << std::endl;
	    for(int i=0; i<12; i++)
	    {
	        std::cout << variables[i] << std::endl;
	    }
	    std::cout <<"H " <<  H << std::endl;
	    std::cout <<"eta " <<  eta << std::endl;
	    std::cout <<"r " <<  r << std::endl;
	    std::cout <<"Hns " <<  Hns << std::endl;
	    std::cout <<"Hs " <<  Hs << std::endl;
	    std::cout <<"Hss " <<  Hss << std::endl;
	    std::cout << "S1Norm " << S1Norm << std::endl;
	    std::cout << "S2Norm " << S2Norm << std::endl;
        std::cout << "LhatNorm " << LhatNorm << std::endl;
	    for(int i=0; i<12; i++)
            std::cout << variables[i] << std::endl;
		abort();
    }

	return H;

}

double SEv2Hamiltonian::nonKeplerianRadius() const
{
	SphSEv2Hamiltonian SphspinH(*this);
	return SphspinH.nonKeplerianRadius();
}


//**********************************************************************//
//
// SphSEv2Hamiltonian member functions:
//
//**********************************************************************//
/**
 * Spherical implementation of the Hamiltonian of Phys.Rev. D81 (2010) 084024
 */
SphSEv2Hamiltonian::SphSEv2Hamiltonian() : SEv2Hamiltonian(), LNhat(3)
{
}

SphSEv2Hamiltonian::SphSEv2Hamiltonian(const SEv2Hamiltonian & targetH)
: SEv2Hamiltonian(), LNhat(3)
{
	int i,j;
	Array rvec(3), pvec(3), S1(3), S2(3), Xhat(3);
	Array Yhat(3,0.), Zhat(3,0.);
	Array rvecprime(3), pvecprime(3), S1prime(3), S2prime(3);
	Array Xprime(3);
	Array Yprime(3), Zprime(3);
	std::vector<double> variables = targetH.get_variables();

// Fixed Cartesian frame
	Xhat[0]=1.;
	Xhat[1]=0.;
	Xhat[2]=0.;

	Yhat[0]=0.;
	Yhat[1]=1.;
	Yhat[2]=0.;

	Zhat[0]=0.;
	Zhat[1]=0.;
	Zhat[2]=1.;

//	std::cout << "Vars wrt fixed Cartesian frame" << std::endl;
//	for(int i=0; i<12;i++)
//		std::cout << variables[i] << std::endl;
//	std::cout << " " << std::endl;

//	Organize the vars
	rvec[0]    = variables[0];
	rvec[1]    = variables[1];
	rvec[2]    = variables[2];

	pvec[0]    = variables[3];
	pvec[1]    = variables[4];
	pvec[2]    = variables[5];

	S1[0]      = variables[6];
	S1[1]      = variables[7];
	S1[2]      = variables[8];

	S2[0]      = variables[9];
	S2[1]      = variables[10];
	S2[2]      = variables[11];

	LNhat = targetH.evaluate_LNhat(targetH.get_variables());
//	LNhat[0] = 1;
//	LNhat[1] = 1;
//	LNhat[2] = 0;
//	LNhat = LNhat/sqrt(innerproduct(LNhat,LNhat));
//	std::cout << "H before rot " << targetH.evaluate(variables) << std::endl;
//	std::cout << "LNhat " << LNhat << std::endl;
	std::vector< std::vector<double> > Rot1(3), Rot2(3);
	Array LNhatTmp(3);
	if (innerproduct(LNhat,Xhat) < 0.9) {
		Xprime = LNhat;
		Yprime = crossproduct(Xprime,Xhat);
		Yprime = Yprime/sqrt(innerproduct(Yprime,Yprime));
		Zprime = crossproduct(Xprime,Yprime);
		Zprime = Zprime/sqrt(innerproduct(Zprime,Zprime));
		Rot1[0].push_back(1.);
		Rot1[0].push_back(0.);
		Rot1[0].push_back(0.);
		Rot1[1].push_back(0.);
		Rot1[1].push_back(1.);
		Rot1[1].push_back(0.);
		Rot1[2].push_back(0.);
		Rot1[2].push_back(0.);
		Rot1[2].push_back(1.);
	}
	else {
		Rot1[0].push_back(1./sqrt(2.));
		Rot1[0].push_back(-1./sqrt(2.));
		Rot1[0].push_back(0.);
		Rot1[1].push_back(1./sqrt(2.));
		Rot1[1].push_back(1./sqrt(2.));
		Rot1[1].push_back(0.);
		Rot1[2].push_back(0.);
		Rot1[2].push_back(0.);
		Rot1[2].push_back(1.);
		for (i=0; i<3; i++)
			for(j=0; j<3; j++)
				LNhatTmp[i] += Rot1[i][j]*LNhat[j];
		Xprime = LNhatTmp;
		Yprime = crossproduct(Xprime,Xhat);
		Yprime = Yprime/sqrt(innerproduct(Yprime,Yprime));
		Zprime = crossproduct(Xprime,Yprime);
		Zprime = Zprime/sqrt(innerproduct(Zprime,Zprime));
	}
	Rot2[0].push_back(innerproduct(Xprime,Xhat));
	Rot2[0].push_back(innerproduct(Xprime,Yhat));
	Rot2[0].push_back(innerproduct(Xprime,Zhat));
	Rot2[1].push_back(innerproduct(Yprime,Xhat));
	Rot2[1].push_back(innerproduct(Yprime,Yhat));
	Rot2[1].push_back(innerproduct(Yprime,Zhat));
	Rot2[2].push_back(innerproduct(Zprime,Xhat));
	Rot2[2].push_back(innerproduct(Zprime,Yhat));
	Rot2[2].push_back(innerproduct(Zprime,Zhat));

//	std::cout << "Rot1" << std::endl;
//	for (i=0; i<3; i++) {
//		for(j=0; j<3; j++)
//			std::cout << Rot1[i][j] << " ";
//		std::cout << std::endl;
//	}
//	std::cout << "Rot2" << std::endl;
//	for (i=0; i<3; i++) {
//		for(j=0; j<3; j++)
//			std::cout << Rot2[i][j] << " ";
//		std::cout << std::endl;
//	}
//	abort();

	Array rvectmp(3), pvectmp(3), S1tmp(3), S2tmp(3);
	Array LNhatprime(3);
	LNhatTmp[0] = 0.;
	LNhatTmp[1] = 0.;
	LNhatTmp[2] = 0.;

	for (i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			rvectmp[i] += Rot1[i][j]*rvec[j];
			pvectmp[i] += Rot1[i][j]*pvec[j];
			S1tmp[i] += Rot1[i][j]*S1[j];
			S2tmp[i] += Rot1[i][j]*S2[j];
			LNhatTmp[i] += Rot1[i][j]*LNhat[j];
		}
	for (i=0; i<3; i++)
		for(j=0; j<3; j++)
		{
			rvecprime[i] += Rot2[i][j]*rvectmp[j];
			pvecprime[i] += Rot2[i][j]*pvectmp[j];
			S1prime[i] += Rot2[i][j]*S1tmp[j];
			S2prime[i] += Rot2[i][j]*S2tmp[j];
			LNhatprime[i] += Rot2[i][j]*LNhatTmp[j];
		}
//	std::cout << "LNhatprime " << LNhatprime << std::endl;
	variables[0] = rvecprime[0];
	variables[1] = rvecprime[1];
	variables[2] = rvecprime[2];

	variables[3] = pvecprime[0];
	variables[4] = pvecprime[1];
	variables[5] = pvecprime[2];

	variables[6]  = S1prime[0];
	variables[7]  = S1prime[1];
	variables[8]  = S1prime[2];

	variables[9]  = S2prime[0];
	variables[10] = S2prime[1];
	variables[11] = S2prime[2];

//	std::cout << "Vars after rotation" << std::endl;
//	for(int i=0; i<12;i++)
//		std::cout << variables[i] << std::endl;
//	std::cout << " " << std::endl;
//	std::cout << "H after rot " << targetH.evaluate(variables) << std::endl;
	//abort();

//	Move to spherical coordinates wrt to a polar axis along {1,0,0} which is LNhat
	variables[0] = sqrt(innerproduct(rvecprime,rvecprime));
	variables[1] = acos(rvecprime[0] / variables[0]);
	variables[2] = atan2(-rvecprime[1], rvecprime[2]);

	variables[3] = innerproduct(rvecprime, pvecprime) / variables[0];
	variables[4] = -innerproduct(crossproduct(crossproduct(rvecprime, Xhat), rvecprime), pvecprime)
					/ variables[0] / sin(variables[1]);
	variables[5] = -innerproduct(crossproduct(rvecprime, Xhat), pvecprime);

//	std::cout << "List of SphH variables:" << std::endl;
//	for (int i = 0; i < 12; i++)
//		std::cout << variables[i] << "   " << std::endl;
//	std::cout << " " << std::endl;

	set_flags_ptr(targetH.get_flags_ptr());
	set_wave_params_ptr(targetH.get_wave_params_ptr());
	set_tuning_params_ptr(targetH.get_tuning_params_ptr());
	set_variables(variables);
	set_saves(targetH.get_saves());
}

SphSEv2Hamiltonian::~SphSEv2Hamiltonian()
{
}

/**
 * Compute r for the stored dynamical configuration
 */
double SphSEv2Hamiltonian::evaluate_radius() const
{
	std::vector<double> variables = get_variables();
	return variables[0];
}

double SphSEv2Hamiltonian::evaluate(std::vector<double> variables) const
{
	int i;
	double r, th, ph, pr, pt, pp;
	std::vector<double> cartvar(variables.size());


	r	= variables[0];
	th	= variables[1];
	ph	= variables[2];
	pr	= variables[3];
	pt	= variables[4];
	pp	= variables[5];
	for (i = 6; i < 12; i++)
		cartvar[i] = variables[i];

//	Assume the polar axis along {1,0,0} and move to Cartesian coordinates:
//	these are Cart coords of the instantaneous orb. plane frame
	cartvar[0] = r*cos(th);
	cartvar[1] =-r*sin(th)*sin(ph);
	cartvar[2] = r*sin(th)*cos(ph);

	if(th==0. || th==Pi)
	{
		if(th==0.)
		{
			th = Pi/2;
			ph = 0.;
		}
		else
		{
			th = Pi/2;
			ph = Pi;
		}
		cartvar[3] = pr*sin(th)*cos(ph) + pt/r*cos(th)*cos(ph) - pp/r/sin(th)*sin(ph);
		cartvar[4] = pr*sin(th)*sin(ph) + pt/r*cos(th)*sin(ph) + pp/r/sin(th)*cos(ph);
		cartvar[5] = pr*cos(th)         - pt/r*sin(th);
	}
	else
	{
		cartvar[3] = pr*cos(th) 		-pt/r*sin(th);
		cartvar[4] =-pr*sin(th)*sin(ph) -pt/r*cos(th)*sin(ph) -pp/r/sin(th)*cos(ph);
		cartvar[5] = pr*sin(th)*cos(ph) +pt/r*cos(th)*cos(ph) -pp/r/sin(th)*sin(ph);
	}

//	Vars in the Cart inst orb frame
//	std::cout << "Vars in the Cart inst orb frame" << std::endl;
//	for(int i=0; i<12;i++)
//		std::cout << cartvar[i] << std::endl;
//	abort();
	return SEv2Hamiltonian::evaluate(cartvar);
}

double SphSEv2Hamiltonian::nonKeplerianRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	return temphptr -> nonKeplerianRadiusSph();
}

/**
 * Compute dHreal/dr at a given radial separation
 */
double SphSEv2Hamiltonian::dHdr_r(double r /**< Radius */) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	variables[0] = r;
	variables[3] = 0.0;
	temphptr->set_variables(variables);
	return temphptr->evaluate_deriv(0);
}

/**
 * Compute dHreal/dpphi at a given angular momentum
 */
double SphSEv2Hamiltonian::dHdr_pp(double pp /**< Angular momentum */) const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	std::vector<double> variables = get_variables();
	variables[3] = 0.0;
	variables[5] = pp;
	temphptr->set_variables(variables);
	return temphptr->evaluate_deriv(0);
}


double SphSEv2Hamiltonian::LightRingRadius() const
{
	copy_ptr<Hamiltonian> temphptr(this->clone());
	return temphptr -> LightRingSph();
}



/* ==============================
 *
 *  Non-member functions
 *
 * ==============================
 */

/**
 * Non-member function that returns the derivative of a Hamiltonian object w.r.t. one variable by Ridders' method of polynomial extrapolation
 */
double dfridr_hamiltonian(const Hamiltonian & targetH,	/**< Hamiltonian object */ 
						  int n_x,						/**< Choose the variable for differentiation */
						  std::vector<double> xvec,		/**< Array of the dynamical configuration */
						  double h,						/**< Step size for the numerical derivative */
						  double *err					/**< OUTPUT: Error estimate */
						  )
{
	const int NTAB = 10;
	const double CON = 1.4, CON2 = CON * CON, BIG = 1.0e30;
	int i, j;
	double errt, fac, hh, ans;
	std::vector<double> pxvec = xvec, mxvec = xvec;
	Matrix a(NTAB, NTAB);

	if (h == 0.0)
	{
		std::cout << "Error! Hamiltonian::dfridr_hamiltonian: ";
		std::cout << "initial step size is zero" << std::endl;
		abort();
	}
	hh = h;
	pxvec[n_x] = xvec[n_x] + hh;
	mxvec[n_x] = xvec[n_x] - hh;
	a[0][0] = targetH.evaluate(pxvec) - targetH.evaluate(mxvec) / (2.0 * hh);

	*err = BIG;
	for (i = 1; i < NTAB; i++)
	{
		hh /= CON;
		pxvec[n_x] = xvec[n_x] + hh;
		mxvec[n_x] = xvec[n_x] - hh;
		a[0][i] = (targetH.evaluate(pxvec)
				-  targetH.evaluate(mxvec)) / (2.0 * hh);
		fac = CON2;
		for (j = 1; j < i; j++)
		{
			a[j][i] = (a[j-1][i]*fac - a[j-1][i-1]) / (fac-1.0);
			fac *= CON2;
			errt = std::max(fabs(a[j][i]-a[j-1][i]),
							fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err)
			{
				*err = errt;
				ans = a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= std::max((*err), 1.0e-16))
		{
			std::cout << "Warning! Hamiltonian::dfridr_hamiltonian: ";
			std::cout << "numerical error is ";
			std::cout << fabs(a[i][i]-a[i-1][i-1]) << std::endl;
		}
	}
	return ans;
}

/**
* Non-member function that takes Hamiltonian objects and an initial guess, and returns a range in pphi (looking outward) that brackets a root of dHreal/dr = 0
**/
int zbrac_hamiltonian(const Hamiltonian & targetH,	/**< Hamiltonian object */
					  double *x1,					/**< INPUT: initial lower bound for the range, OUTPUT: lower bound on the root */
					  double *x2					/**< INPUT: initial upper bound for the range, OUTPUT: upper bound on the root */
					  )
{
	const double	FACTOR = 1.6;
	const int 	NTRY = 50;
	int j;
	double f1, f2;
	if (*x1 == *x2) 
	{
		std::cout << "Warning! Hamiltonian::zbrac_hamiltonian: ";
		std::cout << "Bad initial range in zbrac" << std::endl;
	}
	f1 = targetH.dHdr_pp(*x1);
	f2 = targetH.dHdr_pp(*x2);
	for (j = 1; j <= NTRY; j++)
	{
		if (f1*f2 < 0.0) return 1;
		if (fabs(f1) < fabs(f2))
		{
			if (*x1 + FACTOR * (*x1 - *x2) > 0)
				f1 = targetH.dHdr_pp(*x1 += FACTOR * (*x1 - *x2));
			else
				f1 = targetH.dHdr_pp(*x1 += -*x1);			
		}
		else
			f2 = targetH.dHdr_pp(*x2 += FACTOR * (*x2 - *x1));
	}
	return 0;
}

/**
 * Non-member function that takes Hamiltonian objects and an initial guess, and returns ranges in pphi (looking inward) that bracket multiple roots of dH/dr = 0
 */
void zbrak_hamiltonian(const Hamiltonian & targetH,		/**< Hamiltonian object */
					   double x1,						/**< Lower bound for the search interval */
					   double x2,						/**< Upper bound for the search interval */
					   int n,							/**< This decides the step dx=(x2-x1)/n */
					   std::vector<double> & xb1,		/**< OUTPUT: Array of lower bounds on the found roots */
					   std::vector<double> & xb2,		/**< OUTPUT: Array of upper bounds on the found roots*/
					   int *nb							/**< OUTPUT: Number of roots found */
					   )
{
	int nbb, i;
	double x, fp, fc, dx;
	nbb = 0;
	dx = (x2 - x1)/n;
	fp = targetH.dHdr_pp(x = x1);
	for (i = 1; i <= n; i++)
	{
		fc = targetH.dHdr_pp(x += dx);
		if (fc * fp <= 0.0)
		{
			xb1.push_back(x-dx);
			xb2.push_back(x);
			nbb++;
			if (*nb == nbb) return;
		}
		fp = fc;
	}
	*nb = nbb;
	return;
}

/**
 * Non-member function that takes Hamiltonian objects and returns pphi that satisfies dH/dr = 0
 */
double zbrent_hamiltonian(const Hamiltonian & targetH,	/**< Hamiltonian object */
						  double x1,					/**< Lower bound on the root */
						  double x2,					/**< Upper bound on the root */
						  double tol					/**< Tolerance */
						  )
{
	const int ITMAX = 1000;
	const double EPS = 3.0e-15;
	int iter;
	double a = x1, b = x2, c = x2, d, e, min1, min2;
	double fa = targetH.dHdr_pp(a), fb = targetH.dHdr_pp(b), fc, p, q, r, s, tol1, xm;
	
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
	{
		std::cout << "Error! Hamiltonian::zbrent_hamiltonian: ";
		std::cout << "root must be bracketed. (numerical information: ";
		std::cout << a << " --> " << fa << "   " << b << " --> " << fb;
		std::cout << ")" << std::endl;
		std::cout << "Evolution stopped at r = ";
		std::cout << targetH.evaluate_radius() << std::endl;
		for(int k=0;k<12; k++)
			std::cout << targetH.get_variables()[k] << std::endl;
		abort();
	}
	
	fc = fb;
	for (iter = 0; iter < ITMAX; iter++)
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
		{
			c = a;
			fc = fa;
			e = d = b - a;
		}
		if (fabs(fc) < fabs(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
		xm = 0.5 * (c - b);
		if (fabs(xm) <= tol1 || fb == 0.0)
			return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
		{
			s = fb / fa;
			if (a == c)
			{
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}
			else
			{
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (p > 0.0)
				q = -q;
			p = fabs(p);
			min1 = 3.0 * xm * q - fabs(tol1 * q);
			min2 = fabs(e * q);
			if (2.0 * p < (min1 < min2 ? min1 : min2))
			{
				e = d;
				d = p / q;
			}
			else
			{
				d = xm;
				e = d;
			}
		}
		else
		{
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += (xm >= 0.0 ? fabs(tol1) : -fabs(tol1));
		fb = targetH.dHdr_pp(b);
	}
	std::cout << "Error! Hamiltonian::zbrent_hamiltonian: ";
	std::cout << "maximum number of iterations exceeded." << std::endl;
	return 0.0;
}

/**
 * Non-member function that takes Hamiltonian objects and returns r1 and r2 (look outward) that bracket a root of dH/dr = 0
 */
int zbrac_LR(const Hamiltonian & targetH,	/**< Hamiltonian object */
			 double *x1,					/**< INPUT: initial lower bound for the range, OUTPUT: lower bound on the root */
			 double *x2						/**< INPUT: initial upper bound for the range, OUTPUT: upper bound on the root */
			 )
{
	const double	FACTOR = 1.6;
	const int 	NTRY = 50;
	int j;
	double f1, f2;
	if (*x1 == *x2) 
	{
		std::cout << "Warning! Hamiltonian::zbrac_hamiltonian: ";
		std::cout << "Bad initial range in zbrac" << std::endl;
	}
	f1 = targetH.dHdr_r(*x1);
	f2 = targetH.dHdr_r(*x2);
	for (j = 1; j <= NTRY; j++)
	{
		if (f1*f2 < 0.0) return 1;
		if (fabs(f1) < fabs(f2))
		{
			if (*x1 + FACTOR * (*x1 - *x2) > 0)
				f1 = targetH.dHdr_r(*x1 += FACTOR * (*x1 - *x2));
			else
				f1 = targetH.dHdr_r(*x1 += -*x1);			
		}
		else
			f2 = targetH.dHdr_r(*x2 += FACTOR * (*x2 - *x1));
	}
	return 0;
}

/**
 * Non-member function that takes Hamiltonian objectsand returns ranges in r (looking inward) that bracket multiple roots of dH/dr = 0
 */
void zbrak_LR(const Hamiltonian & targetH,		/**< Hamiltonian object */
			  double x1,						/**< Lower bound for the search interval */
			  double x2,						/**< Upper bound for the search interval */
			  int n,							/**< This decides the step dx=(x2-x1)/n */
			  std::vector<double> & xb1,		/**< OUTPUT: Array of lower bounds on the found roots */
			  std::vector<double> & xb2,		/**< OUTPUT: Array of upper bounds on the found roots*/
			  int *nb							/**< OUTPUT: Number of roots found */
			  )
{
	int nbb, i;
	double x, fp, fc, dx;
	nbb = 0;
	dx = (x2 - x1)/n;
	fp = targetH.dHdr_r(x = x2);
	for (i = n; i >= 1; i--)
	{
		fc = targetH.dHdr_r(x -= dx);
//		std::cout << "x = " << x << std::endl;
		if (fc * fp <= 0.0)
		{
			xb1.push_back(x);
			xb2.push_back(x+dx);
			nbb++;
			return;
//			if (*nb == nbb) return;
		}
		fp = fc;
	}
	*nb = nbb;
	return;
}

/**
 * Non-member function that takes Hamiltonian objects and returns r that satisfies dH/dr = 0
 */
double zbrent_LR(const Hamiltonian & targetH,	/**< Hamiltonian object */
				 double x1,					/**< Lower bound on the root */
				 double x2,					/**< Upper bound on the root */
				 double tol					/**< Tolerance */
				 )
{
	const int ITMAX = 1000;
	const double EPS = 3.0e-15;
	int iter;
	double a = x1, b = x2, c = x2, d, e, min1, min2;
	double fa = targetH.dHdr_r(a), fb = targetH.dHdr_r(b), fc, p, q, r, s, tol1, xm;
	
	
//	for (int k=0; k<12; k++)
//	{
//		std::cout << "In" << std::endl;
//		std::cout << targetH.get_variables()[k] << std::endl;
//	}	
//	std::cout << a << " --> " << fa << "   " << b << " --> " << fb;
//	
	
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
	{
		std::cout << "Error! Hamiltonian::zbrent_LR: ";
		std::cout << "root must be bracketed. (numerical information: ";
		std::cout << a << " --> " << fa << "   " << b << " --> " << fb;
		std::cout << ")" << std::endl;
		std::cout << "Evolution stopped at r = ";
		std::cout << targetH.evaluate_radius() << std::endl;
		abort();
	}
	
	fc = fb;
	for (iter = 0; iter < ITMAX; iter++)
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
		{
			c = a;
			fc = fa;
			e = d = b - a;
		}
		if (fabs(fc) < fabs(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
		xm = 0.5 * (c - b);
		if (fabs(xm) <= tol1 || fb == 0.0)
			return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
		{
			s = fb / fa;
			if (a == c)
			{
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}
			else
			{
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (p > 0.0)
				q = -q;
			p = fabs(p);
			min1 = 3.0 * xm * q - fabs(tol1 * q);
			min2 = fabs(e * q);
			if (2.0 * p < (min1 < min2 ? min1 : min2))
			{
				e = d;
				d = p / q;
			}
			else
			{
				d = xm;
				e = d;
			}
		}
		else
		{
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += (xm >= 0.0 ? fabs(tol1) : -fabs(tol1));
		fb = targetH.dHdr_r(b);
	}
	std::cout << "Error! Hamiltonian::zbrent_LR: ";
	std::cout << "maximum number of iterations exceeded." << std::endl;
	return 0.0;
}
