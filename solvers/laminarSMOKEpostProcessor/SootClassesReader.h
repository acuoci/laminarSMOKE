#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>

class SootClassesReader
{

public:

	void ReadFromFile(rapidxml::xml_document<>& doc) 
	{

		rapidxml::xml_node<>* main_node = doc.first_node("ReacOfProc");
		if (main_node == 0)
		{
		//	ErrorMessage("Error");
			std::cout << "Error1" << std::endl;
			exit(-1);
		}

		// Loop over processes
		for (rapidxml::xml_node<> *proc_node = main_node->first_node("proc"); proc_node; proc_node = proc_node->next_sibling())
		{
			// std::stringstream fInput;
			class_name_.push_back(proc_node->first_attribute("name")->value());

			std::vector<int> reaction_indices;
			std::vector<double> mass_gains;
			std::vector<double> particle_number_gains;
			for (rapidxml::xml_node<> *reaction_node = proc_node->first_node("reac"); reaction_node; reaction_node = reaction_node->next_sibling())
			{
				reaction_indices.push_back(std::atoi(reaction_node->first_attribute("reacnum")->value()));
				mass_gains.push_back(std::atof(reaction_node->first_attribute("massgain")->value()));
				particle_number_gains.push_back(std::atof(reaction_node->first_attribute("partnumgain")->value()));
			}

			reaction_indices_.push_back(reaction_indices);
			mass_gains_.push_back(mass_gains);
			particle_number_gains_.push_back( particle_number_gains);

			number_reactions_per_class_.push_back(reaction_indices.size());
		}

		number_of_classes_ = class_name_.size();

		for (int i = 0; i < number_of_classes_; i++)
			std::cout << class_name_[i] << " " << number_reactions_per_class_[i] << std::endl;
	}


	int number_of_classes() const { return number_of_classes_; }

	std::string	class_name(const int i) const { return class_name_[i]; }

	int	number_reactions_per_class(const int i) const { return number_reactions_per_class_[i]; }
	
	const std::vector<int>&		reaction_indices(const int i) const { return reaction_indices_[i]; }
	const std::vector<double>&	mass_gains(const int i) const { return mass_gains_[i]; }
	const std::vector<double>&	particle_number_gains(const int i) const { return particle_number_gains_[i]; }

private:

	std::vector<std::string>		class_name_;
	std::vector< std::vector<int> >		reaction_indices_;
	std::vector< std::vector<double> >	mass_gains_;
	std::vector< std::vector<double> >	particle_number_gains_;
	int					number_of_classes_;
	std::vector<int>			number_reactions_per_class_;

};
