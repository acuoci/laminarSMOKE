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

	void ReadFromFile(const std::string file_name) 
	{
		std::ifstream fInput(file_name.c_str(), std::ios::in);

		std::string line;
		while (std::getline(fInput, line))
		{
			typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
			boost::char_separator<char> sep(" ");
			tokenizer tok(line,sep);
			
			std::vector<int> indices;

			int count = 1;
			for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
			{
				if (count == 1)
				{
					class_name_.push_back(*it);
				}
				else
				{
					int index = boost::lexical_cast<int>(*it);
					indices.push_back(index);
				}

				count++;
			}

			reaction_indices_.push_back(indices);
		}

		number_of_classes_ = reaction_indices_.size();
		number_reactions_per_class_.resize(number_of_classes_);

		for (int i = 0; i < number_of_classes_; i++)
			number_reactions_per_class_[i] = reaction_indices_[i].size();

		for (int i = 0; i < number_of_classes_; i++)
			std::cout << class_name_[i] << " " << number_reactions_per_class_[i] << std::endl;
	}

	int number_of_classes() const { return number_of_classes_; }

	std::string	class_name(const int i) const { return class_name_[i]; }

	int	number_reactions_per_class(const int i) const { return number_reactions_per_class_[i]; }
	
	const std::vector<int>&	reaction_indices(const int i) const { return reaction_indices_[i]; }

private:

	std::vector<std::string>		class_name_;
	std::vector< std::vector<int> >	reaction_indices_;
	int								number_of_classes_;
	std::vector<int>				number_reactions_per_class_;

};
