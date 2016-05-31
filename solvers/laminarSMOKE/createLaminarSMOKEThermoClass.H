// Create laminarSMOKEthermoClass object
std::vector<double> W(thermodynamicsMapXML->NumberOfSpecies());
std::vector<std::string> species_names(thermodynamicsMapXML->NumberOfSpecies());

for (unsigned int i=0;i<species_names.size();i++)
{
	species_names[i] = thermodynamicsMapXML->NamesOfSpecies()[i];
	W[i] = thermodynamicsMapXML->MW()[i+1];
}

autoPtr<Foam::laminarSMOKEthermoClass> myLaminarSMOKEthermo(Foam::laminarSMOKEthermoClass::New(T,p,Y, W, species_names));
