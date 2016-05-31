#ifndef laminarSMOKEthermoClass_H
#define laminarSMOKEthermoClass_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "fvMatricesFwd.H"
#include "Switch.H"
#include <vector>
#include <string>

namespace Foam
{

	class laminarSMOKEthermoClass : public IOdictionary
	{
	public:

		//- Create IO object if dictionary is present
		IOobject createIOobject(const fvMesh& mesh) const;
		
		//- Return a reference to the selected radiation model
		static autoPtr<laminarSMOKEthermoClass> New(const volScalarField& T, const volScalarField& p, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names);

		// Constructors
		laminarSMOKEthermoClass(const volScalarField& T, const volScalarField& p, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names);

	public:

		const volScalarField& T() const { return T_; }
		const volScalarField& p() const { return p_; }
		const PtrList<volScalarField>& Y() const { return Y_; }
		const volScalarField& Y(const label i) const { return Y_[i]; }
		double W(const label i) const { return W_[i]; }
		unsigned int species_index(const std::string name) const;
	
	private:

		const volScalarField& T_;
		const volScalarField& p_;
		const PtrList<volScalarField>& Y_;
		std::vector<double> W_;
		std::vector<std::string> species_names_;
	};
}

#endif
