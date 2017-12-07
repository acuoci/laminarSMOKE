#include "laminarSMOKEthermoClass.H"

namespace Foam
{

	Foam::autoPtr<laminarSMOKEthermoClass> laminarSMOKEthermoClass::New ( const volScalarField& T, const volScalarField& p, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names )
	{
	    IOobject radIO
	    (
		"laminarSMOKEthermoClass",
		T.time().constant(),
		T.mesh(),
		IOobject::NO_READ,
		IOobject::NO_WRITE,
		true
	    );

	    return autoPtr<laminarSMOKEthermoClass>(new laminarSMOKEthermoClass(T,p,Y,W,species_names));
	}

	Foam::IOobject laminarSMOKEthermoClass::createIOobject ( const fvMesh& mesh ) const
	{
	    IOobject io
	    (
		"laminarSMOKEthermoClass",
		mesh.time().constant(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    );

	    // check if field exists and can be read
            #if DEVVERSION == 1
	    if (io.typeHeaderOk<volScalarField>(true))
	    #else
	    if (io.headerOk())
	    #endif
	    {
		io.readOpt() = IOobject::NO_READ;
		return io;
	    }
	    else
	    {
		io.readOpt() = IOobject::NO_READ;
		return io;
	    }
	}


	laminarSMOKEthermoClass::laminarSMOKEthermoClass(const volScalarField& T, const volScalarField& p, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names) :
		IOdictionary
		    (
			IOobject
			(
			    "laminarSMOKEthermoClass",
			    T.time().constant(),
			    T.mesh(),
			    IOobject::NO_READ,
			    IOobject::NO_WRITE
			)
		    ),
		T_(T),
		p_(p),
		Y_(Y),
		W_(W),
		species_names_(species_names)
	{
		Info << "Created object laminarSMOKEthermoClass" << endl;
	}

	unsigned int laminarSMOKEthermoClass::species_index(const std::string name) const
	{
		for (unsigned int i=0;i<species_names_.size();i++)
		{
			if (name == species_names_[i])
				return i;
		}
	
		Info << "Species " << name << " not available!" << endl;
		abort(FatalError);

		return 0; 
	}
}
