/*-----------------------------------------------------------------------*\
|                                                                         |
|                    ╔═══╦═╗╔═╦═══╦╗╔═╦═══╗                               |
|                    ║╔═╗║║╚╝║║╔═╗║║║╔╣╔══╝                               | 
|   ╔╗╔══╦╗╔╦╦═╗╔══╦═╣╚══╣╔╗╔╗║║ ║║╚╝╝║╚══╗                               |
|   ║║║╔╗║╚╝╠╣╔╗╣╔╗║╔╩══╗║║║║║║║ ║║╔╗║║╔══╝                               |
|   ║╚╣╔╗║║║║║║║║╔╗║║║╚═╝║║║║║║╚═╝║║║╚╣╚══╗                               |
|   ╚═╩╝╚╩╩╩╩╩╝╚╩╝╚╩╝╚═══╩╝╚╝╚╩═══╩╝╚═╩═══╝                               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci                                                     |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of laminarSMOKE solver.                             |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2016, 2015, 2014 A. Cuoci                                |
|   laminarSMOKE is free software: you can redistribute it and/or modify  |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   laminarSMOKE is distributed in the hope that it will be useful,       |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with laminarSMOKE. If not, see <http://www.gnu.org/licenses/>.  |
|                                                                         |
\*-----------------------------------------------------------------------*/

class sparkModel
{
public: 

	sparkModel()
	{
		spark_ = false;
		ignitionTime_ = 0.;
		duration_ = 0.;
		ignitionTemperature_ = 0.;
		ignitionDiameter_ = 0.;
		nIgnCells_ = 0;
	}
	
	bool spark() const { return spark_; }
	Foam::vector position() const { return position_; }
	scalar ignitionTime() const { return ignitionTime_; }
	scalar duration() const { return duration_; }
	scalar ignitionTemperature() const { return ignitionTemperature_; }
	scalar ignitionDiameter() const { return ignitionDiameter_; }
	
	void Read(const dictionary& sparkDictionary)
	{
		spark_ = Switch(sparkDictionary.lookup(word("spark")));
		
		if(spark_ == true)
		{
			position_ = sparkDictionary.lookup("position");
			ignitionTime_ = readScalar(sparkDictionary.lookup("time"));
			duration_ = readScalar(sparkDictionary.lookup("duration"));
			ignitionTemperature_ = readScalar(sparkDictionary.lookup("temperature"));
			ignitionDiameter_ = readScalar(sparkDictionary.lookup("diameter"));   
		}
	}
	
	void Setup(const fvMesh& mesh)
	{
		if (spark_ == true)
		{
			const volVectorField& centres = mesh.C();

			label ignitionCell;
			ignitionCell = mesh.findCell(position_);
			
			if (ignitionCell == -1)
			{
				Info << "Warning: theorical ignition position does not correspond to any cell." << endl;
			}
			else 
			{
				Info << "Ignition cell center: " << centres[ignitionCell] << ". Original position: " << position_ <<endl;

				scalar radius = ignitionDiameter_/2.0;
				ignCellsList_.setSize(1);
				ignCellsList_[0] = ignitionCell;

				forAll(centres, celli)
				{
					scalar dist = mag(centres[celli] - position_);

					if (dist < radius && celli != ignitionCell)
					{
						ignCellsList_.setSize(nIgnCells_+1);
						ignCellsList_[nIgnCells_] = celli;
						nIgnCells_++;
					}
				}

				if (ignCellsList_.size() > 0) 
					Info << "Added " << nIgnCells_ << " cells to ignition point" << endl;
				else
					Info << "The solver wasn't able to find more than one cell in the ignition point" <<endl;
			}
		}
	}
	
	void SetSpark(const double t, volScalarField& T)
	{
		if( (spark_ == true) && (t >= ignitionTime_) && (t <= (ignitionTime_+duration_)) ) 
			for(int i = 0; i<nIgnCells_; i++)
					T[ignCellsList_[i]] = ignitionTemperature_;
	}
	
private:
	Switch spark_;
	Foam::vector position_;
	scalar ignitionTime_;
	scalar duration_;
	scalar ignitionTemperature_;
	scalar ignitionDiameter_;
	labelList ignCellsList_;
    	label nIgnCells_;
};
