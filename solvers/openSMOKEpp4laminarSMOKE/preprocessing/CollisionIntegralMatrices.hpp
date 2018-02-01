/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_CollisionIntegralMatrices_H
#define OpenSMOKE_CollisionIntegralMatrices_H

namespace OpenSMOKE
{
	class CollisionIntegralMatrices
	{
	public:

		static const Eigen::VectorXd deltaStar;
		static const Eigen::VectorXd TStar;
		static const Eigen::MatrixXd O37_8;
		static const Eigen::MatrixXd P37_8;

		static const double FITASTAR_1;
		static const double FITASTAR_2;
		static const double FITASTAR_3;
		static const double FITASTAR_4;
		static const double FITASTAR_5;
		static const double FITASTAR_6;
		static const double FITASTAR_7;

		static const double FITBSTAR_1;
		static const double FITBSTAR_2;
		static const double FITBSTAR_3;
		static const double FITBSTAR_4;
		static const double FITBSTAR_5;
		static const double FITBSTAR_6;
		static const double FITBSTAR_7;

		static const double FITCSTAR_1;
		static const double FITCSTAR_2;
		static const double FITCSTAR_3;
		static const double FITCSTAR_4;
		static const double FITCSTAR_5;
		static const double FITCSTAR_6;
		static const double FITCSTAR_7;
	};

	const Eigen::VectorXd CollisionIntegralMatrices::deltaStar = (Eigen::VectorXd(8) <<
		0., .25, .5, .75, 1., 1.5, 2., 2.5).finished();

	const Eigen::VectorXd CollisionIntegralMatrices::TStar = (Eigen::VectorXd(37) <<
		.1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5,
		3., 3.5, 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18., 20., 25.,
		30., 35., 40., 50., 75., 100.).finished();

	const Eigen::MatrixXd CollisionIntegralMatrices::O37_8 = (Eigen::MatrixXd(37, 8) <<
		4.008, 4.002, 4.655, 5.52, 6.454, 8.214, 9.824, 11.31,
		3.130, 3.164, 3.355, 3.721, 4.198, 5.23, 6.225, 7.160,
		2.649, 2.657, 2.77, 3.002, 3.319, 4.054, 4.785, 5.483,
		2.314, 2.32, 2.402, 2.572, 2.812, 3.386, 3.972, 4.539,
		2.066, 2.073, 2.14, 2.278, 2.472, 2.946, 3.437, 3.918,
		1.877, 1.885, 1.944, 2.06, 2.225, 2.628, 3.054, 3.747,
		1.729, 1.738, 1.79, 1.893, 2.036, 2.388, 2.763, 3.137,
		1.6122, 1.622, 1.67, 1.76, 1.886, 2.198, 2.535, 2.872,
		1.517, 1.527, 1.572, 1.653, 1.765, 2.044, 2.35, 2.657,
		1.44, 1.45, 1.49, 1.564, 1.665, 1.917, 2.196, 2.4780,
		1.3204, 1.33, 1.364, 1.425, 1.51, 1.72, 1.956, 2.199,
		1.234, 1.24, 1.272, 1.324, 1.394, 1.573, 1.777, 1.99,
		1.168, 1.176, 1.202, 1.246, 1.306, 1.46, 1.64, 1.827,
		1.1166, 1.124, 1.146, 1.185, 1.237, 1.372, 1.53, 1.7,
		1.075, 1.082, 1.102, 1.135, 1.181, 1.3, 1.441, 1.592,
		1.0006, 1.005, 1.02, 1.046, 1.08, 1.17, 1.278, 1.397,
		.95, .9538, .9656, .9852, 1.012, 1.082, 1.168, 1.265,
		.9131, .9162, .9256, .9413, .9626, 1.019, 1.09, 1.17,
		.8845, .8871, .8948, .9076, .9252, .972, 1.03, 1.098,
		.8428, .8446, .850, .859, .8716, .9053, .9483, .9984,
		.813, .8142, .8183, .825, .8344, .8598, .8927, .9316,
		.7898, .791, .794, .7993, .8066, .8265, .8526, .8836,
		.7711, .772, .7745, .7788, .7846, .8007, .822, .8474,
		.7555, .7562, .7584, .7619, .7667, .78, .7976, .8189,
		.7422, .743, .7446, .7475, .7515, .7627, .7776, .796,
		.72022, .7206, .722, .7241, .7271, .7354, .7464, .76,
		.7025, .703, .704, .7055, .7078, .7142, .7228, .7334,
		.68776, .688, .6888, .6901, .6919, .697, .704, .7125,
		.6751, .6753, .676, .677, .6785, .6827, .6884, .6955,
		.664, .6642, .6648, .6657, .6669, .6704, .6752, .681,
		.6414, .6415, .6418, .6425, .6433, .6457, .649, .653,
		.6235, .6236, .6239, .6243, .6249, .6267, .629, .632,
		.60882, .6089, .6091, .6094, .61, .6112, .613, .6154,
		.5964, .5964, .5966, .597, .5972, .5983, .600, .6017,
		.5763, .5763, .5764, .5766, .5768, .5775, .5785, .58,
		.5415, .5415, .5416, .5416, .5418, .542, .5424, .543,
		.518, .518, .5182, .5184, .5184, .5185, .5186, .5187).finished();


	const Eigen::MatrixXd CollisionIntegralMatrices::P37_8 = (Eigen::MatrixXd(37, 8) <<
		4.1, 4.266, 4.833, 5.742, 6.729, 8.624, 10.34, 11.890,
		3.263, 3.305, 3.516, 3.914, 4.433, 5.57, 6.637, 7.618,
		2.84, 2.836, 2.936, 3.168, 3.511, 4.329, 5.126, 5.874,
		2.531, 2.522, 2.586, 2.749, 3.004, 3.64, 4.282, 4.895,
		2.284, 2.277, 2.329, 2.46, 2.665, 3.187, 3.727, 4.249,
		2.084, 2.081, 2.13, 2.243, 2.417, 2.862, 3.329, 3.786,
		1.922, 1.924, 1.97, 2.072, 2.225, 2.641, 3.028, 3.435,
		1.7902, 1.795, 1.84, 1.934, 2.07, 2.417, 2.788, 3.156,
		1.682, 1.689, 1.733, 1.82, 1.944, 2.258, 2.596, 2.933,
		1.593, 1.60, 1.644, 1.725, 1.84, 2.124, 2.435, 2.746,
		1.455, 1.465, 1.504, 1.574, 1.67, 1.913, 2.181, 2.45,
		1.355, 1.365, 1.4, 1.461, 1.544, 1.754, 1.989, 2.228,
		1.28, 1.289, 1.321, 1.374, 1.447, 1.63, 1.838, 2.053,
		1.222, 1.231, 1.26, 1.306, 1.37, 1.532, 1.718, 1.912,
		1.176, 1.184, 1.209, 1.25, 1.307, 1.45, 1.618, 1.795,
		1.0933, 1.1, 1.119, 1.15, 1.193, 1.304, 1.435, 1.578,
		1.039, 1.044, 1.06, 1.083, 1.117, 1.204, 1.31, 1.428,
		.9996, 1.004, 1.016, 1.035, 1.062, 1.133, 1.22, 1.32,
		.9699, .9732, .983, .9991, 1.021, 1.08, 1.153, 1.236,
		.9268, .9291, .936, .9473, .9628, 1.005, 1.058, 1.12,
		.8962, .8979, .903, .9114, .923, .9545, .9955, 1.044,
		.8727, .8741, .878, .8845, .8935, .918, .9505, .9893,
		.8538, .8549, .858, .8632, .8703, .890, .9164, .9482,
		.8379, .8388, .8414, .8456, .8515, .868, .8895, .916,
		.8243, .8251, .8273, .8308, .8356, .8493, .8676, .89,
		.8018, .8024, .8039, .8065, .810, .820, .8337, .8504,
		.7836, .784, .7852, .7872, .7899, .7976, .808, .8212,
		.7683, .7687, .7696, .771, .7733, .7794, .788, .7983,
		.7552, .7554, .7562, .7575, .7592, .764, .771, .7797,
		.7436, .7438, .7445, .7455, .747, .7512, .757, .7642,
		.71982, .72, .7204, .7211, .7221, .725, .7289, .7339,
		.701, .7011, .7014, .702, .7026, .7047, .7076, .7112,
		.68545, .6855, .686, .686, .6867, .6883, .6905, .693,
		.6723, .6724, .6726, .673, .6733, .6745, .676, .6784,
		.651, .651, .6512, .6513, .6516, .6524, .6534, .6546,
		.614, .614, .6143, .6145, .6147, .6148, .6148, .6147,
		.5887, .5889, .5894, .59, .5903, .5901, .5895, .5885).finished();

	const double CollisionIntegralMatrices::FITASTAR_1 = 0.1106910525E+01;
	const double CollisionIntegralMatrices::FITASTAR_2 = -0.7065517161E-02;
	const double CollisionIntegralMatrices::FITASTAR_3 = -0.1671975393E-01;
	const double CollisionIntegralMatrices::FITASTAR_4 = 0.1188708609E-01;
	const double CollisionIntegralMatrices::FITASTAR_5 = 0.7569367323E-03;
	const double CollisionIntegralMatrices::FITASTAR_6 = -0.1313998345E-02;
	const double CollisionIntegralMatrices::FITASTAR_7 = 0.1720853282E-03;

	const double CollisionIntegralMatrices::FITBSTAR_1 = 0.1199673577E+01;
	const double CollisionIntegralMatrices::FITBSTAR_2 = -0.1140928763E+00;
	const double CollisionIntegralMatrices::FITBSTAR_3 = -0.2147636665E-02;
	const double CollisionIntegralMatrices::FITBSTAR_4 = 0.2512965407E-01;
	const double CollisionIntegralMatrices::FITBSTAR_5 = -0.3030372973E-02;
	const double CollisionIntegralMatrices::FITBSTAR_6 = -0.1445009039E-02;
	const double CollisionIntegralMatrices::FITBSTAR_7 = 0.2492954809E-03;

	const double CollisionIntegralMatrices::FITCSTAR_1 = 0.8386993788E+00;
	const double CollisionIntegralMatrices::FITCSTAR_2 = 0.4748325276E-01;
	const double CollisionIntegralMatrices::FITCSTAR_3 = 0.3250097527E-01;
	const double CollisionIntegralMatrices::FITCSTAR_4 = -0.1625859588E-01;
	const double CollisionIntegralMatrices::FITCSTAR_5 = -0.2260153363E-02;
	const double CollisionIntegralMatrices::FITCSTAR_6 = 0.1844922811E-02;
	const double CollisionIntegralMatrices::FITCSTAR_7 = -0.2115417788E-03;

	double CollisionIntegral11(double tjk, const double djk)
	{
		// In the future, the OpenSMOKEVector will be removed
		// For this purpose the LocateInSortedVectorFunction is required
		OpenSMOKE::OpenSMOKEVectorDouble deltaStar(8);
		for (unsigned int i = 0; i < 8; i++)
			deltaStar[i + 1] = CollisionIntegralMatrices::deltaStar(i);
		OpenSMOKE::OpenSMOKEVectorDouble TStar(37);
		for (unsigned int i = 0; i < 37; i++)
			TStar[i + 1] = CollisionIntegralMatrices::TStar(i);

		int itStar, idStar;
		double	udx21, udx321, dx31, dx32, dxx1, dxx2, x1, x2, x3, y1, y2, y3, a1, a2, a3;

		if (djk < -0.00001)
			std::cout << "Warning: Diffusivity collision integral undefined (1)" << std::endl;

		if (djk > 2.5)
			std::cout << "Warning: Diffusivity collision integral undefined (2)" << std::endl;

		if (tjk < 0.09)
		{
			//cout << "Diffusivity collision integral undefined (3)" << endl;
			tjk = 0.09;
		}

		if (tjk > 500.)
			std::cout << "Warning: Diffusivity collision integral undefined (4)" << std::endl;

		if (std::fabs(djk)>1.e-5 && tjk>75.)
			std::cout << "Warning: Diffusivity collision integral undefined (5)" << std::endl;


		// 1. TjkStar is larger than the last tabulated element
		if (tjk > TStar[36])
			return (.623 + tjk * (-.136e-2 + tjk * (.346e-5 - tjk * .343e-8)));

		// 2. deltajkStar is smaller than the first tabulated element
		if (std::fabs(djk) <= 1.e-5)
		{
			// interpolates using only the first column
			if (tjk <  TStar[2])
				itStar = 1;
			else
				itStar = TStar.LocateInSortedVector(tjk);

			// iterpolation
			x1 = TStar[itStar];
			x2 = TStar[itStar + 1];
			x3 = TStar[itStar + 2];
			udx21 = 1. / (x2 - x1);
			dx31 = x3 - x1;
			dx32 = x3 - x2;
			dxx1 = tjk - x1;
			dxx2 = tjk - x2;
			udx321 = 1. / (dx31 * dx32);

			// first interpolation with respect to TStar
			a1 = CollisionIntegralMatrices::O37_8(itStar - 1, 0);
			a2 = (CollisionIntegralMatrices::O37_8(itStar, 0) - a1) * udx21;
			a3 = (CollisionIntegralMatrices::O37_8(itStar + 1, 0) - a1 - a2 * dx31) * udx321;

			return (a1 + dxx1 * (a2 + a3 * dxx2));
		}

		// 3. General case

		// Finds the position of *tjk on the TStar vector
		if (tjk <  TStar[2])
			itStar = 1;
		else
			itStar = TStar.LocateInSortedVector(tjk);

		// Finds the position of *djk on the deltaStar vector
		if (djk < deltaStar[2])
			idStar = 1;
		else if (djk > deltaStar[7])
			idStar = 6;
		else
			idStar = deltaStar.LocateInSortedVector(djk);

		// Interpolation
		x1 = TStar[itStar];
		x2 = TStar[itStar + 1];
		x3 = TStar[itStar + 2];
		udx21 = 1. / (x2 - x1);
		dx31 = x3 - x1;
		dx32 = x3 - x2;
		dxx1 = tjk - x1;
		dxx2 = tjk - x2;
		udx321 = 1. / (dx31 * dx32);

		// first interpolation with respect to TStar
		a1 = CollisionIntegralMatrices::O37_8(itStar - 1, idStar - 1);
		a2 = (CollisionIntegralMatrices::O37_8(itStar, idStar - 1) - a1) * udx21;
		a3 = (CollisionIntegralMatrices::O37_8(itStar + 1, idStar - 1) - a1 - a2 * dx31) * udx321;
		y1 = a1 + dxx1 * (a2 + a3 * dxx2);

		// second interpolation with respect to TStar
		a1 = CollisionIntegralMatrices::O37_8(itStar - 1, idStar);
		a2 = (CollisionIntegralMatrices::O37_8(itStar, idStar) - a1) * udx21;
		a3 = (CollisionIntegralMatrices::O37_8(itStar + 1, idStar) - a1 - a2 * dx31) * udx321;
		y2 = a1 + dxx1 * (a2 + a3 * dxx2);

		// third interpolation with respect to TStar
		a1 = CollisionIntegralMatrices::O37_8(itStar - 1, idStar + 1);
		a2 = (CollisionIntegralMatrices::O37_8(itStar, idStar + 1) - a1) * udx21;
		a3 = (CollisionIntegralMatrices::O37_8(itStar + 1, idStar + 1) - a1 - a2 * dx31) * udx321;
		y3 = a1 + dxx1 * (a2 + a3 * dxx2);

		// first interpolation with respect to deltaStar
		x1 = deltaStar[idStar];
		x2 = deltaStar[idStar + 1];
		x3 = deltaStar[idStar + 2];
		udx21 = 1. / (x2 - x1);
		dx31 = x3 - x1;
		dx32 = x3 - x2;
		dxx1 = djk - x1;
		dxx2 = djk - x2;
		udx321 = 1. / (dx31 * dx32);
		a1 = y1;
		a2 = (y2 - a1) * udx21;
		a3 = (y3 - a1 - a2 * dx31) * udx321;

		return (a1 + dxx1 * (a2 + a3 * dxx2));
	}

	void Omega22k(	OpenSMOKE::OpenSMOKEVectorDouble& TkStar, 
					OpenSMOKE::OpenSMOKEVectorDouble& deltakStar,
					OpenSMOKE::OpenSMOKEVectorDouble& omega22k,
					const std::vector<std::string>& names)
	{
		// In the future, the OpenSMOKEVector will be removed
		// For this purpose the LocateInSortedVectorFunction is required
		OpenSMOKE::OpenSMOKEVectorDouble deltaStar(8);
		for (unsigned int i = 0; i < 8; i++)
			deltaStar[i + 1] = CollisionIntegralMatrices::deltaStar(i);
		OpenSMOKE::OpenSMOKEVectorDouble TStar(37);
		for (unsigned int i = 0; i < 37; i++)
			TStar[i + 1] = CollisionIntegralMatrices::TStar(i);

		int itStar;
		int idStar;

		double udx21, dx31, dx32, udx321, dxx1, dxx2, x1, x2, x3, y1, y2, y3, a1, a2, a3;
		double *tk = TkStar.GetHandle();
		double *dk = deltakStar.GetHandle();
		double *ok = omega22k.GetHandle();

		for (int k = 1; k <= TkStar.Size(); k++)
		{
			if (*dk < -0.00001)
				std::cout << "Warning: Viscosity-Conductivity collision integral undefined (1) - Species: " << names[k] << std::endl;

			if (*dk > 2.5)
				std::cout << "Warning: Viscosity-Conductivity collision integral undefined (2) - Species: " << names[k] << std::endl;

			if (*tk < 0.09)
			{
				std::cout << "Warning: Viscosity-Conductivity collision integral undefined (3) - Species: " << names[k] << std::endl;
				*tk = 1.00;		// be careful
			}
			if (*tk > 500.)
				std::cout << "Warning: Viscosity-Conductivity collision integral undefined (4) - Species: " << names[k] << std::endl;

			if (std::fabs(*dk) > 1.e-5 && *tk > 75.)
				std::cout << "Warning: Viscosity-Conductivity collision integral undefined (5) - Species: " << names[k] << std::endl;

			if (*tk > TStar[36])
			{
				*ok = .703 + *tk * (-.146e-2 + *tk * (.357e-5 - *tk * .343e-8));
				tk++;
				dk++;
				ok++;
				continue;
			}
			if (std::fabs(*dk) <= 1.e-5)
			{
				// Interpolates using only the first column
				if (*tk <  TStar[2])
					itStar = 1;
				else
					itStar = TStar.LocateInSortedVector(*tk);

				// Interpolation
				x1 = TStar[itStar];
				x2 = TStar[itStar + 1];
				x3 = TStar[itStar + 2];
				udx21 = 1. / (x2 - x1);
				dx31 = x3 - x1;
				dx32 = x3 - x2;
				dxx1 = *tk - x1;
				dxx2 = *tk - x2;
				udx321 = 1. / (dx31 * dx32);

				// First interpolation with respect to Tstar
				a1 = CollisionIntegralMatrices::P37_8(itStar - 1, 0);
				a2 = (CollisionIntegralMatrices::P37_8(itStar, 0) - a1) * udx21;
				a3 = (CollisionIntegralMatrices::P37_8(itStar + 1, 0) - a1 - a2 * dx31) * udx321;
				*ok = a1 + dxx1 * (a2 + a3 * dxx2);
				tk++;
				dk++;
				ok++;
				continue;
			}

			// Finds the position of *tjk on the TStar vector
			if (*tk <  TStar[2])
				itStar = 1;
			else
				itStar = TStar.LocateInSortedVector(*tk);

			// Finds the position of *djk on the deltaStar vector
			if (*dk < deltaStar[2])
				idStar = 1;
			else if (*dk > deltaStar[7])
				idStar = 6;
			else
				idStar = deltaStar.LocateInSortedVector(*dk);

			// Interpolation
			x1 = TStar[itStar];
			x2 = TStar[itStar + 1];
			x3 = TStar[itStar + 2];
			udx21 = 1. / (x2 - x1);
			dx31 = x3 - x1;
			dx32 = x3 - x2;
			dxx1 = *tk - x1;
			dxx2 = *tk - x2;
			udx321 = 1. / (dx31 * dx32);

			// First interpolation with respect to Tstar
			a1 = CollisionIntegralMatrices::P37_8(itStar - 1, idStar - 1);
			a2 = (CollisionIntegralMatrices::P37_8(itStar, idStar - 1) - a1) * udx21;
			a3 = (CollisionIntegralMatrices::P37_8(itStar + 1, idStar - 1) - a1 - a2 * dx31) * udx321;
			y1 = a1 + dxx1 * (a2 + a3 * dxx2);

			// Second interpolation with respect to Tstar
			a1 = CollisionIntegralMatrices::P37_8(itStar - 1, idStar);
			a2 = (CollisionIntegralMatrices::P37_8(itStar, idStar) - a1) * udx21;
			a3 = (CollisionIntegralMatrices::P37_8(itStar + 1, idStar) - a1 - a2 * dx31) * udx321;
			y2 = a1 + dxx1 * (a2 + a3 * dxx2);

			// Third interpolation with respect to Tstar
			a1 = CollisionIntegralMatrices::P37_8(itStar - 1, idStar + 1);
			a2 = (CollisionIntegralMatrices::P37_8(itStar, idStar + 1) - a1) * udx21;
			a3 = (CollisionIntegralMatrices::P37_8(itStar + 1, idStar + 1) - a1 - a2 * dx31) * udx321;
			y3 = a1 + dxx1 * (a2 + a3 * dxx2);

			// // First interpolation with respect to deltaStar
			x1 = deltaStar[idStar];
			x2 = deltaStar[idStar + 1];
			x3 = deltaStar[idStar + 2];
			udx21 = 1. / (x2 - x1);
			dx31 = x3 - x1;
			dx32 = x3 - x2;
			dxx1 = *dk - x1;
			dxx2 = *dk - x2;
			udx321 = 1. / (dx31 * dx32);
			a1 = y1;
			a2 = (y2 - a1) * udx21;
			a3 = (y3 - a1 - a2 * dx31) * udx321;
			*ok = a1 + dxx1 * (a2 + a3 * dxx2);
			tk++;
			dk++;
			ok++;
		}
	}
}

#endif