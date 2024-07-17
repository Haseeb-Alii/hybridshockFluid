# Low-dissipation hybrid central solver for LES of supersonic turbulent flow.
## Release date: 17 July 2024
## Overview:

This repository contains the OpenFOAM implementation of a hybrid central scheme originally implemented by Lee et al. [1][2], combining linear flux and KNP or KT flux for convective terms based on the modified Bhagatwala and Lele shock sensor [3]. The hybrid central scheme assigns weights to each flux depending on the shock location, using the KNP or KT flux in shock regions and the linear flux in expansion and other regions. More details can be found in the original paper by Lee et al. [1]. This code implements a hybrid central scheme for discretizing the convective terms in the full Navier-Stokes equations for supersonic compressible flows. The code utilizes a modified shock sensor by Bhagatwala and Lele [2] to determine the appropriate flux to use in different regions, ensuring accuracy and stability in the presence of shocks. Also, original solver "shockFluid" has also been modified to handle multiple species by implementing the species transport equation. The current "hybridshockFluid" solver has been developed and tested for OpenFOAM 11.

## Installation

#### 1. Clone the repository
git clone https://github.com/Haseeb-Alii/hybridshockFluid.git

#### 2. Navigate into the cloned repository
cd hybridshockFluid

#### 3. Source OpenFOAM

#### 4. Compile it
./Allwmake

## Contributing

Contributions are welcome! Please follow these steps to contribute:

- Fork the repository.
- Create a new branch: `git checkout -b feature-name`.
- Make your changes and commit them: `git commit -m 'Add new feature'`.
- Push to the branch: `git push origin feature-name`.
- Submit a pull request.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

#### References:

- [1]. Lee, Yachao & Yao, Wei & Fan, Xuejun. (2018). Low-Dissipative Hybrid Compressible Solver Designed for Large-Eddy Simulation of Supersonic Turbulent Flows. AIAA Journal. 56. 1-11. 10.2514/1.J056404. 
- [2]. Yachao, Lee & Yao, Wei & Fan, Xuejun. (2017). A low-dissipation scheme based on OpenFOAM designed for large eddy simulation in compressible flow. 10.2514/6.2017-2444. 
- [3]. Bhagatwala, Ankit & Lele, Sanjiva. (2009). A modified artificial viscosity approach for compressible turbulence simulations. J. Comput. Physics. 228. 4965-4969. 10.1016/j.jcp.2009.04.009. 

