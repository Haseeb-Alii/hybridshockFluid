Hybrid Central Scheme for Convective Terms in Fluid Dynamics

This repository contains the implementation of a hybrid central scheme combining linear flux and KNP or KT flux based on the modified Bagatwala and Lele shock sensor. The scheme assigns weights to each flux depending on the shock location, using the KNP or KT flux in shock regions and linear flux in expansion and other regions.
Table of Contents

    Introduction
    Shock Sensor
    Hybrid Flux
    Installation
    Usage
    Examples
    Contributing
    License
    Acknowledgments

Introduction

This project implements a hybrid central scheme for solving convective terms in fluid dynamics. The scheme utilizes a modified shock sensor by Bagatwala and Lele to determine the appropriate flux to use in different regions, ensuring accuracy and stability in the presence of shocks.
Shock Sensor

The shock sensor by Bagatwala and Lele is defined as:

Θ=12(1−tanh⁡(2.5+10Δc∣∇⋅U∣))×(∇⋅U)2(∇⋅U)2+(∇×U)2+ϵ
Θ=21​(1−tanh(2.5+10∣∇⋅U∣Δc​))×(∇⋅U)2+(∇×U)2+ϵ(∇⋅U)2​

Where:

    ΘΘ is the shock sensor value.
    ΔcΔc is the change in a given quantity cc.
    ∇⋅U∇⋅U and ∇×U∇×U are the divergence and curl of the velocity field UU, respectively.
    ϵϵ is a small positive constant to prevent division by zero.

Hybrid Flux

The hybrid flux for the convective term is given by:

∫V∇⋅[ρUΦ] dV≈ΘfFd+(1−Θf)Fc
∫V​∇⋅[ρUΦ]dV≈Θf​Fd​+(1−Θf​)Fc​

Where:

    FdFd​ is the diffusive flux (KNP or KT flux).
    FcFc​ is the central flux (linear flux).
    ΘfΘf​ is the flux weight determined by the shock sensor.

Installation

Clone the repository:

sh

git clone git@github.com:your-username/your-repository.git
cd your-repository

Install the necessary dependencies:

sh

# Example for Python projects
pip install -r requirements.txt

Usage

Provide instructions on how to use the code. For example:

sh

# Example command to run the code
python main.py --input data/input_file.txt --output results/output_file.txt

Examples

Provide examples to help users understand how to use your project. Include code snippets and expected outputs.

sh

# Example of running a test case
python main.py --input data/test_input.txt --output results/test_output.txt

Contributing

Contributions are welcome! Please follow these steps to contribute:

    Fork the repository.
    Create a new branch: git checkout -b feature-name.
    Make your changes and commit them: git commit -m 'Add new feature'.
    Push to the branch: git push origin feature-name.
    Submit a pull request.

Please ensure your code follows the project's coding standards and passes all tests.
License

This project is licensed under the MIT License. See the LICENSE file for details.
Acknowledgments

    Thanks to Bagatwala and Lele for the shock sensor methodology.
    Thanks to all contributors and users.
