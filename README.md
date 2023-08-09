# QERaman
An open-source program for computing the first-order resonance Raman spectroscopy.
![logo](https://github.com/nguyen-group/QERaman/assets/46996256/9f8f7137-03f1-435a-8e77-3de463bb7afa)

# Requirement
QERaman requires Quantum ESPRESSO version >= 7.0

# Installation
Step 1: Download the latest version of QERaman inside the main Quantum ESPRESSO directory:

    git clone https://github.com/nguyen-group/QERaman.git

Step 2: Go to the source code in the QERaman directory to install the code:

    cd QERaman/src
    make all

It noted that the reader should `make pw pp ph` in the Quantum ESPRESSO directory before `make all` in the QERaman directory. Also, consider adding the `QERaman/bin` to the environment variables `PATH`.

# Documentation
The documentation can be found on the GitHub wiki page: https://github.com/nguyen-group/QERaman/wiki.
It shows the input variables description for `raman.x`, `bands_mat.x`, and `ph_mat.x`, and also the format of the output files.

# References and citing
The theory behind the QERaman code is described in our pre-print:
> QERaman: A open-source program for calculating resonance Raman spectroscopy based on Quantum ESPRESSO,
>
> Nguyen T. Hung, Jianqi Huang, Teng Yang, Riichiro Saito,
>
> https://arxiv.org/abs/xxxx

The implementation of QERaman is built on the Quantum ESPRESSO package and the detail of the input variables can be found in the Quantum ESPRESSO book. If you use the QERaman code in your work, please consider citing the Quantum ESPRESSO package and book as well: 
> https://doi.org/10.1088/0953-8984/21/39/395502 (for Quantum ESPRESSO package)
>
> https://doi.org/10.1201/9781003290964 (for Quantum ESPRESSO book)

# Contributors
- Nguyen Tuan Hung (Tohoku University, Japan)
- Jianqi Huang (Institute of Metal Research, Chinese Academy of Sciences, China)
- Teng Yang (Institute of Metal Research, Chinese Academy of Sciences, China)
- Riichiro Saito (Tohoku University, Japan)

We are open to contributions from anybody who wants to develop the QERaman. Don't hesitate to contact us using the GitHub discussions page: https://github.com/nguyen-group/QERaman/discussions to discuss your ideas.

# License
GNU General Public License (v3)
