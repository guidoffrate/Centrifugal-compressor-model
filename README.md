# Centrifugal-compressor-model
This is a 1-D Matlab model of a centrifugal compressor (impeller + vaneless diffuser) that simulates the compressor performance given the machine's geometrical specifications.
The model is thought explicitly to allow you to:
- optimize the compressor geometry to maximize the total-to-static isentropic efficiency given a prescribed total-to-static pressure ratio (single objective optimization);
- optimize the compressor geometry to simultaneously maximize the total-to-static isentropic efficiency and the total-to-static pressure ratio (multi-objective optimization);
- generate the machine performance maps by simulating off-design operation once the geometry is fixed;

# Content
The repository includes the following files:
- _centrifugal_compressor.m_ simulates the compressor given some machine's geometrical specifications, including blade angles, radiuses, blade passage height, and volute and exhaust cone dimensions;
- _compressor_validation.m_ allows you to visualize the comparison between model predictions and the often-used Eckardt's compressor experimental data, which refers to three different impeller geometries;
- _compressor_optimisation_test.m_ runs two optimizations, a single- and a multi-objective one, allowing you to compare the results between the two approaches. This file requires the presence of the experimental datasets, included in the repositories in the .mat files named _eckardt_impeller_A/B/O.mat_;
- _compressor_map_test.m_ simulates the performance of a compressor geometric design specified in the file for several combinations of rotating speed and mass flow rate. Finally, the file plots a compressor performance maps, including the total-to-static pressure ratio and total-to-static isentropic efficiency.

# Citation
The explanation of the model equations and hypotheses, as well as the validation against experimental data from the literature, are reported in the paper:

- Frate, G.F., Benvenuti, M., Chini, F., and Ferrari, L. (2024). _Optimal design of centrifugal compressors in small-size high-temperature Brayton heat pumps_, Proceedings of 37th International Conference on
Efficiency, Cost, Optimization, Simulation and Environmental Impact of Energy Systems (ECOS), Rodhes, Greece, 30 June - 5 July 2024, doi: 10.52202/077185-0031

Please cite this paper if you use the compressor model in your publications!




