# sootDNSFoam 

The solver is used for sooting or non-sooting laminar flames, which is based on OpenFOAM-7 and Cantera 2.4. 

The thermodynamic properties and reaction kinetics are calculated using Cantera with detailed chemical reaction mechanisms (including quasi-steady state mechanisms). The mixture-averaged approach is used to calculate the molecular diffusion of transported species. Soot particle dynamics is handled using the hybrid method of moments (HMOM).

We are continuously updating the solver, currently we focus on the soot model and the chemical reaction solver. Please read the documentation to install, compile and use.

It was developed by Dr. Junjun Guo (junjun.guo@kaust.edu.sa; King Abdullah University of Science and Technology) and Dr. Yihao Tang (yhtang@umich.edu; University of Michigan).

Publication related to this solver:

[1] Junjun Guo, Prabhu Selvaraj, Yihao Tang, Hong G. Im, and Venkatramanan Raman. "An analysis of soot formation pathways in laminar coflow ethylene flame at higher pressures." In AIAA Scitech 2020 Forum, p. 1660. 2020. https://doi.org/10.2514/6.2020-1660

[2] Junjun Guo, Yihao Tang, Venkat Raman, and Hong G. Im. "Numerical investigation of pressure effects on soot formation in laminar coflow ethylene/air diffusion flames." Fuel 292 (2021): 120176. https://doi.org/10.1016/j.fuel.2021.120176

[3] Hanfeng Jin, Junjun Guo, Tianyu Li, Zhongyue Zhou, Hong G. Im, and Aamir Farooq. "Experimental and numerical study of polycyclic aromatic hydrocarbon formation in ethylene laminar co-flow diffusion flames." Fuel 289 (2021): 119931. https://doi.org/10.1016/j.fuel.2020.119931

[4] Junjun Guo, Peng Liu, Erica Quadarella, Kiran Yalamanchi, Ibraheem Alsheikh, Carson Chu, Fengshan Liu, S. Mani Sarathy, William L. Roberts, and Hong G. Im. "Assessment of physical soot inception model in normal and inverse laminar diffusion flames." Combustion and Flame 246 (2022): 112420. https://doi.org/10.1016/j.combustflame.2022.112420

[5] Peng Liu, JunJun Guo, Erica Quadarella, Anthony Bennett, Sreenivasa R. Gubba, Saumitra Saxena, Obulesu Chatakonda et al. "The effect of preheating temperature on PAH/soot formation in methane/air co-flow flames at elevated pressure." Fuel 313 (2022): 122656. https://doi.org/10.1016/j.fuel.2021.122656

[6] Peng Liu, Junjun Guo, Hong G. Im, and William L. Roberts. "The effects of CO2/CH4 ratio on soot formation for autothermal reforming of methane at elevated pressure." Combustion and Flame (2022): 112379. https://doi.org/10.1016/j.combustflame.2022.112379
