=====================
HAWC2Wrapper Examples
=====================

In the test folder some examples are present.

Example.py
----------

The file *Example.py* contanis two simple optimization examples. These two exaples show how to set-up a workflow that exploits the HAWC2wrapper functionalities to identify a blade stiffness distribution to obtain a specific blade freqeuncy and tip displacement. To switch between the two you can comment and uncomment lines 213 and 214.

The example has different functions to perform specific operations:
    * setup_config_dict: to setup the configuration dictionary that is required by the wrapper;
    * setup_init_structure: to read the initial structure that is going to be used to normalize the stiffnesses;
    * set_top: initialize the problem and select the optimizer;
    * example_blade_tip_target: core part for the blade tip deflection example where design variables and cost functions are added to the problem. Here the number of blade sections used as design variable is specified with the variable n;
    * example_frequency_placement: core part for the freqeuncy placement example where design variables and cost functions are added to the problem. Here the number of blade sections used as design variable is specified with the variable n;

Three components are also defined:
    * NormalizeDesVar: it normalizes the design variable, in this case the blade stiffnesses;
    * CostFunction: it defines the cost function of the tip displacement problem;
    * CostFunction2: it defines the cost function of the freqeuncy placement problem.
    
These examples can be executed also by computing the finite differences in parallel.
To do so you need to change the variable *par_fd* from 1 to a higher number that is lower than the number of design variables *n*.

test_hawc2_wrapper_openmdao.py
------------------------------

This file tests a workflow to perform both HAWC2 and HAWCStab2 computations. It is possible to switch between the two software by changing the variable *H2* at line 17. 

The HAWCStab2 workflow computes the stady states at the wind speeds defined in the configuration dictionary. Both a user defined case, where wind speed, pitch angle and rotor speed are provided, and a general case where a wind sweep is perormed are evaluated. The freqeuncy placement module is also called to evaluate the distance of a modal freqeuncy from a target value.

The HAWC2 workflow launches a set of simulations defined in the folder *DLCs_long* and computes statistics and fatigue.

test_hawc2s_w_fatigue_openmdao.py
---------------------------------

This example is definetely not optimal but it's what I have so far..

This file tests the workflow to compute the fatigue with HAWCStab2 in frequency domain. To run this workflow some steps need to be performed beforehand to generate the wind input files.

    * Uncomment line 69 `compute_wind_input()` and run the script. This will run a function that calls HAWC2 for the cases defined in *DLC2_fatigue*;
    * Once the simulations are executed, it will crash but in your working folder there will be 3 new folders with HAWC2 results;
    * compy the name of the HAWC2 results folder in the variable *cases_list*. Note that you should only change the 8 digit number a the end of the main folder name;
    * comment line 69 `compute_wind_input()` and uncomment line 75 `fatfreq.ComputeWindResp(cases_list, 'wind_structure', Fmax=2.5)`
    * run the script.
    

