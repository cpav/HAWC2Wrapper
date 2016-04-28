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