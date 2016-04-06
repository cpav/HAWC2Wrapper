===========
Usage Guide
===========

The wrapper is devided in the following main blocks:

* Wind turbine **variable trees** (``hawc2_vartrees.py`` and ``vartrees.py``); a definition of a class, composed by several hinerited classes, that defines the turbine;
* A **reader** (``hawc2_inputreader.py``, ``hawc2_inputdict.py``), that reads a reference master file and store all the data;
* A **writer** (``hawc2_inputwriter.py``), that writes the data of a turbine into an htc file;
* An **executer** (``hawc2_wrapper.py``), that executes HAWC2 or HAWC2s for a given htc file;
* A **postprocessor** (``hawc2_output.py``), that reads result files, performs post processing, and stores them;
* A **geometry builder** (``hawc2_geometry.py``), performs some geometrical operationis on the blade geometry;
* An OpenMDAO **workflow** (``hawc2_aeroelasticsolver.py`` and ``hawc2s_aeroelasticsolver.py``), that performs the entire worflow.

The following sections describes each of these blocks.

Variable tree
-------------

The wind turbine is defined in the wrapper through a variable tree.
The variable tree is defined from the class ``HAWC2VarTrees`` and it has the following structure::

    attr: sim(HAWC2Simulation())
        attr: time_stop, solvertype, convergence_limits, on_no_convergence, max_iterations, newmark_deltat, eig_out, logfile
    
    attr: wind(HAWC2Wind())
        attr: sdensity, wsp, tint, horizontal_input, center_pos0, windfield_rotations, shear_type, shear_factor, turb_format, tower_shadow_method, scale_time_start, wind_ramp_t0, wind_ramp_t1, wind_ramp_factor0, wind_ramp_factor1, wind_ramp_abs, iec_gust, iec_gust_type, G_A, G_phi0 ,G_t0, G_T
    
    attr: aero(HAWC2Aero())
        attr: nblades, hub_vec_mbdy_name, hub_vec_coo, links, induction_method, aerocalc_method, aerosections, tiploss_method, dynstall_method, ae_sets, ae_filename, pc_filename, aero_distribution_file, aero_distribution_set

    attr: aerodrag = HAWC2AeroDrag()
        attr: elements (list of HAWC2AeroDragElement)
            attr: mbdy_name, dist, sections, calculation_points

    attr: blade_ae(HAWC2BladeGeometry())
        attr: radius , s, chord, rthick, aeset

    attr: airfoildata(HAWC2AirfoilData())
        attr: nset, desc, pc_sets

    attr: output(HAWC2OutputVT())
        attr: filename, time_start, time_stop, out_format, out_buffer, res_dir

    attr: main_bodies(HAWC2MainBodyList())
        attr: body_name(HAWC2MainBody()) (body_name depends on the name of the body)
            attr: body_name, body_type, st_filename, st_input_type, beam_structure, body_set, nbodies, node_distribution, damping_type, damping_posdef, damping_aniso, copy_main_body, c12axis, concentrated_mass
            attr: orientations (list of HAWC2OrientationBase or HAWC2OrientationRelative)
                HAWC2OrientationBase attr: type, body, inipos, body_eulerang
                HAWC2OrientationRelative attr: type, body1, body2, body2_eulerang, mbdy2_ini_rotvec_d1
            var: constraints (list of classes)
                class: HAWC2ConstraintFix0
                    attr: con_type, mbdy, disable_at
                class: HAWC2ConstraintFix1
                    attr: con_type, mbdy1, mbdy2, disable_at
                class: HAWC2ConstraintFix23
                    attr: con_type, mbdy, dof
                class: HAWC2ConstraintFix4
                    attr: con_type, mbdy1, mbdy2, time
                class: HAWC2ConstraintBearing12
                    attr: name, con_type, mbdy1, mbdy2, bearing_vector, disable_at
                class: HAWC2ConstraintBearing3
                    attr: name, con_type, mbdy1, mbdy2, bearing_vector, omegas
                class: HAWC2ConstraintBearing45
                    attr: name, con_type, mbdy1, mbdy2, bearing_vector
    
    attr: dlls(HAWC2Type2DLLList())
        attr: dll_name(HAWC2Type2DLL()) (dll_name depends on the name of the DLL)
            attr:  name, filename, dll_subroutine_init, dll_subroutine_update, arraysizes_init, arraysizes_update, deltat, output
            attr: output(HAWC2Type2DLLIO())
                attr: out_dic, action_dic
            attr: actions(HAWC2Type2DLLIO())
                attr: out_dic, action_dic
            attr: dll_init(DTUBasicControllerVT() or HAWC2Type2DLLinit())
                DTUBasicControllerVTI(HAWC2Type2DLLinit()) 
                    attr: Vin, Vout, nV, ratedPower = 0 minRPM, maxRPM, gearRatio, designTSR, active, FixedPitch, maxTorque, minPitch, maxPitch, maxPitchSpeed, maxPitchAcc, generatorFreq, generatorDamping, ffFreq, Qg, pgTorque, igTorque, dgTorque, pgPitch, igPitch, dgPitch, prPowerGain, intPowerGain, generatorSwitch, KK1, KK2, nlGainSpeed, softDelay, cutin_t0, stop_t0, TorqCutOff, stop_type, PitchDelay1, PitchVel1, PitchDelay2, PitchVel2, generatorEfficiency, overspeed_limit, minServoPitch, maxServoPitchSpeed, maxServoPitchAcc, poleFreqTorque, poleDampTorque, poleFreqPitch, .poleDampPitch, gainScheduling, prvs_turbine, rotorspeed_gs, Kp2, Ko1, Ko2
                HAWC2Type2DLLinit attr:
                    attr: init_dic

    attr: h2s(HAWC2SVar())
        attr: cases, wsp_cases, wsp_curve, pitch_curve, rpm_curve, operational_data_filename, commands
        attr: ground_fixed(HAWC2SBody()), rotating_axissym(HAWC2SBody()), rotating_threebladed(HAWC2SBody())
            attr: main_body, log_decrements
        attr: second_order_actuator(SecondOrderActuator())
            attr: name, frequency, damping
        attr: options(HAWC2SCommandsOpt())
            attr: include_torsiondeform, bladedeform, tipcorrect, induction, gradients, blade_only, matrixwriteout, eigenvaluewriteout, frequencysorting, number_of_modes, maximum_damping, minimum_frequency, zero_pole_threshold, aero_deflect_ratio, vloc_out, regions, set_torque_limit
        attr: ch_list_in(HAWC2OutputListVT()), ch_list_out(HAWC2OutputListVT())
            attr: sensor_list

    attr: body_order, dlls_order


Reader
------
The module ``hawc2_inputreader`` contains the class ``HAWC2InputReader``. This class can read both HAWC2 and HAWCStab2 htc files and creates a corresponding variable tree.

This class calls internally the class ``HAWC2InputDict`` contained in the module ``hawc2_inputdict``. This second class is the one that actually reads the htc file and it stores it into a dictionary. The class ``HAWC2InputReader`` interprets the dictionary and creates the varibale tree.

An example of the use of the class is::
    
    >>> from hawc2_inputreader import HAWC2InputReader
    
    >>> Reader = HAWC2InputReader('hawc2_master.htc')
    >>> Reader.execute()
    
After the execution the object *Reader* contains and attribute *vartrees* that contains the variable tree describing the turbine.

The module ``hawc2_inputdict`` includes functions to read *st*, *pc*, and *ae* files.

Writer
------
The module ``hawc2_inputwriter`` contains two classes: ``HAWC2InputWriter`` and ``HAWC2SInputWriter``. The class ``HAWC2SInputWriter`` is of type ``HAWC2InputWriter`` therefore it inherits all the attributes of ``HAWC2InputWriter``.
The first class is used to write HAWC2 input files, the second HAWC2s files.

During the initialization of the classes the following attributes can be specified:
    * case_id: for the name of the htc file to be created;
    * vartrees(HAWC2VarTrees()): to initialize the variable that needs to be converted into an htc file;
    * data_directory: for the opath of the data directory;
    * res_directory: for the path of the results directory;
    * turb_directory: for the path of the turbulence files directory;
    * log_directory: for the path of the log files directory;
    * control_directory: for the path of the controller directory.

The class, once executed does not return anything, it only creates a new htc file.

An example of the use of the classes is::

    >>> from hawc2_inputwriter import HAWC2InputWriter, HAWC2SInputWriter
    >>> writer = HAWC2InputWriter()
    >>> writer.case_id = 'new_file_h2'
    >>> writer.vartrees = reader.vartrees
    >>> writer.execute()
    >>> writer = HAWC2SInputWriter()
    >>> writer.case_id = 'new_file_h2s'
    >>> writer.vartrees = reader_h2s.vartrees
    >>> writer.execute()

Executer
--------
The module ``hawc2_wrapper.py`` is the one that is actually launching the runs of HAWC2 and HAWCStab2. The module includes only the class ``HAWC2Wrapper``.  The class performs some cheks of the log files printing on the screen the errors that are found in the logfile. The class also copies the result directory.

An example is::
    
    >>> from hawc2_wrapper import HAWC2Wrapper
    >>> executer = HAWC2Wrapper()
    >>> executer.hawc2bin = 'HAWC2s.exe'
    >>> executer.case_id = 'new_file_h2s'
    >>> executer.execute()

Postprocessor
-------------
The module ``hawc2_output`` is the responsible for reading the output files and perform postprocessing of the data. The module depends on the Wind Energy Toolbox (wetb_) that is used for reading the HAWC2 results files, calculating statistics, fatigue, and load envelopes.

The module containts several classes to perform different types of postprocessing or select and return only specific channels and or results.

The classes are:

    * ``HAWC2OutputBase``: reads HAWC2 output files and computes statistics, fatigue and envelope;
    * ``HAWC2Output(HAWC2OutputBase())``: reorganizes the results read by HAWC2OutputBase into arrays;
    * ``HAWC2SOutputBase``: reads HAWC2s output files;
    * ``HAWC2SOutput(HAWC2SOutputBase())``: reorganizes the results in arrays;
    * ``HAWC2SOutputCompact(HAWC2SOutput())``: reorganizes the results in two compacts arrays;
    * ``FreqDampTarget``: computes a cost funciton for freqeuncy and damping placement;

.. _wetb: https://gitlab.windenergy.dtu.dk/toolbox/WindEnergyToolbox/

The classes can be used as follow::

    >>> output = HAWC2SOutputBase()
    >>> output.case_id = wrapper.case_id
    >>> output.commands = writer.vartrees.h2s.commands
    >>> output.execute()

    >>> case = {}
    >>> case['[case_id]'] = 'dlc12_wsp04_wdir000_s1001'
    >>> case['[res_dir]'] = 'res/dlc12_iec61400-1ed3'
    >>> config = {}
    >>> config['neq'] = 600
    >>> config['no_bins'] = 2**7
    >>> config['m'] = [12]
    >>> output = HAWC2OutputBase(config)
    >>> output.execute(case)

Geometry Builder
----------------
This module contains the class ``HAWC2GeometryBuilder`` and depends on the library PGL_.
The class can be used mainly for two applications:
    * interpolate the initialized c2def values to change the spanwise discretization;
    * define the c2def from the variable tree BladeGeometryVT()
    
.. _PGL: https://gitlab.windenergy.dtu.dk/frza/PGL 

Workflow
--------
The modules ``hawc2_aeroelasticsolver`` and ``hawc2s_aeroelasticsolver`` implement workflows to initialize the variable trees, modify the variable trees according to some inputs, write the htc file, execute HAWC2 or HAWC2s, and perform the postprocessing. The workflows are implemented in the calsses ``HAWC2Workflow`` and ``HAWC2sWorkflow``. These classes are ``Component`` classes. In these classes the variables that need to be changed in the variable trees with respect to the initialization variable tree are added as parameters. Three groups of variables can be added: one for the tip speed ratio, one for the structural properties of the balde, and one for the blade geometry.

The modules contain a ``Group`` class (from OpenMDAO) each, called ``HAWC2AeroElasticSolver`` and ``HAWC2sAeroElasticSolver``. These classes execute the workflow with a ``ParallelGroup``. 
Each of the module also include a class ``Component`` called ``OutputAggregator`` that groups the outputs from the parallel computation into single arrays.

