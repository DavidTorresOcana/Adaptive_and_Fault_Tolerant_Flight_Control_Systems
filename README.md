# Adaptive-and-Fault-Tolerant-Flight-Control-Systems

## Synopsis

For flight control systems, this paper proposes an adaptive control approach based on a framework of Explicit Model Following Direct Adaptive Control scheme. As a first step, a modified F-16 dynamics model is developed to explore control surface redundancies, as well as to enable modelling of dynamics changes result from faults, failures and/or plant deviations. In this modified model, each control surface can be individually controlled. Next, this paper proposes a flight control framework that integrates an Adaptive Neural Network, non-linear dynamic inversion, control allocation, System Identification with Uncested Kalman Filter and Model Reference Following scheme to leverage their synergies. Then, the proposed approach is tested using the F-16 nonlinear model developed and its performance is validated via numerical simulations.

## Docs and publications 
Thesis Doc explain all the theory behind this repository. It can be accessed through:
>   [Thesis PDF](https://cran.ent.sirsidynix.net.uk/client/en_GB/knl/search/detailnonmodal/ent:$002f$002fSD_ILS$002f0$002fSD_ILS:365548/ada)

For more help with the theory please [contact here](mailto:david.torres.ocana@gmail.com).

See here IFAC submition, which explain the the plant model and Adaptive Neural Networks theory and implementation:

>   [IFAC publication PDF](http://www.sciencedirect.com/science/article/pii/S2405896315009404)


## Motivation

Current linear Flight Control Systems (FCS) algorithms are incapable of adapting to sudden changes in terms of aircraft configuration. It is well known that classical control approaches only provide a satisfying performance and robustness if the aircraft is close enough to the model assumed for control design. Any uncertainties or failures lead to degradations in stability and performance. Therefore, linear, model-based control techniques might require complete redesign of control if there are significant changes in aircraft configuration. As a result this tends to restrict the ability to alter the design or carry new equipment or to handle in-flight reconfigurations.

## Requirements

* Requires Matlab 2015a (currently running) or above (some changes may be needed)
* Requires Matlab mex compiler to be installed. See here https://www.mathworks.com/support/sysreq/previous_releases.html
* Requires Flight Gear 2.6 or above
* Requires gamepad or equivalent USB compatible Joystick. Tested with a Logitech G F310

## Installation and execution
### FlightGear (optional)

* Install [FlightGear](https://www.flightgear.org/download/)   (FG). Version tested to work was 2.6 , but other versions may work as well
* Install one of the [models in here](/tools/F16 Model/Fligth Gear/) into FlighGear. *f16_20120812.zip* is the recomended:
    - Uncompress *f16_20120812.zip*, copy ```f16``` folder onto ```C:\Program Files\FlightGear\data\Aircraft``` or wherever you installed FG in your system
* In *f16* folder, open *f16-set.xml* file and modify ```<flight-model>``` field to be ```network```. This will make FG listen ```localhost``` for simulator input
* On the Simulink model (see next section):
    * Configure the block on the image to match your FG configuration
    ![Configure](bat_gen.png)
    * Generate run script (see inside that block)
    * Make sure that block and ```Send net_fdm Packet to FlightGear``` block have same configuration
    * Before executing Simulink run FlightGear by running the generated *.bat* file. You can leave FG running independently of Simulink
    
### FCS

* Add to the Matlab path all libraries in ```lib/```
* Navigate to ```src/7dof FCS Development```
* Run ```RUN_ME.m``` with Matlab and select a flight condition to start with
* A Simulink model will open: press play to start the simulation

In order to have a *real experience*, you could head to *Cockpit* block on the Simulink model and select to use a Joystick as Cockpit inputs. Most gamepads are compatible. 

## Tests

Using the GUI, stop or run the simulation and inject failures, deviations or variations of the F16 model

## Demo Videos

### Tests and demo
[![Failures injections and behaviours](http://img.youtube.com/vi/WFtM5UVIlB4/mqdefault.jpg)](https://www.youtube.com/watch?v=WFtM5UVIlB4)

## Contributors

@David Torres Ocaña
    david.torres.ocana@gmail.com
