# Package_ClusterModel

## Overview
This Java package provides core utilities for implementing various Individual-Based Models (IBM), particularly those involving cluster-based simulations. It supports both network generation and simulation execution, serving as the backbone for multiple IBM implementations.

IBM simulations built with this package typically follow two main steps: [Network generation](#1_network_generation) and [Simulation](#2_simulation). 

## Maintainer and developers
* Ben Hui; ORCiD ID: [0000-0002-6567-5821](https://orcid.org/0000-0002-6567-5821)


## 1. Network Generation
A network defines partnerships within a population and is represented as a CSV file with each partnership specified in the following format:
```
person_id_1, person_id_2, partnership_start, partnership_duration
```
* **person_id_1**, **person_id_2**: Unique identifiers for individuals.
* **partnership_start**: Day the partnership begins.
* **partnership_duration**: Duration of the partnership in days.
  
>[!Note]
>Network generation is independent of the simulation process. Users may generate networks using any method, as long as the format above is respected.

The package included the _Simulation_ClusterModelGeneration_ class to initiate network generation. This class calls objects that inherent the abstract class _Abstract_Runnable_ClusterModel_ContactMap_Generation_ .

Available implementations include:
* Runnable_ClusterModel_ContactMap_Generation_SingleMap(#runnableclusterModelcontactMapgenerationsingleMap)
* Runnable_ClusterModel_ContactMap_Generation_MultiMap(#runnableclusterModelcontactMapgenerationmultiMap)

Please refer to the additional text linked below for more detailed information on each implementation.

### Usage (Contact map generation)
<pre>
java -jar ClusterModel.jar -gen <b><i>Working_Directory</i></b>
</pre>
Arguments:
* <b><i>File_Path_Working_Directory</i></b>: (Required) Path to the working directory where the simulation will run.

### Runnable_ClusterModel_ContactMap_Generation_SingleMap

To be added

### Runnable_ClusterModel_ContactMap_Generation_MultiMap

To be added


## 2. Simulation
Once the network is generated, simulations can be run using model-specific object that inherent the abstract class _Abstract_Runnable_ClusterModel_Transmission_.

In Package_ClusterModel, implementation available includes:
* Runnable_ClusterModel_Transmission, and its extension, Runnable_ClusterModel_Transmission_Map
* Runnable_ClusterModel_MultiTransmission

Please refer to the additional text linked below for more detailed information on each implementation.

> [!NOTE]
> While these implementation can be run on their own (see [Usage](#usage_simulation), it is likely that addtional _Abstract_Runnable_ClusterModel_Transmission_ objects are needed/designed for models with specific purpose. In most cases, an extension of Runnable_ClusterModel_Transmission or Runnable_ClusterModel_MultiTransmission with overriding methods and additional fields are sufficient.   

### Defining simulations setting and intial parameter values 
For both objects, the simulations setting and inital parameter value are specified through an XML configurations file named _simSpecificSim.prop_ in the working directory. Some module also support additional files.

This file follow the Java XML properties format and is read using java.util.Properties#loadFromXML(java.io.InputStream). It has the format as:
```xml
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!DOCTYPE properties SYSTEM "http://java.sun.com/dtd/properties.dtd">
 <properties>
   <entry key="PROP_POP_TYPE">...</entry>
   <entry key="PROP_NUM_SIM_PER_SET">...</entry>
   <entry key="PROP_USE_PARALLEL">...</entry>
   <entry key="PROP_BASESEED">...</entry>
   <entry key="PROP_NUM_SNAP">...</entry>
   <entry key="PROP_SNAP_FREQ">...</entry>
   <entry key="PROP_SIM_SETTING">...</entry>
   <entry key="POP_PROP_INIT_PREFIX_X">...</entry>
   <entry key="POP_PROP_INIT_PREFIX_CLASS_X">...</entry>
   ...
 </properties>
```
Setting that are common to most module includes:
| **Property Name** | **Type** | **Description**  | **Examples** |
| --- | --- | --- | --- |
| `PROP_POP_TYPE`| String | String identifier of simulation type. Some module use this field to specify global parameters. | `<entry key="PROP_POP_TYPE">MultiTransmission_3_4_5</entry>`|  
| `PROP_NUM_SIM_PER_SET` | Integer | Number of simulation to generate. Can be overwritten by custom seed list. | `<entry key="PROP_NUM_SIM_PER_SET">16</entry>` |
| `PROP_USE_PARALLEL` | Integer | Number of simulations run in parallel, or set to less than 1 if parallelisation is not used. | `<entry key="PROP_USE_PARALLEL">16</entry>`| 
|`PROP_BASESEED`| Long | Base seed for the random number generator (RNG), used to derive seeds for each simulation. Can be overwritten by custom seed list| `<entry key="PROP_BASESEED">123456789101112</entry>`|
|`PROP_NUM_SNAP`| Integer | Define the number of snapshots (e.g. prevalence, incidence) to store during simulations. | `<entry key="PROP_NUM_SNAP">30</entry>`|
|`PROP_SNAP_FREQ`|Integer | Define the interval between snapshots (in days). The total simulation duration is therefore `PROP_NUM_SNAP` × `PROP_SNAP_FREQ`| `<entry key="PROP_SNAP_FREQ">365</entry>`|
|`PROP_SIM_SETTING` | Integer | Integer to specify simulation setting. Setting is active if simSetting & 1 << SIM_SETTING_KEY != 0. Note that these setting are Runnable specific and might not be all supported. | `<entry key="PROP_SIM_SETTING">239</entry>`|
|`POP_PROP_INIT_PREFIX_X`|String|Specify the parameter values for the X-th parameter used in model generation. Refer to the respective runnable class for details |
|`POP_PROP_INIT_PREFIX_CLASS_X`|String|Specify the Java classes for the X-th parameter used in model generation. Refer to the respective runnable classes for details| `<entry key="POP_PROP_INIT_PREFIX_CLASS_6">[I</entry>` for integer array (int[]), `<entry key="POP_PROP_INIT_PREFIX_CLASS_16">[[D</entry>` for a double array of two dimension (i.e. double[][]) |


### Runnable_ClusterModel_Transmission and Runnable_ClusterModel_Transmission_Map

To be added

### Runnable_ClusterModel_MultiTransmission

This object is designed to provide the necessary mechanism to support a transmisson model of multiple infections through the same network. 

The key parameters for this object include:


 

#### Usage (Simulation)
Once the network is ready, run simulations using:
```
java -jar ClusterModel.jar -trans Working_Directory
```
* _ClusterModel.jar_: Compiled JAR file.
* _Working_Directory_: Directory containing required files (e.g., XML configurations _simSpecificSim.prop_).

This will then call the launch method in _Simulation_ClusterModelTransmission_, which called upon _generateDefaultRunnable_ method to generate necessary Abstract_Runnable_ClusterModel_Transmission_ object to run a single simulation. 

>[!IMPORTANT]
>The **-trans** command is primarily intended for models investigating bridging dynamics and is retained for backward compatibility. We are actively transitioning toward a modular architecture, where each model is encapsulated in its own module with dedicated parameters and simulation logic. As a result, this generic method may be deprecated in future releases. Users are strongly encouraged to adopt module-specific commands for improved clarity, maintainability, and future support. Please refer to the README within each module for the latest usage instructions.

_Simulation_ClusterModelTransmission_ also designed such that an user-defined _Abstract_Runnable_ClusterModel_Transmission_ can be called upon. This appoached has been used in most Launcher module, where an subclass of _Simulation_ClusterModelTransmission_ are defined and called upon instead. See code snippet below for an example that called a _Runnable_Custom_ object for simulation.

```java
public class Simulation_Custom extends Simulation_ClusterModelTransmission {	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		final String USAGE_INFO = ...;	
		...			
		if (args.length < 1) {
			System.out.println(USAGE_INFO);
			System.exit(0);
		} else {
			Simulation_ClusterModelTransmission.launch(args, new Simulation_Custom());
		}
	}
	
	
	@Override
	public Abstract_Runnable_ClusterModel_Transmission generateDefaultRunnable(long cMap_seed, long sim_seed, 
			Properties loaProperties) {
		...				
		return new Runnable_Custom(...);
	}
	

}


```


