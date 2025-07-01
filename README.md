# Package_ClusterModel

## Overview
**ClusterModel** is a Java package designed to support the development and execution of Individual-Based Models (IBMs), particularly those involving cluster-based simulations. It provides tools for:
- [**Network generation**](#1-network-generation): Creating contact maps that define partnerships within a population.
- [**Simulation execution**](#2-simulation): Running transmission models over the generated networks.

This package serves as the core engine for multiple IBM implementations.

IBM simulations built with this package typically follow two main steps: [Network generation] and [Simulation](#2_simulation). 

## Maintainer and developers
* Ben Hui; ORCiD ID: [0000-0002-6567-5821](https://orcid.org/0000-0002-6567-5821)


## 1. Network Generation
Networks are defined as CSV files with the following format:
```
person_id_1, person_id_2, partnership_start, partnership_duration
```
* **person_id_1**, **person_id_2**: Unique identifiers for individuals.
* **partnership_start**: Day the partnership begins.
* **partnership_duration**: Duration of the partnership in days.
  
>[!Note]
>Network generation is independent of the simulation process. Any method can be used to generate networks, as long as the format is respected.

The package included the _Simulation_ClusterModelGeneration_ class to initiate network generation. This class calls objects that inherent the abstract class _Abstract_Runnable_ClusterModel_ContactMap_Generation_ .

### Available implementations
Network generation is initiated via the _Simulation_ClusterModelGeneration_ class, which uses implementations of the abstract class _Abstract_Runnable_ClusterModel_ContactMap_Generation_

Implementation include in Package_ClusterModel are:

* Runnable_ClusterModel_ContactMap_Generation_SingleMap(#runnable_clusterModel_contactMap_generation_singleMap)
* Runnable_ClusterModel_ContactMap_Generation_MultiMap(#runnable_clusterModel_contactMap_generation_multiMap)

Refer to the additional text linked below for more detailed information on each object.

#### Runnable_ClusterModel_ContactMap_Generation_SingleMap

To be added

#### Runnable_ClusterModel_ContactMap_Generation_MultiMap

To be added

### Usage (Contact map generation)
<pre>
java -jar ClusterModel.jar -gen <b><i>Working_Directory</i></b>
</pre>
Arguments:
* <b><i>Working_Directory</i></b>: (Required) Path to the directory where the contact map are generated.

## 2. Simulation
Simulations are executed using classes that extend _Abstract_Runnable_ClusterModel_Transmission_.


### Available implementations
* [Runnable_ClusterModel_Transmission, and its extension, Runnable_ClusterModel_Transmission_Map](#runnable_clustermodel_transmission-and-runnable_clustermodel_transmission_map)
* [Runnable_ClusterModel_MultiTransmission](#runnable_clustermodel_multitransmission)

> [!NOTE]
> These implementations can be extended for model-specific needs by overriding methods and adding fields.

### Configuration
Simulations are configured using an XML file named _simSpecificSim.prop_ in the working directory. This file follows the Java XML properties format:

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
#### Common properties 
| **Property Name** | **Type** | **Description**  | **Examples** |
| --- | --- | --- | --- |
| `PROP_POP_TYPE`| String | Simulation type identifier | `<entry key="PROP_POP_TYPE">MultiTransmission_3_4_5</entry>`|  
| `PROP_NUM_SIM_PER_SET` | Integer | Number of simulations. Can be overwritten by custom seed list. | `<entry key="PROP_NUM_SIM_PER_SET">16</entry>` |
| `PROP_USE_PARALLEL` | Integer | Number of parallel simulations, or set to less than 1 if parallelisation is not used. | `<entry key="PROP_USE_PARALLEL">16</entry>`| 
|`PROP_BASESEED`| Long | Base random number generator (RNG) seed used to derive seeds for each simulation. Can be overwritten by custom seed list| `<entry key="PROP_BASESEED">123456789101112</entry>`|
|`PROP_NUM_SNAP`| Integer | Number of snapshots (e.g. prevalence, incidence) to store during simulations. | `<entry key="PROP_NUM_SNAP">30</entry>`|
|`PROP_SNAP_FREQ`|Integer | Days between snapshots. The total simulation duration is therefore `PROP_NUM_SNAP` × `PROP_SNAP_FREQ`| `<entry key="PROP_SNAP_FREQ">365</entry>`|
|`PROP_SIM_SETTING` | Integer | Bitmask for simulation settings. Note that these setting are Runnable specific and might not be all supported. | `<entry key="PROP_SIM_SETTING">239</entry>`|
|`POP_PROP_INIT_PREFIX_X`|String|Initial parameter values for the X-th parameter used in model simulations. Refer to the respective runnable class for details |`<entry key="POP_PROP_INIT_PREFIX_6">[0, 0, 104000, 0]</entry>`| 
|`POP_PROP_INIT_PREFIX_CLASS_X`|String|Specify the Java classes for the X-th parameter used in model simulations. Refer to the respective runnable classes for details| `<entry key="POP_PROP_INIT_PREFIX_CLASS_6">[I</entry>` for integer array (int[])<br/> `<entry key="POP_PROP_INIT_PREFIX_CLASS_16">[[D</entry>` for a double array of two dimension (i.e. double[][]) |


### Runnable_ClusterModel_Transmission and Runnable_ClusterModel_Transmission_Map

To be added

### Runnable_ClusterModel_MultiTransmission

This object is designed to provide the necessary mechanism to support a transmisson model of multiple infections through the same network. 

The key parameters for this object include:
| **Property Name** | **Type** | **Description**  | **Examples** |
| --- | --- | --- | --- |

 #### Usage (Simulation)
<pre>
java -jar ClusterModel.jar -trans <b><i>Working_Directory</i></b>
</pre>
* <b><i>Working_Directory</i></b>: (Required) Directory containing _simSpecificSim.prop_ and other required files.

This will then call the launch method in _Simulation_ClusterModelTransmission_, which called upon _generateDefaultRunnable_ method to generate necessary Abstract_Runnable_ClusterModel_Transmission_ object to run a single simulation. 

>[!IMPORTANT]
>The **-trans** command is primarily intended for models investigating bridging dynamics and is retained for backward compatibility. We are actively transitioning toward a modular architecture, where each model is encapsulated in its own module with dedicated parameters and simulation logic. As a result, this generic method may be deprecated in future releases. Users are strongly encouraged to adopt module-specific commands for improved clarity, maintainability, and future support. Please refer to the README within each module for the latest usage instructions.

#### Custom simulation example

You can define a custom simulation by extending _Simulation_ClusterModelTransmission_:

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
		return new Runnable_Custom(...); // for example
	}
}
```


