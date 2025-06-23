# Overview
This Java package provides core utilities for implementing various Individual-Based Models (IBM), particularly those involving cluster-based simulations. It supports both network generation and simulation execution, serving as the backbone for multiple IBM implementations.

# Maintainer and developers
* Ben Hui; ORCiD ID: [0000-0002-6567-5821](https://orcid.org/0000-0002-6567-5821)

# Key Concepts
IBM simulations built with this package typically follow two main steps:

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

The package included the _Simulation_ClusterModelGeneration_ class to initiate network generation. This class calls objects that implement the _Abstract_Runnable_ClusterModel_ContactMap_Generation_ interface.

Available implementations include:
* _Runnable_ClusterModel_ContactMap_Generation_SingleMap_
* _Runnable_ClusterModel_ContactMap_Generation_MultiMap_
  
Refer to the **Javadoc** for each class for detailed usage and configuration.

### Command-Line Usage
```
java -jar ClusterModel.jar -gen Working_Directory
```
* _ClusterModel.jar_: Compiled JAR file.
* _Working_Directory_: Directory containing required files (e.g., XML configurations _simSpecificSim.prop_).

## 2. Simulation Execution
Once the network is generated, simulations can be run using model-specific classes that conform to the _Abstract_Runnable_ClusterModel_ interface.

### Simulation Execution
Once the network is ready, run simulations using:
```
java -jar ClusterModel.jar -trans Working_Directory
```
* _ClusterModel.jar_: Compiled JAR file.
* _Working_Directory_: Directory containing required files (e.g., XML configurations _simSpecificSim.prop_).

>[!IMPORTANT]
>The **-trans** command is primarily intended for models investigating bridging dynamics and is retained for backward compatibility. We are actively transitioning toward a modular architecture, where each model is encapsulated in its own module with dedicated parameters and simulation logic. As a result, this generic method may be deprecated in future releases. Users are strongly encouraged to adopt module-specific commands for improved clarity, maintainability, and future support. Please refer to the README within each module for the latest usage instructions.
