package test;

import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;

import sim.Simulation_ClusterModelTransmission;
import util.Util_Analysis_Compare_Cumulative;
import util.Util_SimulationDirectoryModifications;

public class Test_PostSim_Analysis_Bridging {

	public static void main(String[] args) throws IOException {			
		File[] compareDirs = new File[] { 
				new File("C:\\Users\\bhui\\Documents\\Java_Test\\Results_Bridging\\Baseline_Split_Low_Test_Rate\\Baseline_Split_0"),
				new File("C:\\Users\\bhui\\Documents\\Java_Test\\Results_Bridging\\Baseline_Split_Low_Test_Rate\\Baseline_Split_1"),				
		};
		
		String[] compareFilePrefix = new String[] { "Treatment_Person", "Incidence_Person"};
		File dataDumpCSV_dir = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Results_Bridging\\Baseline_Split_Low_Test_2");		
		String[] time_to_match = new String[] { "6935", "7300"};
		
		for (File f : compareDirs) {
			Util_SimulationDirectoryModifications.generateUniqueZip(f,
					Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("_%d", "") + ".7z",
					Pattern.compile(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("%d",
							"(-{0,1}\\\\d+)")));
		}
		
		Util_Analysis_Compare_Cumulative analysis = new Util_Analysis_Compare_Cumulative(null, compareDirs,
				compareFilePrefix, time_to_match, dataDumpCSV_dir);

		analysis.generateAnalysisCSV();
	}

}
