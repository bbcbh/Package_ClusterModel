package test;

import java.io.File;
import java.io.IOException;

import util.Util_RiskGrpAllocation;

public class Test_Util_RiskGrpAllocation {

	public static void main(String[] args) throws NumberFormatException, IOException, InterruptedException {
	
		// From https://www.unsw.edu.au/research/csrh/our-projects/gay-community-periodic-surveys, Year 2021, 
		// see MSM_Behaviour_Calculation.xlsx
 
		File cMapFolder = new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\Test_cMap");
		float[][] riskCatListAll = new float[][] {
				new float[] { -12, 3, 4, 2, 11, 31, 3, 2, 2, 1, 2, 7, 9, 4, 30, 15, 7, 1 } };
	
		System.out.printf("Starting risk group alloaction from CMap at %s.\n", cMapFolder.getAbsolutePath());
		Util_RiskGrpAllocation.generateRiskGrpAllocationByRiskCat(cMapFolder, riskCatListAll ,1);
		System.out.println("All done.");
	
	}

}
