package test;

import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;

import util.Util_MultDirs_Results;

public class Test_PostSim_Analysis_Viability_Baseline {

	public static void main(String[] args) throws IOException {
		File[] base_dirs = new File[] { new File(
				"C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_Baseline_Resample_00"), };

		Pattern pattern_dir_select = Pattern.compile("MSM_Viability_Baseline_Resample_\\d+_Extra_(\\d+)");

		String[] compareFilePrefix = new String[] { "Treatment_Non_Viable_Site", "Treatment_Person", "Incidence_Person",
				"Infectious_Prevalence_Site" };

		int[][] col_sel_all = new int[][] { //
				new int[] { 2, 3, 4, 5, 6, 7 }, // Time,Inf_0_Gender_2_Site_0,Inf_1_Gender_2_Site_1,Inf_1_Gender_2_Site_2,Inf_1_Gender_2_Site_3,Inf_2_Gender_2_Site_1,Inf_2_Gender_2_Site_2,Inf_2_Gender_2_Site_3
				new int[] { 2, 3 }, // Time,Inf_0_Gender_2,Inf_1_Gender_2,Inf_2_Gender_2
				new int[] { 2, 3 }, // Time,Inf_0_Gender_2,Inf_1_Gender_2,Inf_2_Gender_2
				new int[] { 2, 3, 4, 5, 6, 7 }, // Time,Inf_0_Site_0,Inf_1_Site_1,Inf_1_Site_2,Inf_1_Site_3,Inf_2_Site_1,Inf_2_Site_2,Inf_2_Site_3
		};

		File parmeter_list = null;	
		
//		parmeter_list = new File(
//				"C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\Constriant\\MSM_Viability_Baseline_Sel_Seed_2000\\1729491404705_Parameter_Sel.csv");

		for (File result_dir_base : base_dirs) {
			long tic = System.currentTimeMillis();

			String[] preDefinedSeedOrderList = parmeter_list != null? util.Util_MultDirs_Results
					.generatePreDefinedSeedOrderList(result_dir_base, pattern_dir_select, parmeter_list, 25) : null;
			
			Util_MultDirs_Results.printExtractedResultFromMultiDirs(result_dir_base, pattern_dir_select, compareFilePrefix,
					col_sel_all, preDefinedSeedOrderList);
			
			System.out.printf("Results from %s extracted. Time req = %.3fs\n", result_dir_base.getAbsolutePath(),
					(System.currentTimeMillis() - tic) / 1000f);
		}

	}

}
