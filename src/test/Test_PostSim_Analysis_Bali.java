package test;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import sim.Simulation_ClusterModelTransmission;
import util.Util_7Z_CSV_Entry_Extract_Callable;
import util.Util_MultDirs_Results;

public class Test_PostSim_Analysis_Bali {

	public static void main(String[] arg) throws FileNotFoundException, IOException {

		boolean analyse_combined = false; // Pre 20241115 version

		File baseDir_sims = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Bali");
		Pattern pattern_simDir = Pattern.compile("Bali_.*");
		File baseSim_dir = new File(baseDir_sims, "Bali_Baseline");
		String res_prefix = "Incidence_Person_";
		int numInf = 4;
		double py_conv_factor = 100.0 *(1.0/14000) *(1.0/5);  

		File[] sim_dirs = baseDir_sims.listFiles(new FileFilter() {

			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && !pathname.equals(baseSim_dir)
						&& pattern_simDir.matcher(pathname.getName()).matches();
			}
		});

		System.out.printf("Comparing %d results with those in %s.\n", sim_dirs.length, baseDir_sims.getAbsolutePath());

		HashMap<String, HashMap<String, ArrayList<String[]>>> resMapBase = Util_MultDirs_Results
				.extractedResultsFromSimDirBase(baseSim_dir, res_prefix);
		
		HashMap<String, double[]> baseIncidenceMap = extractedIncidenceMap(numInf, resMapBase);		
	
		
		Percentile ptile = new Percentile();
		double[][] raw_value = new double[numInf][baseIncidenceMap.size()];
		int counter = 0;
		for(Entry<String, double[]> inc : baseIncidenceMap.entrySet()) {
			for(int i = 0; i < numInf; i++) {
				raw_value[i][counter] = inc.getValue()[i+1];
			}						
			counter++;						
		}	
		
		System.out.println("Baseline incidence over 5 years:");
		for(int i = 0; i < numInf; i++) {
			ptile.setData(raw_value[i]);
			System.out.printf("Inf #%d = %.1f (%.1f - %.1f)\n", i,
					ptile.evaluate(50)*py_conv_factor, ptile.evaluate(25)*py_conv_factor, ptile.evaluate(75)*py_conv_factor);
			
		}
		
		for(File sim_dir : sim_dirs) {								
			HashMap<String, HashMap<String, ArrayList<String[]>>> resMap = Util_MultDirs_Results
					.extractedResultsFromSimDirBase(sim_dir, res_prefix);
			if(!resMap.isEmpty()) {			
				System.out.printf("Comparing incidence from %s...\n", sim_dir.getName());
				
				HashMap<String, double[]> incidenceMap = extractedIncidenceMap(numInf, resMap);
				
				double[][] raw_value_sim = new double[numInf][incidenceMap.size()];
				int counter_sim = 0;
				for(Entry<String, double[]> inc : incidenceMap.entrySet()) {
					double[] baseEnt = baseIncidenceMap.get(inc.getKey());					
					for(int i = 0; i < numInf; i++) {
						raw_value_sim[i][counter_sim] = 1 - inc.getValue()[i+1]/ baseEnt[i+1];
					}						
					counter_sim++;						
				}				
				
				System.out.println("Reduction of incidence 5 years:");
				for(int i = 0; i < numInf; i++) {
					ptile.setData(raw_value_sim[i]);
					System.out.printf("Inf #%d = %.1f%% (IQR: %.1f - %.1f%%)\n", i,
							ptile.evaluate(50)*100, ptile.evaluate(25)*100, ptile.evaluate(75)*100);
					
				}
				
				
				
			}
		
		}
		
		
		
		
		

		if (analyse_combined) {

			File baseDir = new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\Bali\\Bali_Baseline");
			// output_analysis(baseDir);
			// System.out.printf("Analysis of %s completed.\n", baseDir.getAbsolutePath());

			/// Compare Results
			File[] compareDirs = new File[] {
					new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\Bali\\Bali_PrEP"),
					new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\Bali\\Bali_PrEP_P30"), };
			File compareBase = new File(baseDir.getParentFile(), "Compare");
			compareBase.mkdirs();

			String[] zipFormat_arr = new String[] { "Incidence_Person_" };
			for (String zipFormat : zipFormat_arr) {
				String zipFileName = String.format("%sAll.csv.7z", zipFormat);
				HashMap<String, ArrayList<String[]>> baseline_ent = new HashMap<>();
				baseline_ent = Util_7Z_CSV_Entry_Extract_Callable.extractedLinesFrom7Zip(new File(baseDir, zipFileName),
						baseline_ent);

				String[] fNameKey = baseline_ent.keySet().toArray(new String[0]);
				Arrays.sort(fNameKey);

				StringBuilder[][] baseline_lines = colSpecificEntriesFromExtractedLines(baseline_ent, fNameKey);

				HashMap<File, StringBuilder[][]> compare_lines_map = new HashMap<>();
				for (File compareDir : compareDirs) {
					HashMap<String, ArrayList<String[]>> compare_ent = Util_7Z_CSV_Entry_Extract_Callable
							.extractedLinesFrom7Zip(new File(compareDir, zipFileName));
					compare_lines_map.put(compareDir, colSpecificEntriesFromExtractedLines(compare_ent, fNameKey));
				}

				for (int c = 1; c < baseline_lines.length; c++) {
					File compareFile = new File(compareBase, String.format("Diff_%sAll_C%d.csv", zipFormat, c));
					PrintWriter pWri = new PrintWriter(compareFile);

					pWri.print("Time");
					for (String fName : fNameKey) {
						pWri.print(',');
						pWri.print(fName.replaceAll(",", "_"));
					}
					pWri.println();
					pWri.println(baseDir.getAbsolutePath());
					for (int r = 1; r < baseline_lines[c].length; r++) {
						pWri.println(baseline_lines[c][r].toString());
					}
					for (File compareDir : compareDirs) {
						pWri.println();
						pWri.println(compareDir.getAbsolutePath());
						StringBuilder[][] compare_lines = compare_lines_map.get(compareDir);

						for (int r = 1; r < compare_lines[c].length; r++) {
							pWri.println(compare_lines[c][r].toString());
						}

					}

					pWri.close();
				}
			}
		}

	}

	public static HashMap<String, double[]> extractedIncidenceMap(int numInf,
			HashMap<String, HashMap<String, ArrayList<String[]>>> resMapBase) {
		HashMap<String, double[]> incidenceMap = new HashMap<>();
		for(Entry<String, HashMap<String, ArrayList<String[]>>> ent : resMapBase.entrySet()) {
			for(Entry<String, ArrayList<String[]>> subEnt : ent.getValue().entrySet()){
				String key = String.format("%s:%s", ent.getKey(), subEnt.getKey());
				double[] val = new double[numInf+1]; // 4 infections + time
				val[0] = Double.NaN;
				for(String[] line : subEnt.getValue()) {					
					if(line[0].equals("7300")) {
						for(int i = 1; i < line.length; i++) {
							val[i] -= Double.parseDouble(line[i]);
						}												
					}else if(line[0].equals("9125")) {
						for(int i = 1; i < line.length; i++) {
							val[i] += Double.parseDouble(line[i]);
						}						
					}										
				}				
				incidenceMap.put(key, val);				
			}			
		}
		return incidenceMap;
	}

	protected static StringBuilder[][] colSpecificEntriesFromExtractedLines(
			HashMap<String, ArrayList<String[]>> fName_csvLine_map, String[] fNameKey) {
		int numRow = fName_csvLine_map.get(fNameKey[0]).size();
		int numCol = fName_csvLine_map.get(fNameKey[0]).get(0).length;

		StringBuilder[][] lines = new StringBuilder[numCol][numRow];

		for (String fNK : fNameKey) {
			ArrayList<String[]> text_map = fName_csvLine_map.get(fNK);
			if (text_map == null) {
				System.err.printf("Warning: %s not found in map!\n", fNK);
			} else {

				for (int r = 1; r < numRow; r++) { // Skip first line
					String[] row = text_map.get(r);
					for (int c = 1; c < numCol; c++) {
						if (lines[c][r] == null) {
							lines[c][r] = new StringBuilder();
							lines[c][r].append(row[0]);
						}
						lines[c][r].append(',');
						lines[c][r].append(row[c]);
					}
				}
			}
		}
		return lines;
	}

	public static void output_analysis(File baseDir) throws IOException, FileNotFoundException {
		// File_Format
		String[] zipFormat = new String[] { "Incidence_Person_", "Incidence_Site_",
				"Infected_All_Stages_Prevalence_Site_", "Infected_Prevalence_Site_", "Infectious_Prevalence_Person_",
				"Infectious_Prevalence_Site_", "Treatment_Person_" };

		String[] stat_file_name = new String[] { "Summary_Incidence_Person_%s.csv", null, null, null,
				"Summary_Infectious_Prevalence_Person_%s.csv", "Summary_Infectious_Prevalence_Site_%s.csv", null };

		Simulation_ClusterModelTransmission.output_analysis_csv(baseDir, zipFormat, stat_file_name);
	}

}
