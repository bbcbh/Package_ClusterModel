package test;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import random.MersenneTwisterRandomGenerator;
import util.Util_MultDirs_Results;
import util.Util_Select_Sim_By_Residue;

public class Test_PostSim_Analysis_Viability_Resample_By_Residue {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		File baseDir = new File(
				"C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_Baseline_Resample");

		Pattern directory_pattern =
				//Pattern.compile("MSM_Viability_Baseline_Resample_\\d\\d\\d_\\d\\d\\d"); // Version 1 only
				//Pattern.compile("MSM_Viability_Baseline_Resample_Directed_.*"); // Direct only
				Pattern.compile("MSM_Viability_Baseline_Resample_.*"); // All

		File paramListSummaryFile = new File(baseDir, "Extract_ParamListSummary");
		paramListSummaryFile.mkdirs();
		paramListSummaryFile = new File(paramListSummaryFile, "Parameter_List.csv");

		boolean reuseParameterList = true;
		boolean genNextSim = !true;
	

		File importDirsBase = new File(
				"C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_Baseline_Resample\\Import_Results");
		// Move result folder
		File[] importDirs = importDirsBase.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory();
			}
		});

		for (File importDir : importDirs) {
			int suc = Util_MultDirs_Results.combineResultDirectories(baseDir, importDir);
			System.out.printf("%d results directories moved from %s\n   to %s.\n", suc, importDir.getAbsolutePath(),
					baseDir.getAbsolutePath());
		}

		// Generate parameter list
		if (reuseParameterList && paramListSummaryFile.exists()) {
			System.out.printf("Reused resample seed based on parameter list from  %s.\n",
					paramListSummaryFile.getAbsolutePath());
		} else {

			System.out.printf("Generating new parameter list...\n",
					paramListSummaryFile.getAbsolutePath());

			double[][] residue_target = new double[][] { //
					// from AW email at 20240905
					null, // Urethral NG culture - almost 100% (not fitted)
					new double[] { 100 - 63.03 }, // Rectal NG culture: Total N = 238, 150 positive (63.03%)
					new double[] { 100 - 31.63 }, // Throat NG culture: Total N = 332, 105 positive (31.63%)
					new double[] { 19, 13, 26 }, // From review paper: 19% non-viable for urtheral
					new double[] { 33, 23, 42 }, // 33% non-viable for rectal
					null, new double[] { 23.9 * 0.9 + 31.5 * 0.1, 23.9, 31.5 }, // NG: Surv Report 2023
					new double[] { 29.2 * 0.9 + 40.9 * 0.1, 29.2, 40.9 }, // CT: Surv Report 2023
//					// Site specific preval from Chow 2019
//					null, // 
//					new double[] { 5.9, 0.2, 24 }, 
//					new double[] { 4.6, 0.5, 16.5 }, 
//					null,
//					new double[] { 8.9, 2.1, 23.0 }, 
//					new double[] { 1.7, 0, 3.6 },
					// Site specific preval from Chow 2015 (Aus. Only)
					new double[] { 2.3 },
					new double[] { 2.9 }, 
					new double[] { 1.7 }, 
					new double[] { 3.0 }, 
					new double[] { 5.6 }, 
					new double[] { 1.7 }											
			};

			double nv_weight = 1;
			double incidence_weight = 1;
			double preval_weight = 10;

			double[] residue_weight = new double[] { 0, nv_weight, nv_weight, nv_weight, nv_weight, nv_weight, 
					incidence_weight, incidence_weight,	preval_weight, preval_weight, preval_weight,
					preval_weight, preval_weight, preval_weight };

			Util_Select_Sim_By_Residue analysis = new Util_Select_Sim_By_Residue(baseDir, directory_pattern) {

				int RM_INDEX_NonViable_NG_U = 0; 
				int RM_INDEX_NonViable_NG_R = RM_INDEX_NonViable_NG_U + 1; 
				int RM_INDEX_NonViable_NG_P = RM_INDEX_NonViable_NG_R + 1;
				int RM_INDEX_NonViable_CT_U = RM_INDEX_NonViable_NG_P + 1;
				int RM_INDEX_NonViable_CT_R = RM_INDEX_NonViable_CT_U + 1;
				int RM_INDEX_NonViable_CT_P = RM_INDEX_NonViable_CT_R + 1;
				int RM_INDEX_Incidence_NG = RM_INDEX_NonViable_CT_P + 1;
				int RM_INDEX_Incidence_CT = RM_INDEX_Incidence_NG + 1;
				int RM_INDEX_Preval_NG_U = RM_INDEX_Incidence_CT + 1;
				int RM_INDEX_Preval_NG_R = RM_INDEX_Preval_NG_U + 1;
				int RM_INDEX_Preval_NG_P = RM_INDEX_Preval_NG_R + 1;
				int RM_INDEX_Preval_CT_U = RM_INDEX_Preval_NG_P + 1;
				int RM_INDEX_Preval_CT_R = RM_INDEX_Preval_CT_U + 1;
				int RM_INDEX_Preval_CT_P = RM_INDEX_Preval_CT_R + 1;
				int LENGTH_RM = RM_INDEX_Preval_CT_P + 1;

				@Override
				public HashMap<String, double[]> generateResidueMapping() {

					HashMap<String, double[]> mapping = new HashMap<>();

					File[] resultDir = baseDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return pathname.isDirectory() && directory_pattern.matcher(pathname.getName()).matches();
						}
					});

					double[] res_entry;

					for (File resDir : resultDir) {
						try {
							int select_row = 19; // T = 6935

							// Incidence
							File[] file_incdence = resDir.listFiles(new FileFilter() {
								@Override
								public boolean accept(File pathname) {
									return Pattern.matches(".*Incidence_Person_.*", pathname.getName());
								}
							});
							HashMap<String, ArrayList<String[]>> entMap_incidence = Util_Select_Sim_By_Residue
									.combineResultMapping(file_incdence);

							for (Entry<String, ArrayList<String[]>> ent : entMap_incidence.entrySet()) {
								String map_key = String.format("%s:%s", resDir.getName(), ent.getKey());
								res_entry = getResidueMappingEntry(mapping, map_key);
								ArrayList<String[]> lines = ent.getValue();
								String[] line_current_line = lines.get(select_row);
								String[] iine_past_line = lines.get(select_row - 1);
								res_entry[RM_INDEX_Incidence_NG] = (Double.parseDouble(line_current_line[2])
										- Double.parseDouble(iine_past_line[2])) / 1040;
								res_entry[RM_INDEX_Incidence_CT] = (Double.parseDouble(line_current_line[3])
										- Double.parseDouble(iine_past_line[3])) / 1040;
							}

							// Site_Specific prevalence
							File[] file_prevalence_site = resDir.listFiles(new FileFilter() {
								@Override
								public boolean accept(File pathname) {
									return Pattern.matches(".*Infectious_Prevalence_Site_.*", pathname.getName());
								}
							});
							HashMap<String, ArrayList<String[]>> entMap_preval = Util_Select_Sim_By_Residue
									.combineResultMapping(file_prevalence_site);

							for (Entry<String, ArrayList<String[]>> ent : entMap_preval.entrySet()) {
								String map_key = String.format("%s:%s", resDir.getName(), ent.getKey());
								res_entry = getResidueMappingEntry(mapping, map_key);
								ArrayList<String[]> lines = ent.getValue();
								String[] line_current_line = lines.get(select_row);
								res_entry[RM_INDEX_Preval_NG_U] = Double.parseDouble(line_current_line[2]) / 1040;
								res_entry[RM_INDEX_Preval_NG_R] = Double.parseDouble(line_current_line[3]) / 1040;
								res_entry[RM_INDEX_Preval_NG_P] = Double.parseDouble(line_current_line[4]) / 1040;
								res_entry[RM_INDEX_Preval_CT_U] = Double.parseDouble(line_current_line[5]) / 1040;
								res_entry[RM_INDEX_Preval_CT_R] = Double.parseDouble(line_current_line[6]) / 1040;
								res_entry[RM_INDEX_Preval_CT_P] = Double.parseDouble(line_current_line[7]) / 1040;
							}

							// Non viability
							File[] file_nv = resDir.listFiles(new FileFilter() {
								@Override
								public boolean accept(File pathname) {
									return Pattern.matches(".*Treatment_Non_Viable_Site_.*", pathname.getName());
								}
							});
							HashMap<String, ArrayList<String[]>> entMap_nv = Util_Select_Sim_By_Residue
									.combineResultMapping(file_nv);

							// Treatment
							File[] file_treatment = resDir.listFiles(new FileFilter() {
								@Override
								public boolean accept(File pathname) {
									return Pattern.matches(".*Treatment_Person_.*", pathname.getName());
								}
							});
							HashMap<String, ArrayList<String[]>> entMap_treatment = Util_Select_Sim_By_Residue
									.combineResultMapping(file_treatment);
							for (Entry<String, ArrayList<String[]>> ent : entMap_nv.entrySet()) {
								String map_key = String.format("%s:%s", resDir.getName(), ent.getKey());
								res_entry = getResidueMappingEntry(mapping, map_key);
								String[] nv_line = ent.getValue().get(select_row);
								String[] treatment_line = entMap_treatment.get(ent.getKey()).get(select_row);
								res_entry[RM_INDEX_NonViable_NG_U] = 100 * Double.parseDouble(nv_line[2])
										/ Double.parseDouble(treatment_line[2]);
								res_entry[RM_INDEX_NonViable_NG_R] = 100 * Double.parseDouble(nv_line[3])
										/ Double.parseDouble(treatment_line[2]);
								res_entry[RM_INDEX_NonViable_NG_P] = 100 * Double.parseDouble(nv_line[4])
										/ Double.parseDouble(treatment_line[2]);
								res_entry[RM_INDEX_NonViable_CT_U] = 100 * Double.parseDouble(nv_line[5])
										/ Double.parseDouble(treatment_line[3]);
								res_entry[RM_INDEX_NonViable_CT_R] = 100 * Double.parseDouble(nv_line[6])
										/ Double.parseDouble(treatment_line[3]);
								res_entry[RM_INDEX_NonViable_CT_P] = 100 * Double.parseDouble(nv_line[7])
										/ Double.parseDouble(treatment_line[3]);
							}

						} catch (IOException ex) {
							ex.printStackTrace(System.err);
						}

					}
					return mapping;
				}

				private double[] getResidueMappingEntry(HashMap<String, double[]> mapping, String map_key) {
					double[] res_entry;
					res_entry = mapping.get(map_key);
					if (res_entry == null) {
						res_entry = new double[LENGTH_RM];
						Arrays.fill(res_entry, Double.NaN);
						mapping.put(map_key, res_entry);
					}
					return res_entry;
				}

			};

			HashMap<String, double[]> mapping = analysis.generateResidueMapping();

			ArrayList<Double> ordered_residue_val = new ArrayList<>();
			ArrayList<String> ordered_keys = new ArrayList<>();

			ArrayList<String> inRangeKeyArray = Util_Select_Sim_By_Residue.select_best_sim(mapping, residue_target,
					residue_weight, ordered_residue_val, ordered_keys);

			System.out.printf("%d results in total, with %d simulation results fitted in selection criteria.\n",
					mapping.size(), inRangeKeyArray.size());

			if (paramListSummaryFile.exists()) {
				File tar_summary = new File(paramListSummaryFile.getParent(),
						String.format("%d_%s", System.currentTimeMillis(), paramListSummaryFile.getName()));
				try {
					Files.move(paramListSummaryFile.toPath(), tar_summary.toPath(), StandardCopyOption.ATOMIC_MOVE);
				} catch (Exception e) {
					e.printStackTrace(System.err);
				}
			}

			Util_Select_Sim_By_Residue.printOrderedParamList(baseDir, ordered_keys, ordered_residue_val,
					inRangeKeyArray, mapping, paramListSummaryFile);

			System.out.printf("Parameter list generated at\n %s.\n", paramListSummaryFile.getAbsolutePath());

		}

		// Resample generate new results
		int genSeedDir_max_result_to_include = 2000; // set to <0 if resample sim are not generated.
		int genSeedDir_MaxSimPerDir = 16;
		int genSeedDir_num_col_include = 25; // Including map and seed id.
		int genSeedDir_in_range_col = 28;		  
		int genSeedDir_num_reseed = 29; 
		int genSeedDir_folder_offset = 20;      
		int genSeedDir_folder_id =  genSeedDir_num_reseed-genSeedDir_folder_offset;


		float genSeedDir_numSimMin = 1000f;
		long genSeedDir_rng_seed = 22519122707291119l;

		File resampleDir = new File(baseDir, "Extract_Resampled_Seed_Sim");

		String resample_dir_name = String.format("MSM_Viability_Baseline_Resample_Directed_%03d", genSeedDir_folder_id); // Adjust
																												// as
																												// needed
		String pattern_PBS_id_line = "#PBS -N MSM_VTD" + genSeedDir_folder_id + "_%03d"; // Adjust as needed

		String pattern_PBS_filename = "ClusterModel_Sim_%03d.pbs";

		File genSeedDir_SIM_DIR = new File(resampleDir, resample_dir_name);
		File genSeedDir_PBS_DIR = new File(resampleDir, "PBS_" + resample_dir_name);
		File genSeedDir_SRC_PBS = new File(
				"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Sim_Template.pbs");
		File genSeedDir_SRC_Template = new File(
				"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Sim_MSM_Viablity");
		String[] genSeedDir_direct_copy = new String[] { "simSpecificSim.prop" };

		if (genNextSim) {

			// Code start
			resampleDir.mkdirs();
			MersenneTwisterRandomGenerator RNG = new MersenneTwisterRandomGenerator(genSeedDir_rng_seed);
			for (int i = 0; i < genSeedDir_num_reseed; i++) {
				RNG = new MersenneTwisterRandomGenerator(RNG.nextLong());
			}

			String[] lines_summary_raw = util.Util_7Z_CSV_Entry_Extract_Callable
					.extracted_lines_from_text(paramListSummaryFile);

			// Select line
			String[] lines_selected = new String[genSeedDir_max_result_to_include + 1];
			System.arraycopy(lines_summary_raw, 0, lines_selected, 0, genSeedDir_max_result_to_include + 1);
			int lastPt = lines_selected.length - 1;

			// Replace less fit line with those that fitted within range
			for (int i = genSeedDir_max_result_to_include + 1; i < lines_summary_raw.length && lastPt > 0; i++) {
				String[] sp = lines_summary_raw[i].split(",");
				if (Boolean.parseBoolean(sp[genSeedDir_in_range_col])) {
					lines_selected[lastPt] = lines_summary_raw[i];
					lastPt--;
				}
			}

			HashMap<Long, HashMap<Integer, ArrayList<Double>>> parameter_map = new HashMap<>(); // Key = mapId, V =
																								// Map<SeedNum,
																								// Parameter>
			ArrayList<String> inRangeLines = new ArrayList<>();

			int includedCounter = 0;

			for (int i = 1; i < lines_selected.length; i++) {
				String[] line_arr = lines_selected[i].split(",");
				Long mapId = Long.parseLong(line_arr[0]);

				HashMap<Integer, ArrayList<Double>> entryMap = parameter_map.get(mapId);
				if (entryMap == null) {
					entryMap = new HashMap<>();
					parameter_map.put(mapId, entryMap);
				}

				for (int c = 2; c < genSeedDir_num_col_include; c++) { // Skip the first two column as they are seed
					ArrayList<Double> entryArr = entryMap.get(c);
					if (entryArr == null) {
						entryArr = new ArrayList<>();
						entryMap.put(c, entryArr);
					}
					entryArr.add(Double.parseDouble(line_arr[c]));
				}
				includedCounter++;

				if (Boolean.parseBoolean(line_arr[genSeedDir_in_range_col])) {
					inRangeLines.add(lines_selected[i]);
				}
			}

			System.out.printf("%d results included for resampling, with %d already in range.\n", includedCounter,
					inRangeLines.size());

			int counter = 0;
			for (Entry<Long, HashMap<Integer, ArrayList<Double>>> entryMap : parameter_map.entrySet()) {
				System.out.printf("Map #%d (seed = %d): Num of resample point = %d.\n", counter, entryMap.getKey(),
						entryMap.getValue().get(2).size());
				counter++;
			}

			Long[] cMap_list = parameter_map.keySet().toArray(new Long[0]);
			Arrays.sort(cMap_list);

			int numDirPerMap = (int) Math.ceil(((genSeedDir_numSimMin / cMap_list.length) / genSeedDir_MaxSimPerDir));

			genSeedDir_SIM_DIR.mkdirs();
			genSeedDir_PBS_DIR.mkdirs();

			// Generate parameter list header
			String[] header_ent = lines_selected[0].split(",");
			StringBuilder header = new StringBuilder(header_ent[0]);
			for (int i = 1; i < genSeedDir_num_col_include; i++) {
				header.append(',');
				header.append(header_ent[i]);
			}

			// Generate seed list for those already in range
			PrintWriter inRangeWriter = new PrintWriter(new File(genSeedDir_SIM_DIR, "Seed_InRange.csv"));
			inRangeWriter.println(header.toString());
			for (String line : inRangeLines) {
				inRangeWriter.println(line);
			}
			inRangeWriter.close();

			// PBS
			String batch_script_filename = "batch_qsub";
			PrintWriter batchScriptWriter = new PrintWriter(new File(genSeedDir_PBS_DIR, batch_script_filename)) {
				@Override
				public void println() {
					write('\n');
				}
			};
			batchScriptWriter.println("#!/bin/bash");
			batchScriptWriter.println("echo \"Batch submiting...\"");

			String[] line_pbs_arr = util.Util_7Z_CSV_Entry_Extract_Callable
					.extracted_lines_from_text(genSeedDir_SRC_PBS);

			int dir_counter = 0;
			LinearInterpolator polator = new LinearInterpolator();

			for (Long cMap_id : cMap_list) {
				HashMap<Integer, ArrayList<Double>> entryMap = parameter_map.get(cMap_id);
				PolynomialSplineFunction[] dist = new PolynomialSplineFunction[genSeedDir_num_col_include];

				for (int p = 2; p < dist.length; p++) { // Offset by cMap_Seed and simSeed
					ArrayList<Double> entArr = entryMap.get(p);
					double[] dist_val = new double[entArr.size()];
					int pt = 0;
					for (Double ent : entArr) {
						dist_val[pt] = ent;
						pt++;
					}
					Arrays.sort(dist_val);

					double[] x_val = new double[dist_val.length];
					for (int x = 1; x < dist_val.length; x++) {
						x_val[x] = x_val[x - 1] + 1.0 / dist_val.length;
					}
					x_val[x_val.length - 1] = 1;
					dist[p] = polator.interpolate(x_val, dist_val);
				}
				for (int d = 0; d < numDirPerMap; d++) {
					File simDir = new File(genSeedDir_SIM_DIR,
							String.format("%s_Extra_%03d", resample_dir_name, dir_counter));
					simDir.mkdirs();

					// Direct copy
					for (String copyFileName : genSeedDir_direct_copy) {
						File copyFile = new File(genSeedDir_SRC_Template, copyFileName);
						if (copyFile.exists()) {
							File tarFile = new File(simDir, copyFile.getName());
							Files.copy(copyFile.toPath(), tarFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
						}
					}

					// Seed File
					File seedFile = new File(simDir, "Seed_List.csv");
					PrintWriter pWri = new PrintWriter(seedFile);
					pWri.println(header.toString());

					for (int r = 0; r < genSeedDir_MaxSimPerDir; r++) {
						pWri.print(cMap_id);
						pWri.print(',');
						pWri.print(RNG.nextLong());
						for (int p = 2; p < dist.length; p++) {
							pWri.print(',');
							double val;
							val = dist[p].value(RNG.nextDouble());
							pWri.print(val);
						}
						pWri.println();
					}

					pWri.close();

					// PBS
					String pBS_FileName = String.format(pattern_PBS_filename, dir_counter);
					PrintWriter pbs_wri = new PrintWriter(new File(genSeedDir_PBS_DIR, pBS_FileName));
					for (int n = 0; n < line_pbs_arr.length; n++) {
						String lineEnt = line_pbs_arr[n];
						switch (n) {
						case 5:
							lineEnt = String.format(pattern_PBS_id_line, dir_counter);
							break;
						case 9:
							lineEnt = String.format("#PBS -l ncpus=%d", genSeedDir_MaxSimPerDir);
							break;
						case 14:
							lineEnt = String.format("#PBS -o ./%s.out", pBS_FileName);
							break;
						case 15:
							lineEnt = String.format("#PBS -e ./%s.err", pBS_FileName);
							break;
						case 23:
							lineEnt = String.format(
									"java -jar \"ClusterModel.jar\" "
											+ "-trans \"/srv/scratch/z2251912/All_Transmission/%s/%s\" "
											+ "-export_skip_backup -printProgress -seedMap=Seed_List.csv",
									genSeedDir_SIM_DIR.getName(), simDir.getName());
						default:

						}

						pbs_wri.println(lineEnt);
					}

					pbs_wri.close();

					batchScriptWriter.printf("qsub %s\n", pBS_FileName);

					dir_counter++;

				}

			}
			batchScriptWriter.close();

			System.out.printf("Resample dir generated at %s.\n", resampleDir.getAbsolutePath());
		}

	}

}
