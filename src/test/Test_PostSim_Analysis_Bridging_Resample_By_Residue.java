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

public class Test_PostSim_Analysis_Bridging_Resample_By_Residue {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		File baseDir = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Bridging\\Baseline_Split_All_Param");

		Pattern directory_pattern = Pattern.compile("SimClusterModel_Transmission_.*");

		File paramListSummaryFile = new File(baseDir, "Extract_ParamListSummary");
		paramListSummaryFile.mkdirs();
		paramListSummaryFile = new File(paramListSummaryFile, "Parameter_List.csv");

		boolean reuseParameterList = !true;
		boolean genNextSim = !true;
		

		File importDirsBase = new File(baseDir, "Import_Results");
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
			System.out.printf("Reused resample seed based on parameter list from  %s. \n",
					paramListSummaryFile.getAbsolutePath());
		} else {

			System.out.printf("Generating new parameter list and store as %s. \n",
					paramListSummaryFile.getAbsolutePath());

			double[][] residue_target = new double[][] { //
					// Site specific preval from Chow 2015 (Aus only)
					new double[] { 2.3 }, // Ure MSM
					new double[] { 2.9 }, // Rect MSM
					new double[] { 1.7 }, // Oro MSM
					null,				  // Ure Female
					null,				  // Ure Male
					new double[] { 41.6 }, // Notification Female 2015
					new double[] { 116.8 }, // Notification All Male 2015 
					new double[] { 24910, 23800, 34900 }, // Incidence MSM 2015 (Assume 10% HIV)
					new double[] { 79.8 }, // Notification Female 2022
					new double[] { 187.3 }, // Notification All Male 2022
					new double[] { 24660, 23900, 31500 }, // Incidence MSM 2022 (Assume 10% HIV)					
			};							
			
			double weight_preval = 1; //1;
			double weight_treatement = 0.1;
			double weight_incidence = 1 / 10000.0;
			
			double[] residue_weight = new double[] { weight_preval, weight_preval, weight_preval, weight_preval, weight_preval,
					weight_treatement, weight_treatement, weight_incidence, weight_treatement, weight_treatement, weight_incidence};

			Util_Select_Sim_By_Residue analysis = new Util_Select_Sim_By_Residue(baseDir, directory_pattern) {
				int RM_INDEX_Preval_MSM_NG_U = 0;
				int RM_INDEX_Preval_MSM_NG_R = RM_INDEX_Preval_MSM_NG_U + 1;
				int RM_INDEX_Preval_MSM_NG_P = RM_INDEX_Preval_MSM_NG_R + 1;
				int RM_INDEX_Preval_FEMALE_NG_U = RM_INDEX_Preval_MSM_NG_P + 1;
				int RM_INDEX_Preval_MALE_HETERO_NG_U = RM_INDEX_Preval_FEMALE_NG_U + 1;
				int RM_INDEX_Treat_FEMALE = RM_INDEX_Preval_MALE_HETERO_NG_U + 1;
				int RM_INDEX_Treat_MALE_ALL = RM_INDEX_Treat_FEMALE + 1;
				int RM_INDEX_Incidence_MSM = RM_INDEX_Treat_MALE_ALL + 1;
				
				int RM_INDEX_Treat_FEMALE_2 = RM_INDEX_Incidence_MSM + 1;
				int RM_INDEX_Treat_MALE_ALL_2 = RM_INDEX_Treat_FEMALE_2 + 1;
				int RM_INDEX_Incidence_MSM_2 = RM_INDEX_Treat_MALE_ALL_2 + 1;
								
				int LENGTH_RM = RM_INDEX_Incidence_MSM_2+1;

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
							int select_row = 15; // T = 5840, or 2015
							int select_row_2 = 22; // T = 8395, or 2022			

							// Site_Specific prevalence
							File[] file_prevalence_site = resDir.listFiles(new FileFilter() {
								@Override
								public boolean accept(File pathname) {
									return Pattern.matches(".*Prevalence_Site_.*", pathname.getName());
								}
							});
							HashMap<String, ArrayList<String[]>> entMap_preval = Util_Select_Sim_By_Residue
									.combineResultMapping(file_prevalence_site);

							for (Entry<String, ArrayList<String[]>> ent : entMap_preval.entrySet()) {
								String map_key = String.format("%s:%s", resDir.getName(), ent.getKey());
								res_entry = getResidueMappingEntry(mapping, map_key);
								ArrayList<String[]> lines = ent.getValue();
								if (select_row < lines.size()) {
									String[] line_current_line = lines.get(select_row);
									res_entry[RM_INDEX_Preval_MSM_NG_U] = 100
											* (Double.parseDouble(line_current_line[10])
													+ Double.parseDouble(line_current_line[14]))
											/ (12940 + 5517);
									res_entry[RM_INDEX_Preval_MSM_NG_R] = 100
											* (Double.parseDouble(line_current_line[11])
													+ Double.parseDouble(line_current_line[15]))
											/ (12940 + 5517);
									res_entry[RM_INDEX_Preval_MSM_NG_P] = 100
											* (Double.parseDouble(line_current_line[12])
													+ Double.parseDouble(line_current_line[16]))
											/ (12940 + 5517);
									res_entry[RM_INDEX_Preval_FEMALE_NG_U] = 100
											* Double.parseDouble(line_current_line[1]) / 500000;
									res_entry[RM_INDEX_Preval_MALE_HETERO_NG_U] = 100
											* Double.parseDouble(line_current_line[6]) / 500000;
								} else {
									// Extinction
									res_entry[RM_INDEX_Preval_MSM_NG_U] = 0;
									res_entry[RM_INDEX_Preval_MSM_NG_R] = 0;
									res_entry[RM_INDEX_Preval_MSM_NG_P] = 0;
									res_entry[RM_INDEX_Preval_FEMALE_NG_U] = 0;
									res_entry[RM_INDEX_Preval_MALE_HETERO_NG_U] = 0;

								}
							}

							// Treatment - Notification
							File[] file_treatment = resDir.listFiles(new FileFilter() {
								@Override
								public boolean accept(File pathname) {
									return Pattern.matches(".*Treatment_Person_.*", pathname.getName());
								}
							});
							HashMap<String, ArrayList<String[]>> entMap_treatment = Util_Select_Sim_By_Residue
									.combineResultMapping(file_treatment);
							for (Entry<String, ArrayList<String[]>> ent : entMap_treatment.entrySet()) {
								String map_key = String.format("%s:%s", resDir.getName(), ent.getKey());
								res_entry = getResidueMappingEntry(mapping, map_key);
								// 2015
								if (select_row < ent.getValue().size()) {
									String[] treatment_line = ent.getValue().get(select_row);
									String[] treatment_last = ent.getValue().get(select_row - 1);
									res_entry[RM_INDEX_Treat_FEMALE] = 100000 * (Double.parseDouble(treatment_line[1])
											- Double.parseDouble(treatment_last[1])) / 500000;
									res_entry[RM_INDEX_Treat_MALE_ALL] = 100000 * (Double.parseDouble(treatment_line[2])
											- Double.parseDouble(treatment_last[2])
											+ Double.parseDouble(treatment_line[3])
											- Double.parseDouble(treatment_last[3])
											+ Double.parseDouble(treatment_line[4])
											- Double.parseDouble(treatment_last[4])) / (500000 + 12940 + 5517);
								} else {
									res_entry[RM_INDEX_Treat_FEMALE] = 0;
									res_entry[RM_INDEX_Treat_MALE_ALL] = 0;
								}
								// 2022
								if (select_row_2 < ent.getValue().size()) {
									String[] treatment_line = ent.getValue().get(select_row_2);
									String[] treatment_last = ent.getValue().get(select_row_2 - 1);
									res_entry[RM_INDEX_Treat_FEMALE_2] = 100000 * (Double.parseDouble(treatment_line[1])
											- Double.parseDouble(treatment_last[1])) / 500000;
									res_entry[RM_INDEX_Treat_MALE_ALL_2] = 100000 * (Double.parseDouble(treatment_line[2])
											- Double.parseDouble(treatment_last[2])
											+ Double.parseDouble(treatment_line[3])
											- Double.parseDouble(treatment_last[3])
											+ Double.parseDouble(treatment_line[4])
											- Double.parseDouble(treatment_last[4])) / (500000 + 12940 + 5517);
								} else {
									res_entry[RM_INDEX_Treat_FEMALE_2] = 0;
									res_entry[RM_INDEX_Treat_MALE_ALL_2] = 0;
								}																
								
							}

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
								if (select_row < ent.getValue().size()) {
									ArrayList<String[]> lines = ent.getValue();
									String[] line_current_line = lines.get(select_row);
									String[] line_past_line = lines.get(select_row - 1);
									res_entry[RM_INDEX_Incidence_MSM] = 100000
											* (Double.parseDouble(line_current_line[3])
													- Double.parseDouble(line_past_line[3])
													+ Double.parseDouble(line_current_line[4])
													- Double.parseDouble(line_past_line[4]))
											/ (12940 + 5517);
								} else {									
									res_entry[RM_INDEX_Incidence_MSM]  = 0;

								}
								
								if (select_row_2 < ent.getValue().size()) {
									ArrayList<String[]> lines = ent.getValue();
									String[] line_current_line = lines.get(select_row_2);
									String[] line_past_line = lines.get(select_row_2 - 1);
									res_entry[RM_INDEX_Incidence_MSM_2] = 100000
											* (Double.parseDouble(line_current_line[3])
													- Double.parseDouble(line_past_line[3])
													+ Double.parseDouble(line_current_line[4])
													- Double.parseDouble(line_past_line[4]))
											/ (12940 + 5517);
								} else {									
									res_entry[RM_INDEX_Incidence_MSM_2]  = 0;

								}

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

			System.out.printf("Parameter list generated at %s.\n", paramListSummaryFile.getAbsolutePath());

		}

		// Resample generate new results

		int genSeedDir_MaxSimPerDir = 16;
		int genSeedDir_num_col_include = 16; // Including map and seed id.
		int genSeedDir_in_range_col = 19;

		float genSeedDir_numSimMin = 1000f;
		long genSeedDir_rng_seed = 22519122707291119l;
		
		int genSeedDir_max_result_to_include = 500; // set to <0 if resample sim are not generated.
		int genSeedDir_num_reseed = 2;
		int genSeedDir_folder_id = 2;
		

		File resampleDir = new File(baseDir, "Extract_Resampled_Seed_Sim");

		String resample_dir_name = String.format("SimClusterModel_Transmission_Extra_%03d", genSeedDir_folder_id); // Adjust
																													// as
																													// needed
		String pattern_PBS_id_line = "#PBS -N SIM_" + genSeedDir_folder_id + "_%03d"; // Adjust as needed

		String pattern_PBS_filename = "ClusterModel_Sim_%03d.pbs";

		File genSeedDir_SIM_DIR = new File(resampleDir, resample_dir_name);
		File genSeedDir_PBS_DIR = new File(resampleDir, "PBS_" + resample_dir_name);
		File genSeedDir_SRC_PBS = new File(
				"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Sim_Template.pbs");
		File genSeedDir_SRC_Template = new File(
				"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Sim");
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
