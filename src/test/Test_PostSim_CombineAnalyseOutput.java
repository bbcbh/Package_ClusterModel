package test;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import sim.Simulation_ClusterModelTransmission;
import util.Util_7Z_CSV_Entry_Extract_Callable;
import util.Util_Analyse_ClusterModel_Transmission_Output_Combined;
import util.Util_SimulationDirectoryModifications;

public class Test_PostSim_CombineAnalyseOutput {

	public static void main(String[] args) throws IOException {
		File tarDirPathBase, extraDir;
		String prefix_dir;
		Pattern tarDirMatch;

		String[][] post_zip_patterns = new String[][] { new String[] {
				Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("_%d", "") + ".7z",
				Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("%d",
						"(-{0,1}\\\\d+)"), },
				new String[] { "Seed_List.csv.zip", "Seed_List(_\\d+).csv", },
				new String[] { "simSpecificSwitch.prop.zip", "simSpecificSwitch(_\\d+).prop", },
				new String[] { "dx_Bali.prop.zip", "dx_Bali(_\\d+).prop", }, };

		final int PROP_TYPE_SIM = 0;
		final int PROP_TYPE_SIM_DOXY_PEP = PROP_TYPE_SIM + 1;
		final int PROP_TYPE_SIM_VIABILITY = PROP_TYPE_SIM_DOXY_PEP + 1;
		final int PROP_TYPE_SIM_BALI = PROP_TYPE_SIM_VIABILITY + 1;

		int propType = PROP_TYPE_SIM_VIABILITY;

		Util_Analyse_ClusterModel_Transmission_Output_Combined analysis;

		switch (propType) {
		case PROP_TYPE_SIM_DOXY_PEP:
			tarDirPathBase = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Result_DoxyPEP\\Resist_Buildup_Eligiblity");
			tarDirMatch = Pattern.compile("MultiTransmission_MSM_S.*");
			prefix_dir = "SimClusterModel_Transmission";
			String replace_string_DoxyPEP = "(-?\\\\d+(?:_-?\\\\d+){0,1})";
			String[] zipFileFormats_DoxyPEP = new String[] {
					Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON_ZIP.replaceAll("%d",
							replace_string_DoxyPEP),
					"Infectious_Prevalence_Site_%d.csv.7z".replaceAll("%d", replace_string_DoxyPEP),
					"Infectious_Prevalence_Person_%d.csv.7z".replaceAll("%d", replace_string_DoxyPEP),
					"Infected_Prevalence_Site_%d.csv.7z".replaceAll("%d", replace_string_DoxyPEP),
					"Treatment_Person_%d.csv.7z".replaceAll("%d", replace_string_DoxyPEP),
					"PEP_Stat_%d.csv.7z".replaceAll("%d", replace_string_DoxyPEP) };
			String[] summaryFileFormat_DoxyPEP = new String[] { "Summary_Incidence_Person_%s.csv",
					"Summary_Infectious_Prevalence_Site_%s.csv", "Summary_Infectious_Prevalence_Person_%s.csv",
					"Summary_Infected_Prevalence_Site_%s.csv", "Summary_Treatment_Stat_%s.csv",
					"Summary_PEP_Stat_%s.csv" };
			boolean[] isCumulData_DoxyPEP = new boolean[] { true, false, false, false, true, true };
			analysis = new Util_Analyse_ClusterModel_Transmission_Output_Combined(zipFileFormats_DoxyPEP,
					summaryFileFormat_DoxyPEP, isCumulData_DoxyPEP);

			// PROP_TYPE_SIM_DOXY_PEP specific analysis

			File[] tarDirs_DOXY_PEP = new File[] { tarDirPathBase };

			if (tarDirMatch != null) {
				tarDirs_DOXY_PEP = tarDirPathBase.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isDirectory() && tarDirMatch.matcher(pathname.getName()).matches();
					}
				});
			}

			for (int fI = 0; fI < tarDirs_DOXY_PEP.length; fI++) {
				tarDirs_DOXY_PEP[fI] = new File(tarDirs_DOXY_PEP[fI], "SimClusterModel_Transmission");
			}

			// Resist Profile

			int[][] resist_col_index = new int[][] { //
					new int[] { 13, 18, 19, 20, 22, 23, 24 }, // Doxy Sensitive
					new int[] { 25, 30, 31, 32, 34, 35, 36 }, // Doxy Resist
			};

			int[][] resist_grouping_index = new int[][] { new int[] { 0 }, new int[] { 1, 2, 3 },
					new int[] { 4, 5, 6 }, };

			String resist_summary_fileName = "Summary_PEP_Sensitivity";

			for (File dir_Doxy_PEP : tarDirs_DOXY_PEP) {

				File[] zip_resist_profile = dir_Doxy_PEP.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.getName().startsWith("PEP_Resist_Profile_");
					}
				});

				HashMap<Integer, HashMap<Integer, ArrayList<Double>>> resist_mapping = new HashMap<>();
				StringBuilder resist_summary_header = new StringBuilder();
				for (File resistFileZip : zip_resist_profile) {
					HashMap<String, ArrayList<String[]>> resMap = util.Util_7Z_CSV_Entry_Extract_Callable
							.extractedLinesFrom7Zip(resistFileZip);
					for (Entry<String, ArrayList<String[]>> ent : resMap.entrySet()) {
						ArrayList<String[]> lines = ent.getValue();
						int lineNum = 0;
						for (String[] line : lines) {
							if (lineNum == 0) {
								if (resist_summary_header.length() == 0) {
									resist_summary_header.append("Time");
									for (int i = 0; i < resist_col_index[0].length; i++) {
										resist_summary_header.append(",");
										resist_summary_header.append(String.format(
												"S%3$d = %1$s / (%1$s + %2$s), S%3$d_Q1, S%3$d_Q3, S%3$d_Disp",
												line[resist_col_index[0][i]], line[resist_col_index[1][i]], i));
									}
									for (int i = 0; i < resist_grouping_index.length; i++) {
										resist_summary_header.append(
												String.format(",Grp_%1$d, Grp_%1$d_Q1, Grp_%1$d_Q3, Grp_%1$d_Disp", i));
									}

								}
							} else {
								int time = Integer.parseInt(line[0]);
								HashMap<Integer, ArrayList<Double>> resist_map = resist_mapping.get(time);

								if (resist_map == null) {
									resist_map = new HashMap<>();
									resist_mapping.put(time, resist_map);
								}
								for (int i = 0; i < resist_col_index[0].length; i++) {
									ArrayList<Double> resist_map_ent = resist_map.get(i);
									if (resist_map_ent == null) {
										resist_map_ent = new ArrayList<>();
										resist_map.put(i, resist_map_ent);
									}
									double rS = Double.parseDouble(line[resist_col_index[0][i]]);
									double rR = Double.parseDouble(line[resist_col_index[1][i]]);
									resist_map_ent.add(100.0 * rS / (rS + rR));
								}

								for (int i = 0; i < resist_grouping_index.length; i++) {
									ArrayList<Double> resist_map_ent = resist_map.get(i + resist_col_index[0].length);
									if (resist_map_ent == null) {
										resist_map_ent = new ArrayList<>();
										resist_map.put(i + resist_col_index[0].length, resist_map_ent);
									}
									double rS_T = 0;
									double rR_T = 0;
									for (int g = 0; g < resist_grouping_index[i].length; g++) {
										rS_T += Double
												.parseDouble(line[resist_col_index[0][resist_grouping_index[i][g]]]);
										rR_T += Double
												.parseDouble(line[resist_col_index[1][resist_grouping_index[i][g]]]);
									}

									resist_map_ent.add(100.0 * rS_T / (rS_T + rR_T));
								}

							}
							lineNum++;
						}
					}

					// Print resistance summary

					if (!resist_mapping.isEmpty()) {
						Percentile percentile = new Percentile();
						File resist_summaryFile = new File(dir_Doxy_PEP, resist_summary_fileName);
						resist_summaryFile.mkdirs();
						resist_summaryFile = new File(resist_summaryFile, resist_summary_fileName + ".csv");
						PrintWriter pWri = new PrintWriter(resist_summaryFile);
						pWri.println(resist_summary_header.toString());

						Integer[] timeArr = resist_mapping.keySet().toArray(new Integer[0]);
						Arrays.sort(timeArr);

						for (Integer rTime : timeArr) {
							pWri.print(rTime);
							HashMap<Integer, ArrayList<Double>> resist_map_rows = resist_mapping.get(rTime);
							for (int i = 0; i < resist_col_index[0].length + resist_grouping_index.length; i++) {
								ArrayList<Double> ent = resist_map_rows.get(i);
								Double[] ent_Double = ent.toArray(new Double[0]);
								double[] ent_arr = new double[ent.size()];
								for (int s = 0; s < ent_arr.length; s++) {
									ent_arr[s] = ent_Double[s];
								}
								// Median,Q1 Q3, Disp
								percentile.setData(ent_arr);
								pWri.printf(",%1$f,%2$f,%3$f,%1$.1f (%2$.1f - %3$.1f)", percentile.evaluate(50),
										percentile.evaluate(25), percentile.evaluate(75));
							}
							pWri.println();
						}
						pWri.close();

					}

				}

			}

			// System.out.println("All done");
			// System.exit(1);

			break;
		case PROP_TYPE_SIM_VIABILITY:
			tarDirPathBase = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_Baseline_Sel_2000_Simple");
			// tarDirPathBase = new
			// File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\All\\MSM_Viability_Baseline");

			tarDirMatch = null;
			prefix_dir = "SimClusterModel_Transmission";

			String replace_string_Viability = "(-?\\\\d+(?:_-?\\\\d+){0,1})";

			String[] zipFileFormats_Viability = new String[] {
					Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON_ZIP.replaceAll("%d",
							replace_string_Viability),
					"Infectious_Prevalence_Site_%d.csv.7z".replaceAll("%d", replace_string_Viability),
					"Infectious_Prevalence_Person_%d.csv.7z".replaceAll("%d", replace_string_Viability),
					"Infected_Prevalence_Site_%d.csv.7z".replaceAll("%d", replace_string_Viability),
					"Treatment_Person_%d.csv.7z".replaceAll("%d", replace_string_Viability),
					"Treatment_Non_Infected_Site_%d.csv.7z".replaceAll("%d", replace_string_Viability),
					"Treatment_Non_Viable_Site_%d.csv.7z".replaceAll("%d", replace_string_Viability), };
			String[] summaryFileFormat_Viability = new String[] { "Summary_Incidence_Person_%s.csv",
					"Summary_Infectious_Prevalence_Site_%s.csv", "Summary_Infectious_Prevalence_Person_%s.csv",
					"Summary_Infected_Prevalence_Site_%s.csv", "Summary_Treatment_Stat_%s.csv",
					"Summary_Non_Infected_Site_%s.csv", "Summary_Treatment_Non_Viable_Site_%s.csv", };
			boolean[] isCumulData_Viability = new boolean[] { true, false, false, false, true, true, true };
			analysis = new Util_Analyse_ClusterModel_Transmission_Output_Combined(zipFileFormats_Viability,
					summaryFileFormat_Viability, isCumulData_Viability);

			break;
		case PROP_TYPE_SIM_BALI:
			tarDirPathBase = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Results_Bali\\");
			tarDirMatch = Pattern.compile("Bali_PREP.*");
			prefix_dir = "SimClusterModel_Transmission";
			String replace_string = "(-?\\\\d+(?:_-?\\\\d+){0,1})";

			String[] zipFileFormats = new String[] {
					Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON_ZIP.replaceAll("%d",
							replace_string),
					"Infectious_Prevalence_Site_%d.csv.7z".replaceAll("%d", replace_string),
					"Infectious_Prevalence_Person_%d.csv.7z".replaceAll("%d", replace_string),
					"Infected_Prevalence_Site_%d.csv.7z".replaceAll("%d", replace_string),
					"PEP_Stat_%d.csv.7z".replaceAll("%d", replace_string), };
			String[] summaryFileFormat = new String[] { "Summary_Incidence_Person_%s.csv",
					"Summary_Infectious_Prevalence_Site_%s.csv", "Summary_Infectious_Prevalence_Person_%s.csv",
					"Summary_Infected_Prevalence_Site_%s.csv", "Summary_PEP_Stat_%s.csv", };
			boolean[] isCumulData = new boolean[] { true, false, false, false, false };
			
			
			analysis = new Util_Analyse_ClusterModel_Transmission_Output_Combined(zipFileFormats, summaryFileFormat,
					isCumulData);

			break;
		case PROP_TYPE_SIM:
		default:
			tarDirPathBase = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Results_Bridging\\" + "Baseline_Split_4");
			tarDirMatch = null;
			prefix_dir = "SimClusterModel_Transmission";
			analysis = new Util_Analyse_ClusterModel_Transmission_Output_Combined();
			analysis.setSkipAnalysis(712);
		}

		File[] targetDirs = new File[] { tarDirPathBase };

		if (tarDirMatch != null) {
			targetDirs = tarDirPathBase.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.isDirectory() && tarDirMatch.matcher(pathname.getName()).matches();
				}
			});
		}

		for (File tarDirPath : targetDirs) {

			Pattern prefix_pattern = Pattern.compile(String.format("%s_Extra_(\\d+)", prefix_dir));
			extraDir = new File(tarDirPath, "Extra");
			extraDir.mkdirs();

			File[] dirSrc = tarDirPath.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.isDirectory() && pathname.getName().startsWith(prefix_dir);
				}
			});

			Arrays.sort(dirSrc, new Comparator<File>() {
				@Override
				public int compare(File o1, File o2) {
					Matcher m1 = prefix_pattern.matcher(o1.getName());
					Matcher m2 = prefix_pattern.matcher(o2.getName());
					m1.matches();
					m2.matches();
					return Integer.compare(Integer.parseInt(m1.group(1)), Integer.parseInt(m2.group(1)));
				}
			});

			if (dirSrc.length > 0) {
				File dir = new File(tarDirPath, prefix_dir);

				if (analysis != null) {

					// Pre 20241023
					Files.move(dirSrc[0].toPath(), new File(tarDirPath, prefix_dir).toPath(),
							StandardCopyOption.ATOMIC_MOVE);
					for (int f = 1; f < dirSrc.length; f++) {
						Files.move(dirSrc[f].toPath(), new File(extraDir, dirSrc[f].getName()).toPath(),
								StandardCopyOption.ATOMIC_MOVE);
					}
					Util_SimulationDirectoryModifications.combineSimOutput(tarDirPath, extraDir);

					System.out.println("Combine sim output completed.");

					
					analysis.setBaseDir(dir);

					System.out.printf("=== %s ===\n", dir.getName());

					analysis.analyse_outputs();
				}

				// Post analysis
				switch (propType) {
				case PROP_TYPE_SIM_VIABILITY:

					String replace_string_Viability = "(-?\\\\d+(?:_-?\\\\d+){0,1})";

					Pattern pattern_treatment = Pattern
							.compile("Treatment_Person_%d.csv.7z".replaceAll("%d", replace_string_Viability));
					Pattern pattern_nv_site = Pattern
							.compile("Treatment_Non_Viable_Site_%d.csv.7z".replaceAll("%d", replace_string_Viability));

					HashMap<String, int[][]> treatment_store = new HashMap<>();
					// Treatment file extract
					File[] treatFileArr = dir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return pattern_treatment.matcher(pathname.getName()).matches();
						}
					});

					for (File treatFile : treatFileArr) {
						Matcher m = pattern_treatment.matcher(treatFile.getName());
						m.find();
						String key_prefix = m.group(1);
						HashMap<String, ArrayList<String[]>> entMap = Util_7Z_CSV_Entry_Extract_Callable
								.extractedLinesFrom7Zip(treatFile);

						for (String keys : entMap.keySet()) {
							Matcher seedMatcher = Pattern.compile("\\[(.*)\\].*").matcher(keys);
							seedMatcher.find();
							String seedId = seedMatcher.group(1);

							ArrayList<String[]> rowsEnt = entMap.get(keys);
							int[][] entArr = null;
							int rowId = 0;
							for (String[] row : rowsEnt) {
								if (entArr == null) {
									entArr = new int[row.length][rowsEnt.size()];
								}
								if (rowId != 0) {
									for (int s = 0; s < row.length; s++) {
										entArr[s][rowId] = Integer.parseInt(row[s]);
									}
								}
								rowId++;
							}
							String key = String.format("%s_%s", key_prefix, seedId);
							if (treatment_store.containsKey(key)) {
								System.err.println("Warning. Entry already exist and will be overwritten.");
							}
							treatment_store.put(key, entArr);
						}
					}
					System.out.printf("Treatment info from %d simulations extracted.\n", treatment_store.size());
					// Non-viable infection extract

					int[] col_sel = new int[] { 2, 3, 4, 5, 6, 7 };
					int[] treat_col = new int[] { 2, 2, 2, 3, 3, 3 };

					File[] nvFileArr = dir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return pattern_nv_site.matcher(pathname.getName()).matches();
						}
					});

					PrintWriter[] nv_pri = new PrintWriter[col_sel.length];

					File nv_stat_dir = new File(dir, "Non_Viable_DX_All");
					nv_stat_dir.mkdirs();
					for (int i = 0; i < nv_pri.length; i++) {
						nv_pri[i] = new PrintWriter(
								new File(nv_stat_dir, String.format("NV_Stat_Col_%d.csv", col_sel[i])));
						nv_pri[i].println("Sim_Id,% NV DX at Time");
						nv_pri[i].println();
					}

					for (File nvFile : nvFileArr) {
						Matcher m = pattern_nv_site.matcher(nvFile.getName());
						m.find();
						String key_prefix = m.group(1);
						HashMap<String, ArrayList<String[]>> entMap = Util_7Z_CSV_Entry_Extract_Callable
								.extractedLinesFrom7Zip(nvFile);
						for (String keys : entMap.keySet()) {
							Matcher seedMatcher = Pattern.compile("\\[(.*)\\].*").matcher(keys);
							seedMatcher.find();
							String seedId = seedMatcher.group(1);
							String key = String.format("%s_%s", key_prefix, seedId);

							if (!treatment_store.containsKey(key)) {
								System.err.println("Warning. Treatment entry not exist and will be skipped.");
							} else {
								int[][] cumul_treatment = treatment_store.get(key);
								ArrayList<String[]> rowsEnt = entMap.get(keys);
								int[] pre_nw_rol = new int[col_sel.length];
								Arrays.fill(pre_nw_rol, -1);

								int tI = 0;
								for (String[] row : rowsEnt) {
									if (tI > 0) {
										for (int cI = 0; cI < col_sel.length; cI++) {
											int cumul_val = Integer.parseInt(row[col_sel[cI]]);
											if (pre_nw_rol[cI] == -1) {
												nv_pri[cI].print(key.replaceAll(",", "_"));
											} else {
												int num_nv = cumul_val - pre_nw_rol[cI];
												int[] cumul_treatment_base = cumul_treatment[treat_col[cI]];
												int num_tr = cumul_treatment_base[tI] - cumul_treatment_base[tI - 1];
												nv_pri[cI].print(',');
												nv_pri[cI].print(100f * num_nv / num_tr);
											}
											pre_nw_rol[cI] = cumul_val;
										}
									}
									tI++;
								}
								for (int cI = 0; cI < col_sel.length; cI++) {
									nv_pri[cI].println();
								}

							}
						}

					}

					for (int i = 0; i < nv_pri.length; i++) {
						nv_pri[i].close();
					}
					System.out.printf("Non-viable DX from %d simulations extracted at %s.\n", treatment_store.size(),
							nv_stat_dir.getAbsolutePath());

					break;
				}

			}

			Files.walk(extraDir.toPath()).sorted(Comparator.reverseOrder()).forEach(path -> {
				try {
					Files.delete(path);
				} catch (IOException e) {
					e.printStackTrace(System.err);
				}

			});

			System.out.printf("==== Combine and analysis completed for %s. ====\n", tarDirPath.getName());
		}

	}

}
