package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import sim.Simulation_ClusterModelTransmission;

public class Util_Compare_ClusterModel_Transmission_Output {

	File ref_result_dir;
	File results_dir;
	int[][] analysis_col_select;
	File[] compare_dir_arr;

	final static int BUFFER = 2048;
	final static String USAGE_INFO = String.format(
			"Usage: java %s -compare RESULTS_DIRECTORY REF_RESULT_DIRNAME ANALSIS_COL_SELECT",
			Util_Compare_ClusterModel_Transmission_Output.class.getName());

	final static int[] ANALYSIS_TYPE_OPTIONS = new int[] {
			Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE,
			Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE,
			Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE,

	};

	public Util_Compare_ClusterModel_Transmission_Output(File results_dir, File ref_result_dir,
			int[][] analysis_col_select) {
		super();
		this.results_dir = results_dir;
		this.ref_result_dir = ref_result_dir;
		this.analysis_col_select = analysis_col_select;

		this.compare_dir_arr = results_dir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && !pathname.getName().equals(ref_result_dir.getName());
			}
		});

		Arrays.sort(this.compare_dir_arr, new Comparator<File>() {
			@Override
			public int compare(File o1, File o2) {
				return o1.getName().compareTo(o2.getName());
			}

		});
	}

	public void startCompare() throws IOException {

		System.out.printf("Comparing %d result(s) from %s, using %s as reference.\n", compare_dir_arr.length,
				results_dir.getAbsolutePath(), ref_result_dir.getName());

		for (int analysis_type : ANALYSIS_TYPE_OPTIONS) {

			switch (analysis_type) {

			case Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE:
				// Prevalence by person
				printComparedOutput(analysis_type, Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON_ZIP,
						"Compare_Output_Prevalence.csv", "%.1f (%.1f - %.1f)", "%.3f (%.3f - %.3f)", null, false);
				break;
			case Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE:
				// Treatment by person
				printComparedOutput(analysis_type,
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_TREATMENT_PERSON_ZIP,
						"Compare_Output_Treatment.csv", "%.3f (%.3f - %.3f)", "%.3f (%.3f - %.3f)", new int[] { -4 },
						true);

				break;
			case Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE:
				// Antibiotic use
				printComparedOutput(analysis_type,
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_ANTIBIOTIC_USAGE_ZIP,
						"Compare_Output_AntiboticUsage.csv", "%.1f (%.1f - %.1f)", "%.3f (%.3f - %.3f)", null, true);

				break;
			default:
				System.out.printf("Analysis type %d not support yet\n", analysis_type);

			}

		}
		System.out.println("Comparsion completed.");

	}

	private void printComparedOutput(int analysis_type, String zip_filename, String output_filename,
			String output_format_ref, String output_format_cmp, int[] col_denominator_offsets, boolean cumul_diff)
			throws IOException {

		Pattern patten_src_zip_file;
		FileFilter zip_filter;
		File[] ref_zips, comp_zips;
		HashMap<String, ArrayList<String[]>> ref_ent_map, comp_ent_map;

		boolean[] col_sel = null;
		int[] time_ent = null;
		String[][] ref_summary = null;
		String[][][] comp_summary = null;
		Percentile percentile = new Percentile();

		int[] col_sel_index = new int[0]; // Length of zero = all, null = non selected

		if (analysis_type < analysis_col_select.length) {
			col_sel_index = analysis_col_select[analysis_type];
		}

		if (col_sel_index != null) {

			patten_src_zip_file = Pattern.compile(zip_filename.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));

			zip_filter = generateFileFilterByPattern(patten_src_zip_file);
			ref_zips = ref_result_dir.listFiles(zip_filter);

			ref_ent_map = new HashMap<>();

			for (File zipFile : ref_zips) {
				ref_ent_map = extractedLinesFrom7Zip(zipFile, ref_ent_map);
			}

			// Filling ref column
			String[] ref_ent_arr = ref_ent_map.keySet().toArray(new String[ref_ent_map.size()]);
			Arrays.sort(ref_ent_arr);

			double[][][] ref_values = null;

			for (int reI = 0; reI < ref_ent_arr.length; reI++) {
				String ref_ent = ref_ent_arr[reI];
				ArrayList<String[]> lines = ref_ent_map.get(ref_ent);

				// Set up initial array
				if (lines != null) {
					int numCol;
					col_sel = new boolean[lines.get(0).length];
					if (col_sel_index.length != 0) {
						numCol = col_sel_index.length;
						col_sel[0] = false;
						for (int c = 0; c < col_sel_index.length; c++) {
							col_sel[col_sel_index[c]] = true;
						}
					} else {
						Arrays.fill(col_sel, true);
						col_sel[0] = false;
						numCol = col_sel.length - 1;
					}

					if (time_ent == null) {
						// Without header line
						ref_summary = new String[lines.size() - 1][numCol];
						comp_summary = new String[compare_dir_arr.length][lines.size() - 1][numCol];
						time_ent = new int[lines.size() - 1];

						ref_values = new double[numCol][lines.size() - 1][ref_ent_arr.length];

						for (int r = 0; r < time_ent.length; r++) {
							String[] line = lines.get(r + 1);
							time_ent[r] = Integer.parseInt(line[0]);
						}

					}

					double[] pre_row_raw_value = new double[col_sel.length];

					for (int t = 0; t < time_ent.length; t++) {
						String[] line = lines.get(t + 1);
						int cPt = 0;
						for (int b = 1; b < col_sel.length; b++) {
							if (col_sel[b]) {
								double raw_value = Double.parseDouble(line[b]);
								double raw_value_denominator = 1;
								if (col_denominator_offsets != null) {
									raw_value_denominator = 0;
									for (int offset : col_denominator_offsets) {
										raw_value_denominator += Double.parseDouble(line[b + offset]);
									}
								}
								double ref = raw_value;
								double ref_denominator = raw_value_denominator;
								if (cumul_diff && t > 0) {
									ref = ref - pre_row_raw_value[b];
									if (col_denominator_offsets != null) {
										for (int offset : col_denominator_offsets) {
											ref_denominator = ref_denominator - pre_row_raw_value[b + offset];
										}
									}
								}
								if (col_denominator_offsets != null) {
									ref = ref / ref_denominator;
								}
								ref_values[cPt][t][reI] = ref;

								pre_row_raw_value[b] = raw_value;
								if (col_denominator_offsets != null) {
									for (int offset : col_denominator_offsets) {
										pre_row_raw_value[b + offset] = Double.parseDouble(line[b + offset]);
									}
								}
								cPt++;
							}
						}
					}
				}
			}

			for (int c = 0; c < ref_values.length; c++) {
				for (int t = 0; t < ref_values[c].length; t++) {
					percentile.setData(ref_values[c][t]);
					ref_summary[t][c] = String.format(output_format_ref, percentile.evaluate(50),
							percentile.evaluate(25), percentile.evaluate(75));

				}
			}

			// Match with other mapping in compare_result_dir
			for (int f = 0; f < compare_dir_arr.length; f++) {
				File compare_dir = compare_dir_arr[f];
				comp_zips = compare_dir.listFiles(zip_filter);
				comp_ent_map = new HashMap<>();
				for (File zipFile : comp_zips) {
					comp_ent_map = extractedLinesFrom7Zip(zipFile, comp_ent_map);
				}
				double[][][] rel_diff_from_ref_values = new double[ref_values.length][time_ent.length][ref_ent_arr.length];

				for (int reI = 0; reI < ref_ent_arr.length; reI++) {
					String ref_ent = ref_ent_arr[reI];
					ArrayList<String[]> lines = comp_ent_map.get(ref_ent);

					if (lines != null) {
						double[] pre_row_raw_val = new double[col_sel.length];

						for (int t = 0; t < time_ent.length; t++) {
							String[] line = lines.get(t + 1);
							int cPt = 0;
							for (int b = 1; b < col_sel.length; b++) {
								if (col_sel[b]) {
									double ref_val = ref_values[cPt][t][reI];
									double raw_cmp_value = Double.parseDouble(line[b]);
									double raw_cmp_value_denominator = 1;

									if (col_denominator_offsets != null) {
										raw_cmp_value_denominator = 0;
										for (int offset : col_denominator_offsets) {
											raw_cmp_value_denominator += Double.parseDouble(line[b + offset]);
										}
									}
									double cmp_val = raw_cmp_value;
									double cmp_val_denominator = raw_cmp_value_denominator;
									if (cumul_diff && t > 0) {
										cmp_val = cmp_val - pre_row_raw_val[b];
										if (col_denominator_offsets != null) {
											for (int offset : col_denominator_offsets) {
												cmp_val_denominator = cmp_val_denominator - pre_row_raw_val[b + offset];
											}
										}
									}
									if (col_denominator_offsets != null) {
										cmp_val = cmp_val / cmp_val_denominator;
									}
									rel_diff_from_ref_values[cPt][t][reI] = cmp_val / ref_val - 1;

									pre_row_raw_val[b] = raw_cmp_value;
									if (col_denominator_offsets != null) {
										for (int offset : col_denominator_offsets) {
											pre_row_raw_val[b + offset] = Double.parseDouble(line[b + offset]);
										}
									}

									cPt++;
								}
							}
						}
					}
				}

				for (int c = 0; c < ref_values.length; c++) {
					for (int t = 0; t < ref_values[c].length; t++) {
						percentile.setData(rel_diff_from_ref_values[c][t]);
						comp_summary[f][t][c] = String.format(output_format_cmp, percentile.evaluate(50),
								percentile.evaluate(25), percentile.evaluate(75));

					}
				}

				// Print results

				File file_compareOutput = new File(results_dir, output_filename);
				PrintWriter pWri = new PrintWriter(new FileWriter(file_compareOutput));
				// First line
				StringBuilder pLine = new StringBuilder();
				pLine.append("Time");
				for (int c = 0; c < ref_values.length; c++) {
					pLine.append(',');
					pLine.append(ref_result_dir.getName());
					if (ref_values.length > 1) {
						pLine.append("_Col_");
						pLine.append(c);
					}
				}
				for (int cdI = 0; cdI < compare_dir_arr.length; cdI++) {
					for (int c = 0; c < ref_values.length; c++) {
						pLine.append(',');
						pLine.append("Rel_to_");
						pLine.append(compare_dir_arr[cdI].getName());
						if (ref_values.length > 1) {
							pLine.append("_Col_");
							pLine.append(c);
						}
					}
				}
				pWri.println(pLine.toString());

				for (int t = 0; t < time_ent.length; t++) {
					pLine = new StringBuilder();
					pLine.append(time_ent[t]);
					for (int c = 0; c < ref_values.length; c++) {
						pLine.append(',');
						pLine.append(ref_summary[t][c]);
					}
					for (int cdI = 0; cdI < compare_dir_arr.length; cdI++) {
						for (int c = 0; c < ref_values.length; c++) {
							pLine.append(',');
							pLine.append(comp_summary[cdI][t][c]);
						}

					}
					pWri.println(pLine.toString());
				}

				pWri.close();

			}
		}
	}

	private static HashMap<String, ArrayList<String[]>> extractedLinesFrom7Zip(File zipFile,
			HashMap<String, ArrayList<String[]>> zip_ent) throws IOException {
		SevenZFile inputZip = new SevenZFile(zipFile);
		SevenZArchiveEntry inputEnt;

		byte[] buf = new byte[BUFFER];
		while ((inputEnt = inputZip.getNextEntry()) != null) {
			String file_name = inputEnt.getName();
			StringBuilder str_builder = new StringBuilder();
			String line;
			ArrayList<String[]> lines = new ArrayList<>();
			int count;
			while ((count = inputZip.read(buf, 0, BUFFER)) != -1) {
				str_builder.append(new String(Arrays.copyOf(buf, count)));
			}
			BufferedReader reader = new BufferedReader(new StringReader(str_builder.toString()));
			while ((line = reader.readLine()) != null) {
				if (line.length() > 0) {
					lines.add(line.split(","));
				}
			}
			zip_ent.put(file_name, lines);
		}
		inputZip.close();
		return zip_ent;
	}

	private static FileFilter generateFileFilterByPattern(Pattern pattern) {
		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pattern.matcher(pathname.getName()).matches();
			}
		};

		return filter;
	}

	public static void launch(String[] arg) throws IOException {
		if (arg.length < 3) {
			System.out.println(Util_Compare_ClusterModel_Transmission_Output.USAGE_INFO);
			System.exit(1);
		} else {
			File results_dir = new File(arg[0]);
			File ref_result_dir = new File(results_dir, arg[1]);
			int[][] analysis_col_select = (int[][]) PropValUtils.propStrToObject(arg[2], int[][].class);

			Util_Compare_ClusterModel_Transmission_Output cmp = new Util_Compare_ClusterModel_Transmission_Output(
					results_dir, ref_result_dir, analysis_col_select);

			cmp.startCompare();
		}

	}

}
