package test;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sim.Simulation_ClusterModelTransmission;
import util.Util_Analysis_Compare_Cumulative;
import util.Util_SimulationDirectoryModifications;

public class Test_PostSim_Analysis_Viability_Test_Performance {
	public static void main(String[] args) throws IOException {
		File baseDir = new File(
				"C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_Baseline_Sel_Seed_2000");

		File[] compareDirs = new File[] { 
				new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_VT_Mix_Sel_Seed_2000_00"),
				new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_VT_Mix_Sel_Seed_2000_01"),
				new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_VT_Mix_Sel_Seed_2000_02"),
				new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_VT_Mix_Sel_Seed_2000_03"),
				new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_VT_Mix_Sel_Seed_2000_04"),
				new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\MSM_Viability_VT_Mix_Sel_Seed_2000_05"), };

		String[] compareFilePrefix = new String[] { "Treatment_Non_Viable_Site", "Treatment_Person",
				"Incidence_Person" };
		File dataDumpCSV_dir = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Viability\\Compare_Results");

		String[] time_to_match = new String[] { "6935", "8760", "10585" };

		Util_SimulationDirectoryModifications.generateUniqueZip(baseDir,
				Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("_%d", "") + ".7z",
				Pattern.compile(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("%d",
						"(-{0,1}\\\\d+)")));

		for (File f : compareDirs) {
			Util_SimulationDirectoryModifications.generateUniqueZip(f,
					Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("_%d", "") + ".7z",
					Pattern.compile(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP.replaceAll("%d",
							"(-{0,1}\\\\d+)")));
		}

		final String[] extra_header = new String[] { //
				// Shared dir output
				"NG_Adj_to_Viabile_Detection", "NG_Adj_to_Non_Viable_Detection", "CT_Adj_to_Viabile_Detection",
				"CT_Adj_to_Non_Viable_Detection",
				// Share simulation output
				"NG_Non_Viable_Percentage_Start_S1", "CT_Non_Viable_Percentage_Start_S1",
				"NG_Non_Viable_Percentage_Start_S2", "CT_Non_Viable_Percentage_Start_S2",
				"NG_Non_Viable_Percentage_Start_S3", "CT_Non_Viable_Percentage_Start_S3",
				"NG_Non_Viable_Percentage_Since_S1", "CT_Non_Viable_Percentage_Since_S1",
				"NG_Non_Viable_Percentage_Since_S2", "CT_Non_Viable_Percentage_Since_S2",
				"NG_Non_Viable_Percentage_Since_S3", "CT_Non_Viable_Percentage_Since_S3",

		};

		Util_Analysis_Compare_Cumulative analysis = new Util_Analysis_Compare_Cumulative(baseDir, compareDirs,
				compareFilePrefix, time_to_match, dataDumpCSV_dir) {

			// Viability linked entry setting
			@Override
			protected String[] generateLinkedEntryHeader() {
				return extra_header;
			}

			@Override
			protected HashMap<String, String[]> generateLinkedEntryMap(File compDir,
					HashMap<String, String[]> comp_entry_map) {
				HashMap<String, String[]> linked_entry_map = new HashMap<>();
				File[] subDirs = compDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isDirectory();
					}
				});

				Pattern line_pattern = Pattern.compile(",\\[\\d+,1,([\\d\\.]+),2,-2,([\\d\\.]+),12,-2\\].*");
				float[] defaultVal = new float[] { 0.99f, 0.99f, 0.99f, 0.99f };
				for (File subDir : subDirs) {
					File switchFile = new File(subDir, sim.Simulation_ClusterModelTransmission.FILENAME_PROP_SWITCH);
					try {
						String[] lines = util.Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(switchFile);
						String entryMapKey = String.format("%s", subDir.getName());
						String[] entryMapEnt = new String[4];

						// Extract entry from switch file
						// NG: lineNum = 54; CT: lineNum = 55;
						int strIndex = 0;
						for (int lineNum : new int[] { 54, 55 }) {
							Matcher line_mat = line_pattern.matcher(lines[lineNum]);
							if (line_mat.matches()) {
								for (int index : new int[] { 0, 1 }) {
									float adj_val = Float.parseFloat(line_mat.group(index + 1));
									entryMapEnt[strIndex + index] = Float
											.toString(adj_val / defaultVal[strIndex + index]);
								}
							}
							strIndex += 2;
						}
						linked_entry_map.put(entryMapKey, entryMapEnt);

						// Baseline non-viablity

						for (Entry<String, String[]> set : comp_entry_map.entrySet()) {
							String[] map_key_sp = set.getKey().split(":");
							if (map_key_sp[0].equals("0") && !map_key_sp[3].equals(time_to_match[0])) {
								Matcher m = Pattern.compile("\\[(.*)\\](.*)_.*_.*.csv").matcher(map_key_sp[2]);
								m.matches();

								String fileName_alterative = new StringBuilder(map_key_sp[2])
										.replace(m.start(2), m.end(2), compareFilePrefix[1]).toString();
								String fileName_common = m.group(1);

								entryMapKey = String.format("%s:%s:%s", map_key_sp[1], fileName_common, map_key_sp[3]);

								String[] at_start_non_viable_ent = comp_entry_map.get(String.format("%s:%s:%s:%s", "0",
										map_key_sp[1], map_key_sp[2], time_to_match[0]));
								String[] at_start_treatment_ent = comp_entry_map.get(String.format("%s:%s:%s:%s", "1",
										map_key_sp[1], fileName_alterative, time_to_match[0]));

								String[] current_non_viable_ent = comp_entry_map.get(
										String.format("%s:%s:%s:%s", "0", map_key_sp[1], map_key_sp[2], map_key_sp[3]));
								String[] current_treatment_ent = comp_entry_map.get(String.format("%s:%s:%s:%s", "1",
										map_key_sp[1], fileName_alterative, map_key_sp[3]));

								// Share simulation output
								// "NG_Non_Viable_Percentage_Start_S1", "CT_Non_Viable_Percentage_Start_S1",
								// "NG_Non_Viable_Percentage_Start_S2", "CT_Non_Viable_Percentage_Start_S2",
								// "NG_Non_Viable_Percentage_Start_S3", "CT_Non_Viable_Percentage_Start_S3",
								// "NG_Non_Viable_Percentage_At_S1", "CT_Non_Viable_Percentage_At_S1",
								// "NG_Non_Viable_Percentage_At_S2",
								// "CT_Non_Viable_Percentage_At_S2", "NG_Non_Viable_Percentage_At_S3",
								// "CT_Non_Viable_Percentage_At_S3",

								String[] ent = new String[2 * 3 * 2];

								for (int inf_pt : new int[] { 0, 1 }) {
									float at_start_treatment_count = Float
											.parseFloat(at_start_treatment_ent[inf_pt + 2]);
									float current_treamtent_count = Float.parseFloat(current_treatment_ent[inf_pt + 2])
											- at_start_treatment_count;
									for (int site_pt : new int[] { 0, 1, 2 }) {
										float at_start_non_viablible_by_site = Float
												.parseFloat(at_start_non_viable_ent[inf_pt * 3 + site_pt + 2]);
										float current_non_viablible_by_site = Float
												.parseFloat(current_non_viable_ent[inf_pt * 3 + site_pt + 2])
												- at_start_non_viablible_by_site;

										ent[site_pt * 2 + inf_pt] = Float.toString(
												100f * at_start_non_viablible_by_site / at_start_treatment_count);
										ent[6 + site_pt * 2 + inf_pt] = Float.toString(
												100f * current_non_viablible_by_site / current_treamtent_count);

									}

								}

								linked_entry_map.put(entryMapKey, ent);

							}
						}

					} catch (IOException ex) {
						ex.printStackTrace(System.err);
					}

				}
				return linked_entry_map;
			}

			@Override
			protected String[] getLinkedEntry(HashMap<String, String[]> linked_entry_map, String map_key) {
				String[] mapKeySp = map_key.split(":");
				String[] shareDirOutput = linked_entry_map.get(String.format("%s", mapKeySp[1]));

				Matcher m = Pattern.compile("\\[(.*)\\].*").matcher(mapKeySp[2]);
				m.find();

				String[] shareSimOutput = linked_entry_map
						.get(String.format("%s:%s:%s", mapKeySp[1], m.group(1), mapKeySp[3]));

				String[] output = new String[extra_header.length];
				System.arraycopy(shareDirOutput, 0, output, 0, shareDirOutput.length);
				System.arraycopy(shareSimOutput, 0, output, shareDirOutput.length, shareSimOutput.length);

				return output;
			}

		};
		analysis.generateAnalysisCSV();

	}

}
