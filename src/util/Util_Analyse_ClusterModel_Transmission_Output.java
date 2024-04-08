package util;

import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;
import org.apache.commons.compress.archivers.sevenz.SevenZOutputFile;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import optimisation.OptTrendFittingFunction;
import population.Population_Bridging;
import sim.Abstract_Runnable_ClusterModel_Transmission;
import sim.SimulationInterface;
import sim.Simulation_ClusterModelTransmission;

public class Util_Analyse_ClusterModel_Transmission_Output {

	private File baseDir;
	private int[] incl_range = new int[] { 2920, 4745 }; // 5 years

	private static final String replace_string = "(-?\\\\d+(?:_-?\\\\d+){0,1})";

	public static final String[] ZIP_FILES_LIST = new String[] {									
			// 0-3
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_TREATMENT_PERSON_ZIP.replaceAll("%d", replace_string),
			Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON_ZIP.replaceAll("%d", replace_string),
			Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE_ZIP.replaceAll("%d", replace_string),
			Simulation_ClusterModelTransmission.FILENAME_INFECTION_HISTORY_ZIP.replaceAll("%d", replace_string),
			// 4-7
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_PERSON_ZIP.replaceAll("%d", replace_string),
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON_ZIP.replaceAll("%d",
					replace_string),
			Simulation_ClusterModelTransmission.FILENAME_VACCINE_COVERAGE_PERSON_ZIP.replaceAll("%d", replace_string),
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_ANTIBIOTIC_USAGE_ZIP.replaceAll("%d", replace_string),
			// 8-9
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON_ZIP.replaceAll("%d", replace_string),
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_SITE_ZIP.replaceAll("%d", replace_string), };

	public static final String[] STAT_FILEFORMAT = new String[] { "Summary_Treatment_Person_%s.csv",
			"Summary_Prevalence_Person_%s.csv", "Summary_Prevalence_Site_%s.csv", "Summary_Infection_History_%s.csv",
			"Summary_DX_Person_%s.csv", "Summary_DX_Sought_Person_%s.csv", "Summary_Vaccine_Person_%s.csv",
			"Summary_Cumul_Antibiotic_Usage_%s.csv", "Summary_Incidence_Person_%s.csv",
			"Summary_Incidence_Site_%s.csv" };

	public static final boolean[] CUMUL_DATA = new boolean[] { true, false, false, false, true, true, false, false,
			true, true };

	public static final boolean[] SKIP_ANALYSIS = new boolean[] { false, false, false, false, false, false, false,
			false, false, false };

	public final static int rISK_GRP_MAP_INDEX_PID = 0;
	public final static int rISK_GRP_MAP_INDEX_NUM_TIME_SPAN = rISK_GRP_MAP_INDEX_PID + 1;
	public final static int rISK_GRP_MAP_INDEX_NUM_INC = rISK_GRP_MAP_INDEX_NUM_TIME_SPAN + 1;
	public final static int rISK_GRP_MAP_INDEX_NUM_NOTIF = rISK_GRP_MAP_INDEX_NUM_INC + 1;
	public final static int rISK_GRP_MAP_LENGTH = rISK_GRP_MAP_INDEX_NUM_NOTIF + 1;

	public Util_Analyse_ClusterModel_Transmission_Output() {

	}

	public void setSkipAnalysis(int key) {
		// Ox1 = skip 1st, Ox11 = skip 1st and 2nd etc.
		for (int i = 0; i < SKIP_ANALYSIS.length; i++) {
			SKIP_ANALYSIS[i] = (key & (1 << i)) != 0;
		}
	}

	public void setIncl_range(int[] ent) {
		incl_range = ent;
	}

	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;
	}

	public void analyse_outputs() throws IOException {
		final int BUFFER = 2048;
		final byte[] buf = new byte[BUFFER];
		File propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);

		if (propFile.isFile()) {
			for (int z = 0; z < ZIP_FILES_LIST.length; z++) {
				String zipFileName = ZIP_FILES_LIST[z];
				String stat_filename_format = STAT_FILEFORMAT[z];
				boolean isCumul = CUMUL_DATA[z];

				// Check if there is any unzipped csv file
				String csvFileName = String.format("\\[.*\\]%s", zipFileName.replaceFirst(".7z", ""));
				Pattern csvFilePattern = Pattern.compile(csvFileName);

				File[] csvFiles = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return csvFilePattern.matcher(pathname.getName()).matches();
					}
				});

				if (csvFiles.length > 0) {

					Matcher file_prefix = Pattern.compile("(.*)_\\(.*\\)(.csv.7z)").matcher(zipFileName);
					file_prefix.matches();
					
					
					File zipFileSingle = new File(baseDir,
							String.format("%s_0%s", file_prefix.group(1), file_prefix.group(2)));

					Util_7Z_CSV_Entry_Extract_Callable.zipFile(csvFiles, zipFileSingle);

				}

				final Pattern pattern_zip = Pattern.compile(zipFileName);

				File[] zipFiles = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pattern_zip.matcher(pathname.getName()).matches();
					}
				});

				if (zipFiles.length == 0) {
					// Try to zip existing file
					Pattern possible_csv = Pattern
							.compile("(?:.*){0,1}" + zipFileName.replaceAll("\\.csv\\.7z", "_(-{0,1}\\\\d+)\\.csv"));

					File[] raw_csvs = baseDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return possible_csv.matcher(pathname.getName()).matches();
						}
					});

					if (raw_csvs.length > 0) {
						ArrayList<Long> cMapSeedArr = new ArrayList<>();

						for (int fI = 0; fI < raw_csvs.length; fI++) {
							// Get cMap_seed
							Matcher m = possible_csv.matcher(raw_csvs[0].getName());
							m.matches();
							long cMap_seed = Long.parseLong(m.group(1));
							int k = Collections.binarySearch(cMapSeedArr, cMap_seed);
							if (k < 0) {
								cMapSeedArr.add(cMap_seed);
							}
						}

						zipFiles = new File[cMapSeedArr.size()];
						int pt = 0;

						for (long cMap_seed : cMapSeedArr) {

							String newZipFilename = zipFileName.replaceFirst("\\(.*\\)", Long.toString(cMap_seed));
							Pattern possible_csv_cMap_seed = Pattern
									.compile(newZipFilename.replaceAll("\\.csv\\.7z", "_(-{0,1}\\\\d+)\\.csv"));

							File[] raw_csvs_cMap_seed = baseDir.listFiles(new FileFilter() {
								@Override
								public boolean accept(File pathname) {
									return possible_csv_cMap_seed.matcher(pathname.getName()).matches();
								}
							});

							zipFiles[pt] = new File(baseDir, newZipFilename);
							SevenZOutputFile outputZip = new SevenZOutputFile(zipFiles[pt]);

							SevenZArchiveEntry entry;
							FileInputStream fIn;

							for (int fI = 0; fI < raw_csvs_cMap_seed.length; fI++) {
								entry = outputZip.createArchiveEntry(raw_csvs_cMap_seed[fI],
										raw_csvs_cMap_seed[fI].getName());
								outputZip.putArchiveEntry(entry);
								fIn = new FileInputStream(raw_csvs_cMap_seed[fI]);
								outputZip.write(fIn);
								outputZip.closeArchiveEntry();
								fIn.close();
							}

							outputZip.close();

							for (File f : raw_csvs_cMap_seed) {
								f.delete();
							}

							pt++;

						}

					}

				}

				if (SKIP_ANALYSIS[z]) {
					zipFiles = new File[0];
				}

				if (zipFiles.length > 0) {
					System.out.printf("Analysing %d file(s) of format \"%s\".\n", zipFiles.length, zipFileName);

					if (zipFileName.equals(Simulation_ClusterModelTransmission.FILENAME_INFECTION_HISTORY_ZIP
							.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"))) {
						// Special case for infection history

						Properties prop = new Properties();
						FileInputStream fIS = new FileInputStream(propFile);
						prop.loadFromXML(fIS);
						fIS.close();

						String popCompositionKey = Simulation_ClusterModelTransmission.POP_PROP_INIT_PREFIX
								+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
						int[] gender_dist = (int[]) PropValUtils.propStrToObject(prop.getProperty(popCompositionKey),
								int[].class);

						int[] cuml_gender_dist = Arrays.copyOf(gender_dist, gender_dist.length);
						// Cumul. gender dist
						for (int i = 1; i < cuml_gender_dist.length; i++) {
							cuml_gender_dist[i] += cuml_gender_dist[i - 1];
						}

						int[][] total_no_incident_reported = new int[Population_Bridging.LENGTH_GENDER][Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE];

						// Key = String.format("%d,%d", g, s)
						HashMap<String, ArrayList<Integer>> inf_history_count_map = new HashMap<>();
						HashMap<String, ArrayList<Long>> inf_history_dur_map = new HashMap<>();
						HashMap<String, ArrayList<Integer>> inf_history_interval_map = new HashMap<>();

						for (File f : zipFiles) {
							SevenZFile resultZip = new SevenZFile(f);
							SevenZArchiveEntry ent;
							while ((ent = resultZip.getNextEntry()) != null) {

								// Assume no incident for all first
								for (int g = 0; g < total_no_incident_reported.length; g++) {
									for (int s = 0; s < total_no_incident_reported[g].length; s++) {
										total_no_incident_reported[g][s] += gender_dist[g];
									}
								}

								StringBuilder txt_entries = new StringBuilder();
								int count;
								while ((count = resultZip.read(buf, 0, BUFFER)) != -1) {
									txt_entries.append(new String(Arrays.copyOf(buf, count)));
								}

								Util_CSV_Table_Map.updateInfectionHistoryMap(inf_history_count_map, inf_history_dur_map,
										inf_history_interval_map, cuml_gender_dist, incl_range,
										total_no_incident_reported, ent.getName(), txt_entries.toString());

							}
							resultZip.close();
						}

						File summary_file = new File(baseDir,
								String.format(stat_filename_format, Arrays.toString(incl_range)));
						PrintWriter pWri = new PrintWriter(summary_file);
						// Incident count
						pWri.println("Gender, Site, # Incidence (bin-size of 1)");
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							for (int s = 0; s < Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
								pWri.print(g);
								pWri.print(',');
								pWri.print(s);

								ArrayList<Integer> history_count_map_ent = inf_history_count_map
										.get(String.format("%d,%d", g, s));
								int[] incident_count;

								if (history_count_map_ent != null) {
									Collections.sort(history_count_map_ent);
									Integer maxIncident = history_count_map_ent.get(history_count_map_ent.size() - 1);
									incident_count = new int[maxIncident + 1];
									incident_count[0] = total_no_incident_reported[g][s];
									for (Integer n : history_count_map_ent) {
										incident_count[n]++;
									}

								} else {
									incident_count = new int[] { total_no_incident_reported[g][s] };
								}

								for (int n : incident_count) {
									pWri.print(',');
									pWri.print(n);
								}

								pWri.println();
							}
						}

						Number[] keys, all_data;
						File histo_data_file;
						PrintWriter pWri_hist_data;
						double[] all_data_double;
						Percentile percent_data;
						File histBaseDir;

						// Infection duration
						histBaseDir = new File(baseDir, String.format("Durations_%s", Arrays.toString(incl_range)));
						histBaseDir.mkdirs();

						pWri.println("Duration of infection");
						pWri.println(
								"Gender, Site, Mean duration, Total infection, Total duration, Median, 25th Quartile, 75th Quartile");
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							for (int s = 0; s < Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
								pWri.print(g);
								pWri.print(',');
								pWri.print(s);

								// Total duration
								all_data = printDuration(pWri, inf_history_dur_map.get(String.format("%d,%d", g, s)));

								HashMap<Long, Long> hist_dur_col = new HashMap<>();
								for (int i = 0; i < all_data.length; i++) {
									Long ent = hist_dur_col.get(all_data[i].longValue());
									if (ent == null) {
										ent = 0l;

									}
									hist_dur_col.put(all_data[i].longValue(), ent + 1);
								}

								keys = hist_dur_col.keySet().toArray(new Long[hist_dur_col.size()]);
								Arrays.sort(keys);

								histo_data_file = new File(histBaseDir,
										String.format("Duration_all_gender_%d_site_%d.csv", g, s));
								pWri_hist_data = new PrintWriter(histo_data_file);
								for (Number key : keys) {
									pWri_hist_data.printf("%d,%d\n", key.longValue(), hist_dur_col.get(key));
								}
								pWri_hist_data.close();

							}
						}

						pWri.println("Duration of treated infection");
						pWri.println(
								"Gender, Site, Mean duration, Total infection, Total duration, Median, 25th Quartile, 75th Quartile");
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							for (int s = 0; s < Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
								pWri.print(g);
								pWri.print(',');
								pWri.print(s);
								printDuration(pWri, inf_history_dur_map.get(String.format("T_%d,%d", g, s)));
							}
						}

						histBaseDir = new File(baseDir,
								String.format("Infection_intervals_%s", Arrays.toString(incl_range)));
						histBaseDir.mkdirs();

						pWri.println("Interval between infections ");
						pWri.println("Gender, Site, Mean, Median, 25th Quartile, 75th Quartile");

						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							for (int s = 0; s < Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {

								pWri.print(g);
								pWri.print(',');
								pWri.print(s);

								// Infection interval
								ArrayList<Integer> history_interval_map_ent = inf_history_interval_map
										.get(String.format("%d,%d", g, s));

								all_data = history_interval_map_ent.toArray(new Integer[0]);
								all_data_double = new double[all_data.length];

								double interval_sum = 0;
								int interval_count = 0;
								for (int i = 0; i < all_data.length; i++) {
									all_data_double[i] = all_data[i].intValue();
									interval_count++;
									interval_sum += all_data_double[i];
								}

								pWri.print(',');
								pWri.print(interval_sum / interval_count);

								percent_data = new Percentile();
								percent_data.setData(all_data_double);

								pWri.print(',');
								pWri.print(percent_data.evaluate(50));

								pWri.print(',');
								pWri.print(percent_data.evaluate(25));

								pWri.print(',');
								pWri.print(percent_data.evaluate(75));

								pWri.println();

								HashMap<Integer, Integer> hist_interval_col = new HashMap<>();
								for (Integer interval : history_interval_map_ent) {
									Integer ent = hist_interval_col.get(interval);
									if (ent == null) {
										ent = 0;

									}
									hist_interval_col.put(interval, ent + 1);
								}

								keys = hist_interval_col.keySet().toArray(new Integer[hist_interval_col.size()]);
								Arrays.sort(keys);
								histo_data_file = new File(histBaseDir,
										String.format("Infection_interval_all_gender_%d_site_%d.csv", g, s));
								pWri_hist_data = new PrintWriter(histo_data_file);
								for (Number key : keys) {
									pWri_hist_data.printf("%d,%d\n", key.intValue(), hist_interval_col.get(key));
								}
								pWri_hist_data.close();

							}

						}

						pWri.close();

					} else {

						Util_CSV_Table_Map csvTableMapping = null;
						ArrayList<Util_CSV_Table_Map> csvTableExtra = new ArrayList<>();
						ArrayList<int[][]> csvTableExtra_colSel = new ArrayList<>();
						ArrayList<String> csvTableExtra_filename = new ArrayList<>();

						String replace_str = "(-{0,1}\\\\d+(?:_\\\\d+){0,1})";

						if (zipFileName.equals(Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON_ZIP
								.replaceAll("%d", replace_str))) {
							csvTableExtra.add(new Util_CSV_Table_Map("Time,Heterosexual,MSM,All"));
							csvTableExtra_colSel.add(
									new int[][] { new int[] { 1, 2 }, new int[] { 3, 4 }, new int[] { 1, 2, 3, 4 }, });
							csvTableExtra_filename.add("Summary_Prevalence_Person_BehavGrp_%s.csv");
						}

						if (zipFileName.equals(Simulation_ClusterModelTransmission.FILENAME_VACCINE_COVERAGE_PERSON_ZIP
								.replaceAll("%d", replace_str))) {
							csvTableExtra.add(new Util_CSV_Table_Map("Time,All_Active,All_Partial,All_Expired"));
							csvTableExtra_colSel.add(new int[][] { new int[] { 1, 5, 9, 13 },
									new int[] { 2, 6, 10, 14 }, new int[] { 3, 7, 11, 15 }, });
							csvTableExtra_filename.add("Summary_Vaccine_Person_BehavGrp_%s.csv");
						}

						if (zipFileName.equals(Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE_ZIP
								.replaceAll("%d", replace_str))) {
							csvTableExtra.add(new Util_CSV_Table_Map("Time,Site_1,Site_2,Site_3"));
							csvTableExtra_colSel.add(
									new int[][] { new int[] { 10, 14 }, new int[] { 11, 15 }, new int[] { 12, 16 }, });
							csvTableExtra_filename.add("Summary_Prevalence_MSM_Site_%s.csv");
						}

						if (zipFileName.equals(Simulation_ClusterModelTransmission.FILENAME_CUMUL_TREATMENT_PERSON_ZIP
								.replaceAll("%d", replace_str))) {
							Util_CSV_Table_Map table = new Util_CSV_Table_Map("Time,Female,Male");
							table.setCumulative(isCumul);
							csvTableExtra.add(table);
							csvTableExtra_colSel.add(new int[][] { new int[] { 1 }, new int[] { 2, 3, 4 } });
							csvTableExtra_filename.add("Summary_Treatment_BehavGrp_%s.csv");
						}

						if (zipFileName.equals(Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_PERSON_ZIP
								.replaceAll("%d", replace_str))) {
							Util_CSV_Table_Map table = new Util_CSV_Table_Map("Time,Female,Male");
							table.setCumulative(isCumul);
							csvTableExtra.add(table);
							csvTableExtra_colSel.add(new int[][] { new int[] { 1 }, new int[] { 2, 3, 4 } });
							csvTableExtra_filename.add("Summary_DX_BehavGrp_%s.csv");
						}

						if (zipFileName.equals(Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON_ZIP
								.replaceAll("%d", replace_str))) {
							Util_CSV_Table_Map table = new Util_CSV_Table_Map("Time,Female,All_Male,MSM");
							table.setCumulative(isCumul);
							csvTableExtra.add(table);
							csvTableExtra_colSel
									.add(new int[][] { new int[] { 1 }, new int[] { 2, 3, 4 }, new int[] { 3, 4 } });
							csvTableExtra_filename.add("Summary_Incident_BehavGrp_%s.csv");
						}

						for (File f : zipFiles) {
							HashMap<String, ArrayList<String[]>> file_ent = new HashMap<>();
							file_ent = Util_7Z_CSV_Entry_Extract_Callable.extractedLinesFrom7Zip(f, file_ent);

							// System.out.printf("From %s: %d simulations\n", f.getName(), file_ent.size());

							for (String zipEntName : file_ent.keySet()) {
								ArrayList<String[]> data = file_ent.get(zipEntName);
								if (csvTableMapping == null) {
									csvTableMapping = new Util_CSV_Table_Map(data.get(0));
									csvTableMapping.setCumulative(isCumul);
								}
								for (int r = 1; r < data.size(); r++) {								
									if (data.get(r).length > 0) {
										try {
											csvTableMapping.addRow(data.get(r));
											if (!csvTableExtra.isEmpty()) {
												String[] lineAtt = data.get(r);
												Integer time = Integer.parseInt(lineAtt[0]);
												for (int i = 0; i < csvTableExtra.size(); i++) {
													int[][] colSel = csvTableExtra_colSel.get(i);
													double[] col_selSum = new double[colSel.length];
													for (int c = 0; c < colSel.length; c++) {
														for (int cI : colSel[c]) {
															col_selSum[c] += Double.parseDouble(lineAtt[cI]);
														}
													}
													csvTableExtra.get(i).addRow(time, col_selSum);
												}
											}
										} catch (Exception ex) {
											System.err.printf("Error in adding row from %s\n", zipEntName);
										}
									}
								}
							}
						}

						if (csvTableMapping != null) {
							String summaryFileFormat = stat_filename_format;
							csvTableMapping.printSummaryFile(summaryFileFormat, baseDir);
						}
						if (!csvTableExtra.isEmpty()) {
							for (int i = 0; i < csvTableExtra.size(); i++) {
								csvTableExtra.get(i).printSummaryFile(csvTableExtra_filename.get(i), baseDir);
							}
						}

					}

				}

			}

			System.out.println("Output analysis completed.");
		} else

		{
			System.out.printf("Properties file %s NOT found. Exiting...\n", propFile.getAbsolutePath());
			System.exit(-1);
		}

	}

	private Number[] printDuration(PrintWriter pWri, ArrayList<Long> history_dur_map_ent) {
		Number[] all_data;
		double[] all_data_double;
		Percentile percent_data;
		pWri.print(',');
		pWri.print(1f * history_dur_map_ent.get(1) / history_dur_map_ent.get(0));

		pWri.print(',');
		pWri.print(history_dur_map_ent.get(0));

		pWri.print(',');
		pWri.print(history_dur_map_ent.get(1));

		all_data = history_dur_map_ent.subList(2, history_dur_map_ent.size()).toArray(new Long[0]);

		all_data_double = new double[all_data.length];
		for (int i = 0; i < all_data.length; i++) {
			all_data_double[i] = all_data[i].longValue();
		}

		percent_data = new Percentile();
		percent_data.setData(all_data_double);

		pWri.print(',');
		pWri.print(percent_data.evaluate(50));

		pWri.print(',');
		pWri.print(percent_data.evaluate(25));

		pWri.print(',');
		pWri.print(percent_data.evaluate(75));

		pWri.println();
		return all_data;
	}

	public static int[] calculateIncNoticationFromHistory(Integer pid, ArrayList<String[]> ent, int[] time_range,
			int incNotifRow) {

		int numInc = 0;
		int numNotif = 0;
		int[] time_span = new int[] { -1, -1 };

		if (incNotifRow < ent.size()) {
			String[] inc_notif_row = ent.get(incNotifRow);
			for (int s = 2; s < inc_notif_row.length; s++) { // Offset: pid, site...
				String time_str = inc_notif_row[s];
				Integer nextTime = Integer.parseInt(time_str);
				time_range[0] = Math.min(time_range[0], Math.abs(nextTime));
				time_range[1] = Math.max(time_range[1], Math.abs(nextTime));
				if (nextTime < 0) {
					numNotif++;
				} else {
					numInc++;
				}
			}
		} else {
			// Backward compatibly for result prior to 20240206

			int[] timing_pt = new int[ent.size()];
			Arrays.fill(timing_pt, 2);
			int inf_stat = 0;
			int nextTimeRow;

			while ((nextTimeRow = getNextTimeRowFromHistory(ent, timing_pt)) >= 0) {
				int nextTime = Integer.parseInt(ent.get(nextTimeRow)[timing_pt[nextTimeRow] - 1]);

				time_range[0] = Math.min(time_range[0], Math.abs(nextTime));
				time_range[1] = Math.max(time_range[1], Math.abs(nextTime));

				if (time_span[0] < 0) {
					time_span[0] = nextTime;
				}
				if (time_span[1] < nextTime) {
					time_span[1] = nextTime;
				}

				if (nextTime < 0) {
					numNotif++;
					inf_stat = 0; // Treatment for all
				} else {
					if ((inf_stat & 1 << nextTimeRow) == 0) { // New infection at site

						if (inf_stat == 0) { // new infection at any site
							numInc++;
						}

						inf_stat |= 1 << nextTimeRow;
					} else {
						if ((inf_stat & 1 << nextTimeRow) != 0) {
							inf_stat -= 1 << nextTimeRow;
						}
					}
				}
			}
		}

		int[] mappingEnt = new int[rISK_GRP_MAP_LENGTH];
		mappingEnt[rISK_GRP_MAP_INDEX_PID] = pid;
		mappingEnt[rISK_GRP_MAP_INDEX_NUM_INC] = numInc;
		mappingEnt[rISK_GRP_MAP_INDEX_NUM_NOTIF] = numNotif;
		mappingEnt[rISK_GRP_MAP_INDEX_NUM_TIME_SPAN] = time_span[1] - time_span[0];
		return mappingEnt;
	}

	private static int getNextTimeRowFromHistory(ArrayList<String[]> ent, int[] timing_pt) {
		int nextTimeRow = -1;
		int nextTime = Integer.MAX_VALUE;
		int nextTime_ent = Integer.MAX_VALUE;
		for (int i = 0; i < timing_pt.length; i++) {
			String[] ent_s = ent.get(i);
			if (timing_pt[i] < ent_s.length && Math.abs(Integer.parseInt(ent_s[timing_pt[i]])) < nextTime) {
				nextTime = Math.abs(Integer.parseInt(ent_s[timing_pt[i]]));
				nextTime_ent = Integer.parseInt(ent_s[timing_pt[i]]);
				nextTimeRow = i;

			}
		}
		if (nextTimeRow != -1) {
			for (int i = 0; i < timing_pt.length; i++) {
				String[] ent_s = ent.get(i);
				if (timing_pt[i] < ent_s.length && nextTime_ent == Integer.parseInt(ent_s[timing_pt[i]])) {
					timing_pt[i]++;
				}

			}
		}

		return nextTimeRow;
	}

	public static void cleanUpOutputDir(File cleanupDir) {
		System.out.printf("Cleaning up outputs in %s\n", cleanupDir.getAbsolutePath());
		Pattern outputZipsPattern = Pattern.compile("(.+_)(-{0,1}(?!0)\\d+)\\.csv\\.7z");
		HashMap<Long, ArrayList<Long>> sim_listing = new HashMap<>();

		File[] archiveFiles = cleanupDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return outputZipsPattern.matcher(pathname.getName()).matches();
			}
		});

		for (File archiveFile : archiveFiles) {
			Matcher mArchive = outputZipsPattern.matcher(archiveFile.getName());
			mArchive.matches();

			String filePrefix = mArchive.group(1);
			long cMapSeed = Long.parseLong(mArchive.group(2));
			Pattern ent_patten = Pattern.compile(String.format("%s%d_(-{0,1}(?!0)\\d+)\\.csv", filePrefix, cMapSeed));
			Pattern remove_pattern = Pattern
					.compile(String.format("%s%d\\.csv\\.7z_(-{0,1}(?!0)\\d+)\\.7z", filePrefix, cMapSeed));

			ArrayList<Long> simSeedArr = sim_listing.get(cMapSeed);
			if (simSeedArr == null) {
				simSeedArr = new ArrayList<>();
				sim_listing.put(cMapSeed, simSeedArr);
			}

			try {
				SevenZFile archive7Z = new SevenZFile(archiveFile);
				SevenZArchiveEntry ent;

				while ((ent = archive7Z.getNextEntry()) != null) {
					Matcher mEnt = ent_patten.matcher(ent.getName());
					mEnt.matches();
					long simSeed = Long.parseLong(mEnt.group(1));
					int res = Collections.binarySearch(simSeedArr, simSeed);
					if (res < 0) {
						simSeedArr.add(~res, simSeed);
					}
				}
				archive7Z.close();

				// Remove backup 7z files
				File[] toBeRemoved = cleanupDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return remove_pattern.matcher(pathname.getName()).matches();
					}
				});

				for (File removed : toBeRemoved) {
					FileUtils.delete(removed);
				}

				// Unzipped CSV files
				File[] extra_CSVs = cleanupDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return ent_patten.matcher(pathname.getName()).matches();
					}
				});

				if (extra_CSVs.length > 0) {
					if (archiveFile.length() > 100 * 1024) {
						System.out.printf(
								"%d unzipped CSV found by the archive file %s is too large (%d bytes) to rezip automatically.\n",
								extra_CSVs.length, archiveFile.getName(), archiveFile.length());

					} else {

						SevenZArchiveEntry entry, inputEnt;
						FileInputStream fIn;

						File preZip = new File(cleanupDir,
								archiveFile.getName() + "_" + Long.toString(System.currentTimeMillis()) + ".7z");
						Files.copy(archiveFile.toPath(), preZip.toPath(), StandardCopyOption.COPY_ATTRIBUTES);

						SevenZOutputFile outputZip = new SevenZOutputFile(archiveFile);
						SevenZFile inputZip = new SevenZFile(preZip);

						final int BUFFER = 2048000;
						byte[] buf = new byte[BUFFER];
						while ((inputEnt = inputZip.getNextEntry()) != null) {
							outputZip.putArchiveEntry(inputEnt);
							int count;
							while ((count = inputZip.read(buf, 0, BUFFER)) != -1) {
								outputZip.write(Arrays.copyOf(buf, count));
							}
							outputZip.closeArchiveEntry();
						}
						inputZip.close();

						for (int fI = 0; fI < extra_CSVs.length; fI++) {
							entry = outputZip.createArchiveEntry(extra_CSVs[fI], extra_CSVs[fI].getName());
							outputZip.putArchiveEntry(entry);
							fIn = new FileInputStream(extra_CSVs[fI]);
							outputZip.write(fIn);
							outputZip.closeArchiveEntry();
							fIn.close();
						}
						outputZip.close();

						preZip.delete();
						for (File f : extra_CSVs) {
							f.delete();
						}

					}
				}

			} catch (IOException ex) {
				System.err.printf("File clean up for %s failed. File not removed.\n", archiveFile.getAbsolutePath());
				ex.printStackTrace(System.err);

			}
		}

		// Adding SeedIndexCases file if needed
		for (Long cMapSeed : sim_listing.keySet()) {

			System.out.printf("\tCMap_Seed = %d, # sim = %d\n", cMapSeed, sim_listing.get(cMapSeed).size());

			String seedFileName = String.format(Simulation_ClusterModelTransmission.FILENAME_INDEX_CASE_LIST_ZIP,
					cMapSeed);
			File seedFile = new File(cleanupDir, seedFileName);
			try {
				if (!seedFile.exists()) {
					ArrayList<Long> simSeedArr = sim_listing.get(cMapSeed);
					SevenZOutputFile outputZip = new SevenZOutputFile(seedFile);
					SevenZArchiveEntry entry;
					FileInputStream fIn;

					for (Long simSeed : simSeedArr) {
						String genSeedFileName = String.format(
								Simulation_ClusterModelTransmission.FILENAME_INDEX_CASE_LIST, cMapSeed, simSeed);
						File dummmyFile = new File(cleanupDir, genSeedFileName);
						PrintWriter pWri = new PrintWriter(dummmyFile);
						pWri.printf("Dummy entry created at %d", System.currentTimeMillis());
						pWri.close();
						entry = outputZip.createArchiveEntry(dummmyFile, genSeedFileName);
						outputZip.putArchiveEntry(entry);
						fIn = new FileInputStream(dummmyFile);
						outputZip.write(fIn);
						outputZip.closeArchiveEntry();
						fIn.close();
						FileUtils.delete(dummmyFile);
					}
					outputZip.close();

				} else {
					// Remove backup seed archive
					Pattern remove_pattern_seed = Pattern
							.compile(String.format("%s_(-{0,1}(?!0)\\d+)\\.7z", seedFileName));
					File[] toBeRemovedSeed = cleanupDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return remove_pattern_seed.matcher(pathname.getName()).matches();
						}
					});
					for (File removed : toBeRemovedSeed) {
						FileUtils.delete(removed);
					}
				}
			} catch (IOException e) {
				System.err.printf("File operation for %s failed.\n", seedFile.getAbsolutePath());
				e.printStackTrace(System.err);
			}

		}
	}

	@SuppressWarnings("rawtypes")
	public static void extractTrendResults(File basedir, File fittingTarDir) throws IOException {

		final Pattern sEED_FILE_PATTERN = Pattern
				.compile(Simulation_ClusterModelTransmission.FILENAME_INDEX_CASE_LIST_ZIP.replaceAll("%d",
						"(-{0,1}(?!0)\\\\d+)"));

		final int[] tIME_POINTS = new int[] { 1095, 1460, 1825, 2190, 2555, 2920, 3285, 3650, 4015, 4380, 4745, 5110,
				5475, 5840, 6205, 6570, 6935, 7300 };

		final int tIME_POINTS_FITTING_PT = 7;
		final int nUM_GROUP = 4;
		final int nUM_SITE = 4;

		File[] fittingTargets = fittingTarDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.getName().endsWith(".csv");
			}
		});

		PrintWriter pri_seed_key = new PrintWriter(new File(basedir, "Trent_ExtractMapping.csv"));
		PrintWriter[] pri_trend_output = new PrintWriter[fittingTargets.length];

		// Load trend CSV
		// Key: Path,type,tar_grp,weight
		HashMap<String, double[][]> target_trend_collection = OptTrendFittingFunction.loadTrendCSV(fittingTargets);

		target_trend_collection.remove(OptTrendFittingFunction.OPT_TREND_CSV_RANGE);
		String[] trend_target_key = target_trend_collection.keySet()
				.toArray(new String[target_trend_collection.size()]);

		String[][] trend_target_key_split = new String[trend_target_key.length][];

		// Seed list

		File[] seedFile = basedir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return sEED_FILE_PATTERN.matcher(pathname.getName()).matches();
			}
		});
		String[] seedList = new String[seedFile.length];
		for (int s = 0; s < seedList.length; s++) {
			Matcher m = sEED_FILE_PATTERN.matcher(seedFile[s].getName());
			m.matches();
			seedList[s] = m.group(1);
		}

		ArrayList<String> csvFilePatternArr = new ArrayList<>();
		ArrayList<int[]> entryPairArr = new ArrayList<>();

		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] trend_incl = new ArrayList[pri_trend_output.length];

		for (int i = 0; i < pri_trend_output.length; i++) {
			String[] trend_keys = trend_target_key[i].split(",");
			trend_target_key_split[i] = trend_keys;

			pri_trend_output[i] = new PrintWriter(new File(basedir, String.format("Trent_Extract_%s",
					new File(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_PATH]).getName())));

			String type_key = trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_TYPE];
			int incl_grp = Integer.parseInt(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_TARGET_GRP]);
			int site = -1;

			Matcher m = OptTrendFittingFunction.OPT_TREND_TYPE_FORMAT_BY_SITE
					.matcher(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_TYPE]);
			if (m.find()) {
				site = Integer.parseInt(m.group(2));
				type_key = m.group(1);
			}

			String[] csv_file_pattern_str_arr = new String[0];

			if (OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_NUMINF.equals(type_key)
					|| OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_INCID.equals(type_key)) {
				String csv_file_pattern_str;
				if (site < 0) {
					if (OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_NUMINF.equals(type_key)) {
						csv_file_pattern_str = Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON;
					} else {
						csv_file_pattern_str = Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON;
					}
				} else {
					if (OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_NUMINF.equals(type_key)) {
						csv_file_pattern_str = Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE;
					} else {
						csv_file_pattern_str = Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_SITE;
					}
				}

				csv_file_pattern_str_arr = new String[] { csv_file_pattern_str };
				int index = csvFilePatternArr.indexOf(csv_file_pattern_str);
				if (index < 0) {
					index = csvFilePatternArr.size();
					csvFilePatternArr.add(csv_file_pattern_str);

				}
			} else if (OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_DX.equals(type_key)) {
				csv_file_pattern_str_arr = new String[] {
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_PERSON,
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON };
			}

			for (String csv_file_pattern_str : csv_file_pattern_str_arr) {
				int index = csvFilePatternArr.indexOf(csv_file_pattern_str);
				if (index < 0) {
					index = csvFilePatternArr.size();
					csvFilePatternArr.add(csv_file_pattern_str);
				}
				for (int g = 0; g < nUM_GROUP; g++) {
					if ((incl_grp & 1 << g) != 0) {
						int row_index = entryPairArr.size();
						if (trend_incl[i] == null) {
							trend_incl[i] = new ArrayList<>();
						}
						trend_incl[i].add(row_index);
						if (site < 0) {
							// Time,Gender_0,Gender_1,Gender_2,Gender_3
							// Time,Total_Positive_DX_Gender_0,...
							entryPairArr.add(new int[] { index, g + 1 });
						} else {
							// Time,Gender_0_Site_0,Gender_0_Site_1,Gender_0_Site_2,Gender_0_Site_3,...Gender_1_Site_0
							entryPairArr.add(new int[] { index, 1 + g * nUM_SITE + site });
						}
					}
				}
			}

			pri_trend_output[i].println("Time-value pairing");

		}

		// Extract values
		pri_seed_key.println("Residue,CMap_Seed,Sim_Seed");
		String[] csv_Filenames_Pattern = csvFilePatternArr.toArray(new String[csvFilePatternArr.size()]);
		int[][] entry_Pairs = entryPairArr.toArray(new int[entryPairArr.size()][]);

		int numProcess = Runtime.getRuntime().availableProcessors();
		ExecutorService exec = Executors.newFixedThreadPool(numProcess);
		HashMap<String, Future<Map<String, long[]>>> res_collection = new HashMap<>();

		for (String mapSeed : seedList) {
			Util_7Z_CSV_Entry_Extract_Callable callable = new Util_7Z_CSV_Entry_Extract_Callable(basedir, mapSeed,
					csv_Filenames_Pattern, entry_Pairs, tIME_POINTS);
			res_collection.put(mapSeed, exec.submit(callable));
		}

		exec.shutdown();
		try {
			if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
				System.err.println("Thread time-out!");
			}

			ArrayList<Comparable[]> resKeyPrintMapArr = new ArrayList<>();
			HashMap<String, long[][]> keyPrintMap = new HashMap<>();

			for (String mapSeed : seedList) {
				Future<Map<String, long[]>> resultFuture = res_collection.get(mapSeed);
				Map<String, long[]> result = resultFuture.get();
				String[] keySet = result.keySet().toArray(new String[result.size()]);

				for (String key : keySet) {
					String[] key_s = key.split("_");
					String key_print = String.format("%s,%s", key_s[0], key_s[1]);
					long[][] entry = keyPrintMap.get(key_print);
					if (entry == null) {
						entry = new long[tIME_POINTS.length][];
						keyPrintMap.put(key_print, entry);
					}

					int time = Integer.parseInt(key_s[2]);
					long[] ent = result.get(key);

					int index = Arrays.binarySearch(tIME_POINTS, time);
					if (index >= 0) {
						entry[index] = ent;
					}
				}
			}

			// Calculate residue
			for (String keyPrint : keyPrintMap.keySet()) {
				long[][] val_all = keyPrintMap.get(keyPrint);

				Comparable[] resKeyPrintMap = new Comparable[2]; // Residue, keyPrint

				UnivariateFunction[] interpolation = new PolynomialSplineFunction[fittingTargets.length];
				UnivariateInterpolator interpolator = new LinearInterpolator();

				for (int i = 0; i < pri_trend_output.length; i++) {
					double[] t_val = new double[tIME_POINTS_FITTING_PT];
					double[] y_val = new double[tIME_POINTS_FITTING_PT];
					for (int t = 0; t < t_val.length; t++) {
						t_val[t] = tIME_POINTS[t];
						for (int c : trend_incl[i]) {
							y_val[t] += val_all[t][c];
						}
					}
					interpolation[i] = interpolator.interpolate(t_val, y_val);
				}
				double residue = 0;
				for (int i = 0; i < pri_trend_output.length; i++) {
					double[][] target_val = target_trend_collection.get(trend_target_key[i]);
					double offset = 0;

					if (trend_target_key_split[i][OptTrendFittingFunction.OPT_TREND_MAP_KEY_TYPE].startsWith("Cumul")) {
						offset = interpolation[i].value(tIME_POINTS[0]);
					}

					double weight = Double
							.parseDouble(trend_target_key_split[i][OptTrendFittingFunction.OPT_TREND_MAP_KEY_WEIGHT]);

					for (int t = 0; t < target_val[0].length; t++) {
						double tar2modelTime = target_val[0][t] + tIME_POINTS[0];
						double target_y = target_val[1][t];
						double model_y = interpolation[i].value(tar2modelTime);
						residue += weight * Math.pow((model_y - offset) - target_y, 2);
					}
				}

				resKeyPrintMap[0] = residue;
				resKeyPrintMap[1] = keyPrint;
				resKeyPrintMapArr.add(resKeyPrintMap);

			}

			Comparable[][] resKeyPrintMap = resKeyPrintMapArr.toArray(new Comparable[resKeyPrintMapArr.size()][]);

			Arrays.sort(resKeyPrintMap, (o1, o2) -> {
				int res = Double.compare((Double) o1[0], (Double) o2[0]);

				if (res == 0) {
					res = ((String) o1[1]).compareTo((String) o2[1]);
				}
				return res;
			});

			for (Comparable[] resKeyPrint : resKeyPrintMap) {
				Double res = (Double) resKeyPrint[0];
				String keyPrint = (String) resKeyPrint[1];
				pri_seed_key.printf("%s,%s\n", res.toString(), keyPrint);

				long[][] val_all = keyPrintMap.get(keyPrint);

				for (int i = 0; i < pri_trend_output.length; i++) {
					long[] val = new long[tIME_POINTS.length];
					for (int t = 0; t < val.length; t++) {
						for (int c : trend_incl[i]) {
							val[t] += val_all[t][c];
						}
						if (t > 0) {
							pri_trend_output[i].print(',');

						}
						pri_trend_output[i].print(tIME_POINTS[t]);
					}
					for (int t = 0; t < val.length; t++) {
						pri_trend_output[i].print(',');
						pri_trend_output[i].print(val[t]);
					}
					pri_trend_output[i].println();
				}

			}

		} catch (InterruptedException | ExecutionException e) {
			e.printStackTrace(System.err);
		}

		pri_seed_key.close();
		for (PrintWriter pWri : pri_trend_output) {
			pWri.close();
		}

	}

}
