package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Properties;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;

import population.Population_Bridging;
import sim.Runnable_ClusterModel_Transmission;
import sim.SimulationInterface;
import sim.Simulation_ClusterModelTransmission;

public class Util_Analyse_ClusterModel_Transmission_Output {

	private File baseDir;

	public static final String[] ZIP_FILES_LIST = new String[] {
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_TREATMENT_PERSON_ZIP.replaceAll("%d",
					"(-{0,1}(?!0)\\\\d+)"),
			Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON_ZIP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"),
			Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE_ZIP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"),
			Simulation_ClusterModelTransmission.FILENAME_INFECTION_HISTORY_ZIP.replaceAll("%d",
					"(-{0,1}(?!0)\\\\d+)"),
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_PERSON_ZIP.replaceAll("%d",
					"(-{0,1}(?!0)\\\\d+)"),
			Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON_ZIP.replaceAll("%d",
					"(-{0,1}(?!0)\\\\d+)"),
			Simulation_ClusterModelTransmission.FILENAME_VACCINE_COVERAGE_PERSON_ZIP.replaceAll("%d",
					"(-{0,1}(?!0)\\\\d+)"),
			};

	public static final String[] STAT_FILEFORMAT = new String[] { "Summary_Treatment_Person_%s.csv",
			"Summary_Prevalence_Person_%s.csv", "Summary_Prevalence_Site_%s.csv", "Summary_Infection_History.csv",  
			"Summary_DX_Person_%s.csv", "Summary_DX_Sought_Person_%s.csv", "Summary_Vaccine_Person_%s.csv" };

	public static final boolean[] CUMUL_DATA = new boolean[] { true, false, false, false, true, true, false };

	public Util_Analyse_ClusterModel_Transmission_Output() {

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

				final Pattern pattern_zip = Pattern.compile(zipFileName);

				File[] zipFiles = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pattern_zip.matcher(pathname.getName()).matches();
					}
				});

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
						int[] incl_range = new int[] { 2920, 4745 }; // 5 years
						
						int[] cuml_gender_dist = Arrays.copyOf(gender_dist, gender_dist.length);
						// Cumul. gender dist
						for (int i = 1; i < cuml_gender_dist.length; i++) {
							cuml_gender_dist[i] += cuml_gender_dist[i - 1];
						}
						
						int[][] total_no_incident_reported = new int[Population_Bridging.LENGTH_GENDER][Runnable_ClusterModel_Transmission.LENGTH_SITE];

						HashMap<String, ArrayList<Integer>> inf_history_map = new HashMap<>();
						
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


								Util_CSV_Table_Map.updateInfectionHistoryMap(inf_history_map, 
										cuml_gender_dist, incl_range, total_no_incident_reported,
										 ent.getName(), txt_entries.toString());

							}
							resultZip.close();
						}

						File summaryFile = new File(baseDir, stat_filename_format);
						PrintWriter pWri = new PrintWriter(summaryFile);
						pWri.println("Gender, Site, # Incident");
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							for (int s = 0; s < Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
								pWri.print(g);
								pWri.print(',');
								pWri.print(s);

								ArrayList<Integer> history_map_ent = inf_history_map.get(String.format("%d,%d", g, s));
								int[] incident_count;

								if (history_map_ent != null) {
									Collections.sort(history_map_ent);
									Integer maxIncident = history_map_ent.get(history_map_ent.size() - 1);
									incident_count = new int[maxIncident + 1];
									incident_count[0] = total_no_incident_reported[g][s];
									for (Integer n : history_map_ent) {
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
						pWri.close();

					} else {

						Util_CSV_Table_Map csvTableMapping = null;
						for (File f : zipFiles) {
							SevenZFile resultZip = new SevenZFile(f);
							SevenZArchiveEntry ent;
							while ((ent = resultZip.getNextEntry()) != null) {
								StringBuilder txt_entries = new StringBuilder();
								int count;
								while ((count = resultZip.read(buf, 0, BUFFER)) != -1) {
									txt_entries.append(new String(Arrays.copyOf(buf, count)));
								}

								BufferedReader lines = new BufferedReader(new StringReader(txt_entries.toString()));
								String line = lines.readLine(); // First line
								if (csvTableMapping == null) {
									csvTableMapping = new Util_CSV_Table_Map(line);
									csvTableMapping.setCumulative(isCumul);
								}

								while ((line = lines.readLine()) != null) {
									if (line.length() > 0) {
										try {
											csvTableMapping.addRow(line);
										} catch (Exception ex) {
											System.err.printf("Error in adding row from %s (%s)\n", ent.getName(),
													line);
										}
									}
								}
							}
							resultZip.close();
						}

						if (csvTableMapping != null) {
							String[] headers = csvTableMapping.getHeader();
							for (int s = 1; s < headers.length; s++) {
								String summary = csvTableMapping.displayStat(s);
								File summaryFile = new File(baseDir, String.format(stat_filename_format, headers[s]));
								PrintWriter pWri = new PrintWriter(summaryFile);
								pWri.println(summary);
								pWri.close();
							}

						}
					}

				}

			}
			
			System.out.println("Output analysis completed.");
		} else {
			System.out.printf("Properties file %s NOT found. Exiting...\n", propFile.getAbsolutePath());
			System.exit(-1);
		}

	}

	/*
	public static void main(String[] args) throws IOException {

		Util_Analyse_ClusterModel_Transmission_Output analysis = new Util_Analyse_ClusterModel_Transmission_Output();

		if (args.length == 0) {
			System.out.printf("Usage: java %s PROP_FILE_DIRECTORY\n", analysis.getClass().toString());
			System.exit(0);
		} else {
			analysis.setBaseDir(new File(args[0]));
			analysis.analyse_outputs();
		}

	}
	*/

}
