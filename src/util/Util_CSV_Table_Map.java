package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import sim.Runnable_ClusterModel_Transmission;

public class Util_CSV_Table_Map extends HashMap<String, ArrayList<Double>> {
	private static final long serialVersionUID = 3940051551837596423L;

	private String[] header;
	private ArrayList<Integer> time_pt;
	boolean isCumulative = false;

	public Util_CSV_Table_Map(String[] headerRow) {
		this(concatStr(headerRow));
	}

	public Util_CSV_Table_Map(String headerRow) {
		header = headerRow.split(",");
		time_pt = new ArrayList<>();
	}

	public String[] getHeader() {
		return header;
	}

	public ArrayList<Integer> getTime_pt() {
		return time_pt;
	}

	public void setCumulative(boolean isCumulative) {
		this.isCumulative = isCumulative;
	}

	public void addRow(Integer time, double[] values) {
		int time_index = Collections.binarySearch(time_pt, time);
		if (time_index < 0) {
			time_pt.add(~time_index, time);
		}
		for (int s = 0; s < values.length; s++) {
			String key = String.format("%d,%d", time, s + 1);
			double value = values[s];
			ArrayList<Double> map_ent = this.get(key);
			if (map_ent == null) {
				map_ent = new ArrayList<>();
				this.put(key, map_ent);
			}
			map_ent.add(value);
		}

	}

	public void addRow(String[] csv_line_sp) {
		addRow(concatStr(csv_line_sp));
	}

	public void addRow(String csv_line) {
		String[] entry = csv_line.split(",");
		int time = Integer.parseInt(entry[0]);
		int time_index = Collections.binarySearch(time_pt, time);
		if (time_index < 0) {
			time_pt.add(~time_index, time);
		}
		for (int s = 1; s < entry.length; s++) {
			String key = String.format("%d,%d", time, s);
			double value = Double.parseDouble(entry[s]);
			ArrayList<Double> map_ent = this.get(key);
			if (map_ent == null) {
				map_ent = new ArrayList<>();
				this.put(key, map_ent);
			}
			map_ent.add(value);
		}
	}

	public String displayStat(int colNum) {
		StringBuilder str = new StringBuilder();
		Percentile per = new Percentile();
		str.append("Time, N, Mean, Median, Q1, Q3\n");

		boolean hasData = false;

		double[] pre_row = null;
		for (Integer time : time_pt) {
			String map_key = String.format("%d,%d", time, colNum);
			ArrayList<Double> map_values = this.get(map_key);
			if (map_values != null) {
				double[] val = new double[map_values.size()];
				int p = 0;
				for (Double d : map_values) {
					val[p] = d.doubleValue();
					hasData |= val[p] != 0;					
					p++;
				}																	
				if (!isCumulative) {
					per.setData(val);
				} else {
					if (pre_row == null) {
						pre_row = val;
						val = null;
					} else {
						for (int i = 0; i < val.length; i++) {
							double currentVal = val[i];
							val[i] = val[i] - pre_row[i];
							pre_row[i] = currentVal;
						}																						
						per.setData(val);
					}
				}
				if (val != null) {
					double sum = 0;
					for (double singleVal : val) {
						sum += singleVal;
					}
					str.append(time);
					str.append(',');
					str.append(map_values.size());
					str.append(',');
					str.append(sum / map_values.size());
					str.append(',');
					str.append(per.evaluate(50));
					str.append(',');
					str.append(per.evaluate(25));
					str.append(',');
					str.append(per.evaluate(75));
					str.append('\n');
				}

			}
		}

		if (hasData) {
			return str.toString();
		} else {
			return null;
		}
	}

	public void printSummaryFile(String summaryFileFormat, File baseDir) throws FileNotFoundException {
		String dirName = summaryFileFormat;
		int subIndex = summaryFileFormat.indexOf("_%s");
		if (subIndex > 0) {
			dirName = summaryFileFormat.substring(0, subIndex);
		}
		File resultsDir = new File(baseDir, dirName);
		resultsDir.mkdirs();

		String[] headers = getHeader();
		for (int s = 1; s < headers.length; s++) {
			String summary = displayStat(s);
			if (summary != null) {
				File summaryFile = new File(resultsDir, String.format(summaryFileFormat, headers[s]));
				PrintWriter pWri = new PrintWriter(summaryFile);
				pWri.println(summary);
				pWri.close();
			}
		}
	}

	private static String concatStr(String[] splitedStr) {
		StringBuilder res = new StringBuilder();
		for (String hE : splitedStr) {
			if (res.length() > 0) {
				res.append(',');
			}
			res.append(hE);
		}
		return res.toString();
	}

	public static void updateInfectionHistoryMap(HashMap<String, ArrayList<Integer>> infection_history_count_map,
			HashMap<String, ArrayList<Long>> infection_history_duration_map,
			HashMap<String, ArrayList<Integer>> inf_history_interval_map, int[] cumulative_gender_distribution,
			int[] incl_time_range, int[][] total_no_incident_reported_so_far, String src_file_name,
			String src_file_lines) throws IOException {
		BufferedReader lines = new BufferedReader(new StringReader(src_file_lines));
		String line;
		while ((line = lines.readLine()) != null) {
			if (line.length() > 0) {
				try {
					String[] entries = line.split(",");
					int id = Integer.parseInt(entries[0]);
					int gender = Runnable_ClusterModel_Transmission.getPersonGrp(id, cumulative_gender_distribution);
					int site = Integer.parseInt(entries[1]);

					String key = String.format("%d,%d", gender, site);
					String key_treatment = String.format("T_%d,%d", gender, site);

					ArrayList<Integer> history_map_count_ent = infection_history_count_map.get(key);
					if (history_map_count_ent == null) {
						history_map_count_ent = new ArrayList<>();
						infection_history_count_map.put(key, history_map_count_ent);
					}

					ArrayList<Long> history_map_dur_map_ent = infection_history_duration_map.get(key);
					if (history_map_dur_map_ent == null) {
						history_map_dur_map_ent = new ArrayList<>();
						infection_history_duration_map.put(key, history_map_dur_map_ent);

						history_map_dur_map_ent.add(0l);
						history_map_dur_map_ent.add(0l);
					}

					ArrayList<Long> history_map_dur_map_treatment_ent = infection_history_duration_map
							.get(key_treatment);
					if (history_map_dur_map_treatment_ent == null) {
						history_map_dur_map_treatment_ent = new ArrayList<>();
						infection_history_duration_map.put(key_treatment, history_map_dur_map_treatment_ent);
						history_map_dur_map_treatment_ent.add(0l);
						history_map_dur_map_treatment_ent.add(0l);
					}

					ArrayList<Integer> history_map_interval_map_ent = inf_history_interval_map.get(key);
					if (history_map_interval_map_ent == null) {
						history_map_interval_map_ent = new ArrayList<>();
						inf_history_interval_map.put(key, history_map_interval_map_ent);
					}

					int numInfections = 0;
					int lastInfectionTimeStamp = -1;

					for (int i = 2; i < entries.length; i += 2) {
						int inf_start = Integer.parseInt(entries[i]);
						if (inf_start >= incl_time_range[0] && inf_start < incl_time_range[1]) {
							numInfections++;
							if (i + 1 < entries.length) {
								int inf_end = Math.abs(Integer.parseInt(entries[i + 1]));
								int dur = inf_end - inf_start;
								history_map_dur_map_ent.set(0, history_map_dur_map_ent.get(0) + 1);
								history_map_dur_map_ent.set(1, history_map_dur_map_ent.get(1) + dur);
								history_map_dur_map_ent.add((long) dur);

								if (Integer.parseInt(entries[i + 1]) < 0) {
									history_map_dur_map_treatment_ent.set(0,
											history_map_dur_map_treatment_ent.get(0) + 1);
									history_map_dur_map_treatment_ent.set(1,
											history_map_dur_map_treatment_ent.get(1) + dur);
									history_map_dur_map_treatment_ent.add((long) dur);
								}

							}
							if (lastInfectionTimeStamp > 0) {
								history_map_interval_map_ent.add(inf_start - lastInfectionTimeStamp);
							}
							lastInfectionTimeStamp = inf_start;
						}
					}
					history_map_count_ent.add(numInfections);
					total_no_incident_reported_so_far[gender][site]--;

				} catch (Exception ex) {
					System.err.printf("Error in adding line from %s (%s)\n", src_file_name, line);
				}
			}
		}
	}

}