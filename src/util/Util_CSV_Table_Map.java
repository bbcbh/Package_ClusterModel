package util;

import java.io.BufferedReader;
import java.io.IOException;
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
		for (int s =0; s < values.length; s++) {
			String key = String.format("%d,%d", time, s+1);
			double value = values[s];
			ArrayList<Double> map_ent = this.get(key);
			if (map_ent == null) {
				map_ent = new ArrayList<>();
				this.put(key, map_ent);
			}
			map_ent.add(value);
		}
		
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

		double[] pre_row = null;
		for (Integer time : time_pt) {
			String map_key = String.format("%d,%d", time, colNum);
			ArrayList<Double> map_values = this.get(map_key);
			if (map_values != null) {
				double[] val = new double[map_values.size()];
				int p = 0;				
				for (Double d : map_values) {
					val[p] = d.doubleValue();					
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
					for(double singleVal : val) {
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

		return str.toString();
	}

	public static void updateInfectionHistoryMap(HashMap<String, ArrayList<Integer>> infection_history_map,
			int[] cumulative_gender_distribution, int[] incl_time_range, int[][] total_no_incident_reported_so_far,
			String src_file_name, String src_file_lines) throws IOException {
		BufferedReader lines = new BufferedReader(new StringReader(src_file_lines));
		String line;
		while ((line = lines.readLine()) != null) {
			if (line.length() > 0) {
				try {
					String[] entries = line.split(",");
					int id = Integer.parseInt(entries[0]);
					int gender = Runnable_ClusterModel_Transmission.getGenderType(id,
							cumulative_gender_distribution);
					int site = Integer.parseInt(entries[1]);
	
					int numInfections = 0;
					for (int i = 2; i < entries.length; i += 2) {
						int inf_start = Integer.parseInt(entries[i]);
						if (inf_start >= incl_time_range[0] && inf_start < incl_time_range[1]) {
							numInfections++;
						}
					}
	
					ArrayList<Integer> history_map_ent = infection_history_map
							.get(String.format("%d,%d", gender, site));
					if (history_map_ent == null) {
						history_map_ent = new ArrayList<>();
						infection_history_map.put(String.format("%d,%d", gender, site),
								history_map_ent);
					}
					history_map_ent.add(numInfections);
	
					total_no_incident_reported_so_far[gender][site]--;
	
				} catch (Exception ex) {
					System.err.printf("Error in adding line from %s (%s)\n", src_file_name,
							line);
				}
			}
		}
	}

}