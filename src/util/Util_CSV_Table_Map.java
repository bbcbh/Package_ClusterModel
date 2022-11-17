package util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

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
				double sum = 0;
				for (Double d : map_values) {
					val[p] = d.doubleValue();
					sum += val[p];
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

}