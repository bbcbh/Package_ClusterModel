package util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

public abstract class Util_Select_Sim_By_Residue {

	protected File baseDir;
	protected Pattern directory_pattern;

	public Util_Select_Sim_By_Residue(File baseDir, Pattern directory_pattern) {
		this.baseDir = baseDir;
		this.directory_pattern = directory_pattern;
	}

	public abstract HashMap<String, double[]> generateResidueMapping(); // Key= DirectoryName:SeedList_Number Val=
																		// value of interest (or residue)

	public static UnivariateFunction extractedWeightedInterpolateFunction(UnivariateInterpolator polator,
			double[] dist_val) {
		Arrays.sort(dist_val);
		if (Double.isNaN(dist_val[dist_val.length - 1])) {
			return null;
		}

		double[] x_val = new double[dist_val.length];
		for (int x = 1; x < dist_val.length; x++) {
			x_val[x] = x_val[x - 1] + 1.0 / dist_val.length;
		}
		x_val[x_val.length - 1] = 1;

		UnivariateFunction resFunc = polator.interpolate(x_val, dist_val);
		return resFunc;
	}

	public static void printOrderedParamList(File baseDir, ArrayList<String> ordered_keys,
			ArrayList<Double> ordered_residue_val, ArrayList<String> inRangeKeyArray,
			HashMap<String, double[]> residue_mapping, File paramListSummaryFile)
			throws FileNotFoundException, IOException {

		int pt = 0;

		double lastResidue = Double.NaN;
		ArrayList<String> sameResidueParamList = new ArrayList<>();
		HashMap<String, String[]> extract_seedList = new HashMap<>();

		ArrayList<StringBuilder> lines_collections = new ArrayList<>();

		for (String key : ordered_keys) {
			Double residue = ordered_residue_val.get(pt);
			String[] key_sp = key.split(":");
			String dir_name = key_sp[0];
			int seedNum = Integer.parseInt(key_sp[1]);

			String[] extract_seed = extract_seedList.get(dir_name);
			double[] residue_val = residue_mapping.get(key);

			if (extract_seed == null) {
				File seedFile = new File(baseDir, dir_name);
				seedFile = new File(seedFile, "Seed_List.csv");
				extract_seed = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(seedFile);

				extract_seedList.put(dir_name, extract_seed);

				// Check for header length consistency
				if (extract_seed.length > 1) {
					int seedLength = extract_seed[1].split(",").length;
					String[] extract_seed_ext = Arrays.copyOf(extract_seed[0].split(","), seedLength);

					StringBuilder str = null;
					for (int i = 0; i < extract_seed_ext.length; i++) {
						if (str == null) {
							str = new StringBuilder();
						} else {
							str.append(',');
						}
						str.append(extract_seed_ext[i]);
					}

					extract_seed[0] = str.toString();

				}

			}

			if (lines_collections.size() == 0) {
				StringBuilder header = new StringBuilder();
				lines_collections.add(header);
				// Print header
				header.append(extract_seed[0]);
				header.append(',');
				header.append(',');
				header.append("Dir");
				header.append(',');
				header.append("Residue");
				header.append(',');
				header.append("In target range");
				for (int i = 0; i < residue_val.length; i++) {
					header.append(String.format(",Outcome_%d", i));
				}
			}

			boolean addNewRow = !residue.equals(lastResidue);

			if (addNewRow) {
				lastResidue = residue;
				sameResidueParamList.clear();
			}

			int sameResiduePt = Collections.binarySearch(sameResidueParamList, extract_seed[seedNum + 1]);
			if (sameResiduePt < 0) {
				addNewRow |= sameResiduePt < 0; // Same residue but difference parameter value.
				sameResidueParamList.add(~sameResiduePt, extract_seed[seedNum + 1]);
			}

			if (addNewRow) {
				StringBuilder line = new StringBuilder();
				lines_collections.add(line);
				line.append(extract_seed[seedNum + 1]); // Seed starts at line 1
				line.append(',');
				line.append(',');
				line.append(key);
				line.append(',');
				line.append(residue);
				line.append(',');
				line.append(Collections.binarySearch(inRangeKeyArray, key) >= 0);
				for (int i = 0; i < residue_val.length; i++) {
					line.append(',');
					line.append(residue_val[i]);
				}
			}
			pt++;
		}

		// Post processing
		StringBuilder[] lines_collections_arr = lines_collections.toArray(new StringBuilder[0]);
		Arrays.sort(lines_collections_arr, 1, lines_collections_arr.length, new Comparator<StringBuilder>() {
			@Override
			public int compare(StringBuilder o1, StringBuilder o2) {
				int res = 0;
				String[] c1 = o1.toString().split(",");
				String[] c2 = o2.toString().split(",");
				for (int c = 0; c < Math.min(c1.length, c2.length) && res == 0; c++) {
					if (c < 2) {
						res = Long.compare(Long.parseLong(c1[c]), Long.parseLong(c2[c]));
					} else {
						res = Double.compare(Double.parseDouble(c1[c]), Double.parseDouble(c2[c]));
					}
				}
				return res;
			}

		});

		PrintWriter pWri_summary = new PrintWriter(paramListSummaryFile);
		for (StringBuilder builder : lines_collections_arr) {
			pWri_summary.println(builder.toString());
		}
		pWri_summary.close();
	}

	/**
	 * Select best simulation/parameter set based on square sum difference from
	 * targeted residue
	 * 
	 * @param residue_map          Residue mapping, most likely from
	 *                             generateResidueMapping() method
	 * @param residue_target_range Residue target and range - either as {single
	 *                             value}, {single value, lower limit, upper limit}
	 *                             or null (i.e. not used)
	 * @param residue_weight       - weight in calculation residue
	 * @param order_residue_val    - input/output of order_resiude_value
	 * @param order_resiude_key    - input/output of order_resiude_key
	 * @return An ArrayList that shows all resiude_key that falls in range
	 */

	public static ArrayList<String> select_best_sim(HashMap<String, double[]> residue_map,
			double[][] residue_target_range, double[] residue_weight, ArrayList<Double> order_residue_val,
			ArrayList<String> order_resiude_key) {

		ArrayList<String> inRangeKeyArr = new ArrayList<>();

		for (Entry<String, double[]> ent : residue_map.entrySet()) {
			double[] data = ent.getValue();
			String key = ent.getKey();
			double residue = 0;
			boolean inRange = true;

			if (residue_target_range != null && residue_weight != null) {
				for (int i = 0; i < data.length && !Double.isNaN(residue); i++) {
					if (residue_target_range[i] != null && residue_weight[i] != 0) {
						inRange &= residue_target_range[i].length == 1
								|| (residue_target_range[i][1] <= data[i] && data[i] <= residue_target_range[i][2]);
						residue += residue_weight[i] * Math.pow(data[i] - residue_target_range[i][0], 2);
					}
				}

				if (!Double.isNaN(residue)) {
					int pt = Collections.binarySearch(order_residue_val, residue);
					if (pt < 0) {
						pt = ~pt;
					}
					order_residue_val.add(pt, residue);
					order_resiude_key.add(pt, key);
					if (inRange) {
						inRangeKeyArr.add(key);
					}
				}
			} else {
				order_residue_val.add(Double.NaN); // Not used
				order_resiude_key.add(key);
				inRangeKeyArr.add(key);
			}
		}

		Collections.sort(inRangeKeyArr);

		return inRangeKeyArr;
	}

	public static HashMap<String, ArrayList<String[]>> combineResultMapping(File[] resFileList)
			throws IOException, FileNotFoundException {
		return combineResultMapping(resFileList, false);

	}

	public static HashMap<String, ArrayList<String[]>> combineResultMapping(File[] resFileList, boolean incl_all)
			throws IOException, FileNotFoundException {
		File[] files = resFileList;
		Pattern keyPattern = Pattern.compile("\\[Seed_List.csv,(\\d+)\\].*");

		HashMap<String, ArrayList<String[]>> entMap = new HashMap<>();
		for (File f : files) {
			if (f.getName().endsWith("7z")) {
				entMap = util.Util_7Z_CSV_Entry_Extract_Callable.extractedLinesFrom7Zip(f, entMap, keyPattern);
			} else {
				String[] lines = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(f);
				ArrayList<String[]> line_arr = new ArrayList<>();
				for (String s : lines) {
					line_arr.add(s.split(","));
				}

				Matcher m = keyPattern.matcher(f.getName());
				m.find();
				entMap.put(m.group(1), line_arr);
			}
		}

		if (!incl_all) {
			ArrayList<String> non_matched = new ArrayList<>();
			for (String key : entMap.keySet()) {
				if (!Pattern.matches("\\d+", key)) {
					non_matched.add(key);
				}
			}
			for (String key : non_matched) {
				System.err.printf("Warning: Key %s not matched with default format and is removed from reuslt format\n",
						key);
				entMap.remove(key);
			}

		}

		return entMap;
	}

}
