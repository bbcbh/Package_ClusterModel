package util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public abstract class Util_Select_Sim_By_Residue {

	protected File baseDir;
	protected Pattern directory_pattern;

	public Util_Select_Sim_By_Residue(File baseDir, Pattern directory_pattern) {
		this.baseDir = baseDir;
		this.directory_pattern = directory_pattern;
	}

	public abstract HashMap<String, double[]> generateResidueMapping(); // Key= DirectoryName:SeedList_Number Val=
																		// value of interest (or residue)
	public static void printOrderedParamList(File baseDir,
			ArrayList<String> ordered_keys,	ArrayList<Double> ordered_residue_val,
			ArrayList<String> inRangeKeyArray,
			HashMap<String, double[]> residue_mapping, File paramListSummaryFile)
			throws FileNotFoundException, IOException {
		PrintWriter pWri_summary = null;
		HashMap<String, String[]> extract_seedList = new HashMap<>();
		int pt = 0;
		for (String key : ordered_keys) {
			Double residue = ordered_residue_val.get(pt);
			String[] key_sp = key.split(":");
			String dir_name = key_sp[0];
			int seedNum = Integer.parseInt(key_sp[1]);
			if(seedNum == 0) {
				int k = 1;
			}
				

			String[] extract_seed = extract_seedList.get(dir_name);
			double[] residue_val = residue_mapping.get(key);

			if (extract_seed == null) {
				File seedFile = new File(baseDir, dir_name);
				seedFile = new File(seedFile, "Seed_List.csv");
				extract_seed = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(seedFile);
				extract_seedList.put(dir_name, extract_seed);
			}

			if (pWri_summary == null) {
				pWri_summary = new PrintWriter(paramListSummaryFile);
				// Print header
				pWri_summary.print(extract_seed[0]);
				pWri_summary.print(',');
				pWri_summary.print(',');
				pWri_summary.print("Dir");
				pWri_summary.print(',');
				pWri_summary.print("Residue");
				pWri_summary.print(',');
				pWri_summary.print("In target range");
				for (int i = 0; i < residue_val.length; i++) {
					pWri_summary.print(String.format(",Outputcome_%d", i));
				}
				pWri_summary.println();
			}

			pWri_summary.print(extract_seed[seedNum+1]); // Seed starts at line 1
			pWri_summary.print(',');
			pWri_summary.print(',');
			pWri_summary.print(key);
			pWri_summary.print(',');
			pWri_summary.print(residue);
			pWri_summary.print(',');
			pWri_summary.print(Collections.binarySearch(inRangeKeyArray, key) >=0);
			for (int i = 0; i < residue_val.length; i++) {
				pWri_summary.print(',');
				pWri_summary.print(residue_val[i]);
			}
			pWri_summary.println();
			pt++;
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
		}		
		
		Collections.sort(inRangeKeyArr);
					
		return inRangeKeyArr;
	}

	public static HashMap<String, ArrayList<String[]>> combineResultMapping(File[] resFileList)
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
		return entMap;
	}

}
