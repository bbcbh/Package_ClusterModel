package util;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Pattern;

public class Util_Analysis_Compare_Cumulative {

	File baseDir;
	File[] compareDirs;
	String[] compareFilePrefix;
	File dataDumpCSV_dir;
	String[] time_to_match;

	final String dataDump_file_name = "cumul_data_compare_%s_%s.csv";
	final String linesMap_key_format = "%s:%s:%s:%s"; // pattern index, dir_name, entry_name, time_to_match
	final String pWriMap_key_format = "%s_%s"; /// pattern index,time_to_match

	public Util_Analysis_Compare_Cumulative(File baseDir, File[] compareDirs, String[] compareFilePrefix,
			String[] time_to_match, File dataDumpCSV_dir) {
		this.baseDir = baseDir;
		this.compareDirs = compareDirs;
		this.compareFilePrefix = compareFilePrefix;
		this.time_to_match = time_to_match;
		this.dataDumpCSV_dir = dataDumpCSV_dir;
	}

	public void generateAnalysisCSV() throws IOException {

		File srcDir = baseDir;
		String[][] header = new String[compareFilePrefix.length][];

		HashMap<String, String[]> base_entry_map = null;
		HashMap<String, float[]> base_entry_diff = null;

		if (baseDir != null) {
			base_entry_map = generateEntryMapFromDirectory(srcDir, compareFilePrefix, header);
			System.out.printf("Reading results from %s completed. # entries in total = %d.\n", srcDir.getAbsolutePath(),
					base_entry_map.size());
			base_entry_diff = generateEntryDiffMap(base_entry_map);
		}

		HashMap<String, PrintWriter> printWriterMap = new HashMap<>();
		boolean headerPrinted = false;
		dataDumpCSV_dir.mkdirs();
				
		for (int pIndex = 0; pIndex < compareFilePrefix.length; pIndex++) {
			for (int t = 1; t < time_to_match.length; t++) {
				String time = time_to_match[t];
				String pWriMapKey = String.format(pWriMap_key_format, Integer.toString(pIndex), time);
				PrintWriter pWri = new PrintWriter(
						new File(dataDumpCSV_dir, String.format(dataDump_file_name, compareFilePrefix[pIndex], time)));				
				printWriterMap.put(pWriMapKey, pWri);
			}
		}

		// Go through each directory and compare value with entry map
		Arrays.sort(compareDirs);

		for (File compDir : compareDirs) {
			HashMap<String, String[]> comp_entry_map = generateEntryMapFromDirectory(compDir, compareFilePrefix, header);							
			HashMap<String, float[]> comp_diff = generateEntryDiffMap(comp_entry_map);
			HashMap<String, String[]> linked_entry_map = generateLinkedEntryMap(compDir, comp_entry_map);
			
			if(!headerPrinted) {
				for (int pIndex = 0; pIndex < compareFilePrefix.length; pIndex++) {
					for (int t = 1; t < time_to_match.length; t++) {
						String pWriMapKey = String.format(pWriMap_key_format, Integer.toString(pIndex), time_to_match[t]);						
						PrintWriter pWri = printWriterMap.get(pWriMapKey);						
						// Print Header
						String[] headerEnt = header[pIndex];
						pWri.print("Result_Index");
						for (int i = 1; i < headerEnt.length; i++) { // Skip time column
							if (base_entry_map != null) {
								pWri.print(',');
								pWri.print(String.format("Abs_Diff[%s]", headerEnt[i]));
								pWri.print(',');
								pWri.print(String.format("Rel_Change[%s]", headerEnt[i]));
							} else {
								pWri.print(',');
								pWri.print(String.format("%s(%s - %s)", headerEnt[i], time_to_match[t], time_to_match[0]));
							}
						}
						String[] extra_header = generateLinkedEntryHeader();
						for (String extraH : extra_header) {
							pWri.print(',');
							pWri.print(extraH);
						}
						pWri.println();												
					}
				}														
				headerPrinted = true;
			}		

			ArrayList<Entry<String, float[]>> entrtSet_arr = new ArrayList<>(comp_diff.entrySet());

			Collections.sort(entrtSet_arr, new Comparator<Entry<String, float[]>>() {
				@Override
				public int compare(Entry<String, float[]> o1, Entry<String, float[]> o2) {
					return o1.getKey().compareTo(o2.getKey());
				}

			});

			for (Entry<String, float[]> diff_map_ent : entrtSet_arr) {
				String diff_key = diff_map_ent.getKey();
				float[] val = Arrays.copyOf(diff_map_ent.getValue(), diff_map_ent.getValue().length);
				
				

				if (base_entry_diff != null) {
					float[] b_val = base_entry_diff.get(diff_key);
					if (b_val == null || val == null) {
						System.err.printf("Warning: Result for <%s> not found in both set.\n", diff_key);
					} else {
						String[] diff_key_split = diff_key.split(":");
						String pWriMapKey = String.format(pWriMap_key_format, diff_key_split[0],
								diff_key_split[diff_key_split.length - 1]);
						PrintWriter pWri = printWriterMap.get(pWriMapKey);
						pWri.print((diff_key_split[1] + ":" + diff_key_split[2]).replaceAll(",", "_"));
						for (int i = 1; i < val.length; i++) { // Skip time column
							pWri.print(',');
							pWri.print(val[i] - b_val[i]);
							pWri.print(',');
							pWri.print(100f * (val[i] - b_val[i]) / b_val[i]);
						}
						String[] linked_entry = getLinkedEntry(linked_entry_map, diff_key);
						if (linked_entry != null) {
							for (String lE : linked_entry) {
								pWri.print(',');
								pWri.print(lE);
							}
						}
						pWri.println();
					}
				} else {
					String[] diff_key_split = diff_key.split(":");
					String pWriMapKey = String.format(pWriMap_key_format, diff_key_split[0],
							diff_key_split[diff_key_split.length - 1]);
					PrintWriter pWri = printWriterMap.get(pWriMapKey);
					pWri.print((diff_key_split[1] + ":" + diff_key_split[2]).replaceAll(",", "_"));	
					for (int i = 1; i < val.length; i++) { // Skip time column
						pWri.print(',');
						pWri.print(val[i]);						
					}				
					String[] linked_entry = getLinkedEntry(linked_entry_map, diff_key);
					if (linked_entry != null) {
						for (String lE : linked_entry) {
							pWri.print(',');
							pWri.print(lE);
						}
					}
					pWri.println();

				}
			}
			System.out.printf("Compare results with %s completed.\n", compDir.getAbsolutePath());
		}

		for (PrintWriter pWri : printWriterMap.values()) {
			pWri.close();
		}

		System.out.printf("Compare results for all (# directories = %d) completed.\n", compareDirs.length);

	}

	protected String[] getLinkedEntry(HashMap<String, String[]> linked_entry_map, String map_key) {
		return null;
	}

	protected HashMap<String, String[]> generateLinkedEntryMap(File compDir, HashMap<String, String[]> comp_entry_map) {
		return new HashMap<>();
	}

	protected String[] generateLinkedEntryHeader() {
		return new String[0];
	}

	protected HashMap<String, float[]> generateEntryDiffMap(HashMap<String, String[]> entry_map) {
		HashMap<String, float[]> entry_diff = new HashMap<>();

		for (Entry<String, String[]> entrySet : entry_map.entrySet()) {
			String[] key_split = entrySet.getKey().split(":");
			if (key_split[key_split.length - 1].equals(time_to_match[0])) {
				String[] entryT0 = entrySet.getValue();
				for (int t = 1; t < time_to_match.length; t++) {
					float[] diffEnt = new float[entryT0.length];
					diffEnt[0] = Float.parseFloat(time_to_match[t]);
					String diffKey = String.format(linesMap_key_format, key_split[0], key_split[1], key_split[2],
							time_to_match[t]);
					String[] entryT = entry_map.get(diffKey);
					for (int c = 1; c < diffEnt.length; c++) {
						if(entryT != null && entryT0 != null) {
							diffEnt[c] = Float.parseFloat(entryT[c]) - Float.parseFloat(entryT0[c]);
						}else {
							diffEnt[c] = 0;							
						}
					}
					entry_diff.put(diffKey, diffEnt);
				}
			}
		}

		return entry_diff;
	}

	protected HashMap<String, String[]> generateEntryMapFromDirectory(File srcDir, String[] compareFilePrefix,
			String[][] header) throws IOException {
		HashMap<String, String[]> entry_map = new HashMap<>();
		File[] resultDirs = srcDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory();
			}
		});

		for (File dir : resultDirs) {
			for (int pIndex = 0; pIndex < compareFilePrefix.length; pIndex++) {
				Pattern p = Pattern.compile(String.format("%s_.*.csv.7z", compareFilePrefix[pIndex]));
				File[] zipFiles = dir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return p.matcher(pathname.getName()).matches();
					}
				});
				for (File zf : zipFiles) {
					HashMap<String, ArrayList<String[]>> entry = util.Util_7Z_CSV_Entry_Extract_Callable
							.extractedLinesFrom7Zip(zf);
					for (String key : entry.keySet()) {
						ArrayList<String[]> lines = entry.get(key);
						for (String[] line : lines) {
							if (header != null && header[pIndex] == null) {
								header[pIndex] = line;
							}

							for (String time : time_to_match) {
								if (time.equals(line[0])) {
									String linemap_key = String.format(linesMap_key_format, Integer.toString(pIndex),
											dir.getParentFile().getName() +"->"+ dir.getName() , key, time);
									entry_map.put(linemap_key, line);
								}
							}
						}
					}
				}
			}

			// Raw file version
			Pattern[] compareFile_CSV_Pattern = new Pattern[compareFilePrefix.length];

			for (int pIndex = 0; pIndex < compareFile_CSV_Pattern.length; pIndex++) {
				Pattern p = Pattern.compile(String.format("\\[.*\\]%s_-?\\d+_.*.csv", compareFilePrefix[pIndex]));
				File[] csvFiles = dir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return p.matcher(pathname.getName()).matches();
					}
				});
				for (File csv : csvFiles) {
					String[] lines = util.Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(csv);
					for (String lineTotal : lines) {
						String[] line = lineTotal.split(",");
						if (header != null && header[pIndex] == null) {
							header[pIndex] = line;
						}
						for (String time : time_to_match) {
							if (time.equals(line[0])) {
								String linemap_key = String.format(linesMap_key_format, Integer.toString(pIndex),
									dir.getParentFile().getName() +"->"+ dir.getName() , csv.getName(), time);
								entry_map.put(linemap_key, line);
							}
						}
					}
				}
			}
		}
		return entry_map;
	}

}
