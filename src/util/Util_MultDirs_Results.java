package util;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Util_MultDirs_Results {

	public static int combineResultDirectories(File combinedDir, File inputDir) {
		File[] dir_to_move = inputDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory();
			}
		});

		int sucCounters = 0;		

		for (File dir : dir_to_move) {
			String dirName = dir.getName();
			File targetPath = new File(combinedDir, dirName);

			// Regenerate new target path until it is new
			int counter = 1;
			String org_dir = dirName;
			
			Matcher m = Pattern.compile("\\d+").matcher(dirName);
			m.find();
			int str_length = m.end() - m.start();
			
			while (targetPath.exists()) {													
				dirName = dirName.replaceFirst("\\d+", String.format("%0" + str_length + "d", counter));				
				targetPath = new File(combinedDir, dirName);
				counter++;
			}			
			if (org_dir.equals(dirName)) { // No number
				dirName = String.format("%s_%s", dirName, String.format("%02d", counter));
			}

			if (!targetPath.exists()) {
				try {
					Files.move(dir.toPath(), targetPath.toPath(), StandardCopyOption.ATOMIC_MOVE);
					sucCounters++;
				} catch (Exception e) {
					e.printStackTrace(System.err);
				}
			} else {
				System.err.printf("combineResultDirectories: Cannot found a new name for %s. Diretory NOT moved.\n",
						dir.getAbsolutePath());

			}

		}
		
		return sucCounters;
		

	}

	public static String[] generatePreDefinedSeedOrderList(File result_dir_base, Pattern pattern_dir_select,
			File param_list, int num_col) throws IOException {
		String[] lines = util.Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(param_list);

		for (int i = 0; i < lines.length; i++) {
			String[] lineSp = lines[i].split(",");
			StringBuilder strBd = new StringBuilder(lineSp[0]);
			for (int c = 1; c < num_col; c++) {
				strBd.append(',');
				strBd.append(lineSp[c]);
			}
			lines[i] = strBd.toString();
		}

		HashMap<String, String> entryMap = new HashMap<>(); // Key = lines, Val = DirName_SeedList:EntNum

		File[] result_dirs = result_dir_base.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && pattern_dir_select.matcher(pathname.getName()).matches();
			}
		});

		for (File resDir : result_dirs) {
			String[] seedLine = util.Util_7Z_CSV_Entry_Extract_Callable
					.extracted_lines_from_text(new File(resDir, "Seed_List.csv"));

			for (int i = 1; i < seedLine.length; i++) {
				entryMap.put(seedLine[i], String.format("%s:%d", resDir.getName(), i - 1));
			}
		}

		String[] key = new String[lines.length]; // DirName:SeedList_EntNum

		for (int i = 1; i < key.length; i++) {
			key[i] = entryMap.get(lines[i]);

			if (key[i] == null) {
				System.err.printf("generatePreDefinedSeedOrderList: Warning. The entry correspond to <%s> not found!\n",
						lines[i]);
			}
		}

		return key;

	}

	public static void printExtractedResultFromMultiDirs(File result_dir_base, Pattern pattern_dir_select,
			String[] compareFilePrefix, int[][] col_sel_all) throws IOException, FileNotFoundException {
		printExtractedResultFromMultiDirs(result_dir_base, pattern_dir_select, compareFilePrefix, col_sel_all, null);
	}

	public static void printExtractedResultFromMultiDirs(File result_dir_base, Pattern pattern_dir_select,
			String[] compareFilePrefix, int[][] col_sel_all, String[] param_list_order)
			throws IOException, FileNotFoundException {
		HashMap<String, PrintWriter> pWri_collection = new HashMap<>(); // Key: file_prefex_col_number
		File extract_values_base = new File(result_dir_base, "Extract_Values");
		extract_values_base.mkdirs();

		File[] result_dirs = result_dir_base.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && pattern_dir_select.matcher(pathname.getName()).matches();
			}
		});

		for (int pId = 0; pId < compareFilePrefix.length; pId++) {			
			final String prefix = compareFilePrefix[pId];
			
			HashMap<String, HashMap<String, ArrayList<String[]>>> resMapByDirectory = extractedResultsFromResultDirs(
					result_dirs, prefix);

			String[] order_key_list = param_list_order;

			if (order_key_list == null) {
				ArrayList<String> order_key_arr = new ArrayList<>();
				for (Entry<String, HashMap<String, ArrayList<String[]>>> entry : resMapByDirectory.entrySet()) {
					for (String zipKey : entry.getValue().keySet()) {
						order_key_arr.add(String.format("%s:%s", entry.getKey(), zipKey));
					}
				}
				order_key_list = order_key_arr.toArray(new String[0]);
			}
			for (String key : order_key_list) {
				if (key != null) {
					String[] key_arr = key.split(":");
					HashMap<String, ArrayList<String[]>> resMap = resMapByDirectory.get(key_arr[0]);
					ArrayList<String[]> lines = resMap.get(key_arr[1]);
					if (lines == null) {
						System.err.printf("extractedResultFromMultiDirs: Result for <%s> not found.\n", key);
					} else {
						for (int col : col_sel_all[pId]) {
							String pWri_k = String.format("%s_%d", prefix, col);
							PrintWriter pWri = pWri_collection.get(pWri_k);
							if (pWri == null) {
								pWri = new PrintWriter(new File(extract_values_base, String.format("%s.csv", pWri_k)));
								pWri_collection.put(pWri_k, pWri);
								pWri.print("Entry_Name");
								// Time line
								for (String[] ent : lines) {
									pWri.print(',');
									pWri.print(ent[0]);
								}
								pWri.println();
							}
							// Entry name
							pWri.print(key);
							for (String[] ent : lines) {
								pWri.print(',');
								pWri.print(ent[col]);

							}
							pWri.println();
						}
					}
				}
			}
		}

		for (PrintWriter pWri : pWri_collection.values()) {
			pWri.close();
		}
	}
	
	
	public static HashMap<String, HashMap<String, ArrayList<String[]>>> extractedResultsFromSimDirBase(File simDirBase, final String result_file_prefix) 
			throws IOException {			
		File[] simDirs = simDirBase.listFiles(new FileFilter() {			
			@Override
			public boolean accept(File pathname) {				
				return pathname.isDirectory();
			}
		});				
		return extractedResultsFromResultDirs(simDirs,result_file_prefix);
	}	
	
	private static HashMap<String, HashMap<String, ArrayList<String[]>>> extractedResultsFromResultDirs(File[] result_dirs,
			final String result_file_prefix) throws IOException, FileNotFoundException {
		
		HashMap<String, HashMap<String, ArrayList<String[]>>> resMapByDirectory = new HashMap<>();
		for (File result_dir : result_dirs) {

			File[] res_file = result_dir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.getName().startsWith(result_file_prefix)
							|| Pattern.matches("\\[.*\\]" + result_file_prefix + "_.*.csv", pathname.getName());
				}
			});
			HashMap<String, ArrayList<String[]>> resMap = new HashMap<>();
			Pattern keyPattern = Pattern.compile("\\[Seed_List.csv,(\\d+)\\].*");

			for (File f : res_file) {
				if (f.getName().endsWith("7z")) {
					resMap = util.Util_7Z_CSV_Entry_Extract_Callable.extractedLinesFrom7Zip(f, resMap, keyPattern);
				} else {
					String[] lines = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(f);
					ArrayList<String[]> line_arr = new ArrayList<>();
					for (String s : lines) {
						line_arr.add(s.split(","));
					}

					Matcher m = keyPattern.matcher(f.getName());
					m.find();
					resMap.put(m.group(1), line_arr);
				}
			}
			resMapByDirectory.put(result_dir.getName(), resMap);
		}
		return resMapByDirectory;
	}

}
