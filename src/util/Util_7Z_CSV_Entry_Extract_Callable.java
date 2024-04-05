package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;

public class Util_7Z_CSV_Entry_Extract_Callable implements Callable<Map<String, long[]>> {

	final private File BaseDir;
	final private String[] CSV_Filenames_Pattern;
	final private String CMapId;
	final private int[][] Entry_Pairs; // {[File_id, colNum}
	final private int[] TimePoint; // Read last line if null

	public Util_7Z_CSV_Entry_Extract_Callable(File baseDir, String cMapId, String[] csv_Filenames_Pattern,
			int[][] entry_Pairs, int[] timePoint) {
		super();
		BaseDir = baseDir;
		CMapId = cMapId;
		CSV_Filenames_Pattern = csv_Filenames_Pattern;
		Entry_Pairs = entry_Pairs;
		TimePoint = timePoint;
	}

	public Util_7Z_CSV_Entry_Extract_Callable(File baseDir, String cMapId, String[] csv_Filenames_Pattern,
			int[][] entry_Pairs) {
		this(baseDir, cMapId, csv_Filenames_Pattern, entry_Pairs, null);
	}

	@Override
	public Map<String, long[]> call() throws Exception {
		final int BUFFER = 2048;
		final byte[] buf = new byte[BUFFER];
		HashMap<String, long[]> res = new HashMap<>(); // key = cMapSeed_SimSeed, value =

		for (int file_id = 0; file_id < CSV_Filenames_Pattern.length; file_id++) {
			Pattern pattern_csv_file = Pattern
					.compile(CSV_Filenames_Pattern[file_id].replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));

			File archiveFile = new File(BaseDir,
					(CSV_Filenames_Pattern[file_id].replaceFirst("_%d", "") + ".7z").replaceAll("%d", CMapId));

			SevenZFile archive7Z = new SevenZFile(archiveFile);
			SevenZArchiveEntry ent;
			while ((ent = archive7Z.getNextEntry()) != null) {
				Matcher m = pattern_csv_file.matcher(ent.getName());
				if (m.matches()) {
					String res_key = String.format("%s_%s", m.group(1), m.group(2));

					if (TimePoint != null) {
						for (int i : TimePoint) {
							String res_key_time = String.format("%s_%d", res_key, i);
							long[] res_ent = res.get(res_key_time);
							if (res_ent == null) {
								res_ent = new long[Entry_Pairs.length];
								Arrays.fill(res_ent, 0);
								res.put(res_key_time, res_ent);
							}
						}

					} else {
						long[] res_ent = res.get(res_key);
						if (res_ent == null) {
							res_ent = new long[Entry_Pairs.length];
							Arrays.fill(res_ent, -1);
							res.put(res_key, res_ent);
						}
					}

					StringBuilder txt_entries = new StringBuilder();
					int count;
					while ((count = archive7Z.read(buf, 0, BUFFER)) != -1) {
						txt_entries.append(new String(Arrays.copyOf(buf, count)));
					}
					String line, lastline = "", headerLine = null;
					BufferedReader lines = new BufferedReader(new StringReader(txt_entries.toString()));

					while ((line = lines.readLine()) != null) {
						if (line.length() > 0) {
							if (headerLine == null) {
								headerLine = line;
							} else {
								if (TimePoint != null) {
									String[] lineEnt = line.split(",");
									int time = Integer.parseInt(lineEnt[0]);
									int key = Arrays.binarySearch(TimePoint, time);
									if (key >= 0) {
										String res_key_time = String.format("%s_%d", res_key, TimePoint[key]);
										long[] res_ent = res.get(res_key_time);
										for (int k = 0; k < Entry_Pairs.length; k++) {
											if (Entry_Pairs[k][0] == file_id) {
												try {
													res_ent[k] = Long.parseLong(lineEnt[Entry_Pairs[k][1]]);
												} catch (ArrayIndexOutOfBoundsException ex) {
													System.err.printf("Error reading line with time = %d in %s : %s\n",
															TimePoint[key], BaseDir.getName(), ent.getName());
													if (res_ent[k] != -1) {
														res_ent[k] = -1;
													}
												}
											}
										}
									}
								}
								lastline = line;
							}
						}
					}
					lines.close();

					if (TimePoint == null) {
						long[] res_ent = res.get(res_key);

						String[] lastlineEnt = lastline.split(",");
						for (int k = 0; k < Entry_Pairs.length; k++) {
							if (Entry_Pairs[k][0] == file_id) {
								try {
									res_ent[k] = Long.parseLong(lastlineEnt[Entry_Pairs[k][1]]);
								} catch (ArrayIndexOutOfBoundsException ex) {
									System.err.printf("Error in last line in %s : %s\n", BaseDir.getName(),
											ent.getName());
									if (res_ent[k] != -1) {
										res_ent[k] = -1;
									}

								}
							}
						}
					}
				}
			}
			archive7Z.close();

		}
		return res;
	}

	public static String[] extracted_lines_from_text(File srcTxt) throws FileNotFoundException, IOException {
		ArrayList<String> lines = new ArrayList<>();
		BufferedReader reader = new BufferedReader(new FileReader(srcTxt));
		String line;
		int line_counter = 0;
		while ((line = reader.readLine()) != null) {
			lines.add(line);
			line_counter++;
		}		
		reader.close();
		String[] line_pbs_arr = lines.toArray(new String[line_counter]);
		return line_pbs_arr;
	}

	public static HashMap<String, ArrayList<String[]>> extractedLinesFrom7Zip(File zipFile) throws IOException {
		return extractedLinesFrom7Zip(zipFile, new HashMap<String, ArrayList<String[]>>());
	}

	public static HashMap<String, ArrayList<String[]>> extractedLinesFrom7Zip(File zipFile,
			HashMap<String, ArrayList<String[]>> zip_ent) throws IOException {
		SevenZFile inputZip = new SevenZFile(zipFile);
		SevenZArchiveEntry inputEnt;

		byte[] buf = new byte[Util_Compare_ClusterModel_Transmission_Output.BUFFER];
		while ((inputEnt = inputZip.getNextEntry()) != null) {
			String file_name = inputEnt.getName();
			StringBuilder str_builder = new StringBuilder();
			String line;
			ArrayList<String[]> lines = new ArrayList<>();
			int count;
			while ((count = inputZip.read(buf, 0, Util_Compare_ClusterModel_Transmission_Output.BUFFER)) != -1) {
				str_builder.append(new String(Arrays.copyOf(buf, count)));
			}
			BufferedReader reader = new BufferedReader(new StringReader(str_builder.toString()));
			while ((line = reader.readLine()) != null) {
				if (line.length() > 0) {
					lines.add(line.split(","));
				}
			}
			zip_ent.put(file_name, lines);
		}
		inputZip.close();
		return zip_ent;
	}

	public static StringBuilder extractCSVStringFromZip(SevenZFile inputZip) throws IOException {
		StringBuilder str = new StringBuilder();
		int count;
		final int BUFFER = 2048;
		byte[] buf = new byte[BUFFER];

		while ((count = inputZip.read(buf, 0, BUFFER)) != -1) {
			str.append(new String(Arrays.copyOf(buf, count)));
		}
		return str;
	}

}
